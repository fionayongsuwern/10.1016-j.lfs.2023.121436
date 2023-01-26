%% Data Preparation
% Make sure mean data has headers in the format 'Mean1', 'Mean2'...'Mean30', etc
% Make sure coordinate data has headers 'X' for column with x-coordinates and 'Y' for column with y-coordinates
close all
clear
clc


%% Inputs
file = 'exp1-i4';
ext = '.xlsx';
filename = horzcat(file,ext);
coordfile = horzcat('coords_',file,ext);
leadersfile = horzcat('wave originator_',file,ext);


%% Import data
data = readtable(filename);
coordata = readtable(coordfile);
leadersdata = readtable(leadersfile);
numOfCols = width(data);
numOfRows = height(data);


%% Denoising
newData = table;

for var = 1:1:numOfCols
    Mean = horzcat('Mean',num2str(var));
    newMean = zeros(numOfRows,1);
    
    % Huang-Hilbert empirical mode decomposition
    IMF = emd(data.(Mean));
    IMF = IMF';
    [HS,f,t,imfinsf,imfinse] = hht(IMF,1);
    t = t+1;
    x = IMF.*t.*(exp(1i*(cumtrapz(imfinsf)))); 
    x = x(:,2:end-1);
    
    % Setting 20% threshold
    Baseline = 1.2*real(sum(x,2));  
    baselineData = data.(Mean)-Baseline;
    
    % Binarising signal
    for row = 1:1:numOfRows
        if baselineData(row) > 0
            newMean(row) = 1;
        elseif baselineData(row) < 0
            newMean(row) = 0;
        end
    end
    
    newData = [newData array2table(newMean,'VariableNames',{(Mean)})];
end

% Data set after HHT-EMD
newData = [array2table(t, 'VariableNames',{'Time(s)'}) newData];


%% Coactivity for 3D Data

for S = 1:1:numOfCols
    Coactivity = table;
    S
        
    for T1 = 1:1:numOfCols
        Ti = horzcat('Mean',num2str(T1));
        C = horzcat('Mean',num2str(T1));
        CoAct.(C) = table;
        
        for T2 = 1:1:numOfCols
            Tj = horzcat('Mean',num2str(T2));
            
            Cij = (sum((newData.(Ti) == 1) & (newData.(Tj) == 1)))/(sqrt(sum(newData.(Ti)) * sum(newData.(Tj))));
            CoAct.(C) = [CoAct.(C) array2table(Cij,'VariableNames',{Tj})];
        end
        
        Coactivity = [Coactivity; CoAct.(C)];
    end
end

save(file);


%% Monte Carlo Simulation for 3D Data
n = 10000; % number of iterations
for j = 1:1:numOfCols
    M = horzcat('Mean',num2str(j));
    MonteCarloCoactivity.(M) = table;
end


for i = 7001:1:n
    MCS = table;
            
    % Generate Dataset
    for row = 1:1:numOfCols
        yName = horzcat('Mean',num2str(row));
        y = randi([0,1],numOfRows,1);
        MCS = [MCS array2table(y,'VariableNames',{yName})];
    end
     
    % Generate Monte Carlo Coactivity
    for T1 = 1:1:numOfCols
        Ti = horzcat('Mean',num2str(T1));
        Set.(Ti) = table;
        
        for T2 = 1:1:numOfCols
            Tj = horzcat('Mean',num2str(T2));
            
            Cij = (sum((MCS.(Ti) == 1) & (MCS.(Tj) == 1)))/(sqrt(sum(MCS.(Ti)) * sum(MCS.(Tj))));
            Set.(Ti) = [Set.(Ti) array2table(Cij,'VariableNames',{Tj})];
        end
    end
    
    for j = 1:1:numOfCols
        M = horzcat('Mean',num2str(j));
        
        MonteCarloCoactivity.(M) = [MonteCarloCoactivity.(M); Set.(M)];
    end
    i
    
    if rem(i,100) == 0
        save(file);
    end
end


%% Comparison to Monte Carlo 3D dataset
for S = 1:1:numOfCols
    Significance = table;
    S
        
    for i = 1:1:numOfCols
        A = horzcat('Mean',num2str(i));
        Sign.(A) = table;
        
        for j = 1:1:numOfCols
            BO = horzcat('Mean',num2str(j));
            BMC = horzcat('Mean',num2str(j));
            
            Observed = CoAct.(A).(BO);
            MonteCarlo = MonteCarloCoactivity.(A).(BMC);
            Sg = ttest2(Observed,MonteCarlo,'Alpha',0.001);
            
            Sign.(A) = [Sign.(A) array2table(Sg,'VariableNames',{BMC})];
        end
        
        Significance = [Significance; Sign.(A)];
    end
end

% Significant Coactivity
numOfRows = height(Significance);
SignificantCoactivity = Significance;

for i = 1:1:numOfCols
    A = horzcat('Mean',num2str(i));
        
    for j = 1:1:numOfRows
        if Significance.(A)(j) == 1
           SignificantCoactivity.(A)(j) = Coactivity.(A)(j);
        elseif Significance.(A)(j) == 0
           SignificantCoactivity.(A)(j) = 0;
        elseif isnan(Significance.(A)(j))
           SignificantCoactivity.(A)(j) = NaN;
        end
    end
end

save(file);


%% Checking for Power-Law Relationship
threshold = 0.8; % Setting threshold at minimum % coactivity
threshold2 = 0.6;
threshold3 = 0.4;
connections = 30; % Setting connectivity % threshold for number of hubs

Percent = sum(table2array(SignificantCoactivity)>=threshold)/(numOfCols-1)*100;
Percent2 = sum(table2array(SignificantCoactivity)>=threshold2)/(numOfCols-1)*100;
Percent3 = sum(table2array(SignificantCoactivity)>=threshold3)/(numOfCols-1)*100;
Hub = Percent>=connections;
Hub2 = Percent2>=connections;
Hub3 = Percent3>=connections;
Hubs = num2cell(Hub);
Hubs = ['Hub Cells',Hubs];
numOfHubs = sum(Percent>=connections);
numOfHubs = numOfHubs/numOfCols*100;
numOfHubs = horzcat('Percentage Hub Cells = ',num2str(numOfHubs),'%');
PercentConnect = num2cell(Percent);
PercentConnect = ['Percent of Connections',PercentConnect];

X = unique(Percent);
Y = zeros(1,length(X));

for i = 1:1:length(X)
    Y(i) = (sum(Percent==X(i)))/numOfCols*100;
end

PowerLaw = figure('Name','Power Law','NumberTitle','off');

try
    [slope, intercept, MSE, R2, S] = logfit(X,Y,'loglog');
    R = horzcat('R² = ',num2str(R2));
    Threshold = horzcat('Coactivity Threshold = ',num2str(threshold));
    MinCons = horzcat('Min Percentage Connections = ',num2str(connections),'%');
    txt = [R;cellstr(numOfHubs);cellstr(Threshold);cellstr(MinCons)];
catch
    txt = ('Not enough data points to plot a graph.');
    txt = txt';
end

xlabel('Number of Connections','FontSize',12,'Color','k','FontName','Arial');
ylabel('Probability (%)','FontSize',12,'Color','k','FontName','Arial');
yticks([1 10 100]);
yticklabels({'1' '10' '100'});
ax.YAxis.Exponent = 0;
xticks([0 100]);
xticklabels({'1' '100'});
ax.XAxis.Exponent =0;
set(gca,'TickDir','out','box','off','XMinorTick','off','YMinorTick','off');

PowerLawWidth = PowerLaw.Position(3);
PowerLawWidth = PowerLawWidth*1.55;

PowerLaw.Position(3) = PowerLawWidth;

ax = gca;
ax.Position(3) = 0.5;
str = txt';
annotation('textbox',[0.65,0.1,0.3,0.4],'String',str,'FitBoxToText','on','FontSize',10,'Color','k','FontWeight','bold');

file1 = horzcat(file,'_Power Law');
saveas(gcf,file1,'jpg');
saveas(gcf,file1,'bmp');


%% Plot Calcium Traces
CalciumTraces = figure('Name','Calcium Traces','NumberTitle','off');

for var = 1:1:numOfCols
    Mean = horzcat('Mean',num2str(var));
    
    fluorescenceY = data.(Mean);
    plot(fluorescenceY);
    
    hold on
end

ax = gca;
ax.XAxis.Label.String = 'Time';
ax.XAxis.Label.FontSize = 14;
ax.XAxis.Label.FontWeight = 'bold';
ax.YAxis.Label.String = 'Fluorescence';
ax.YAxis.Label.FontWeight = 'bold';
ax.YAxis.Label.FontSize = 14;
ax.YAxis.Exponent = 0;
ax.YGrid = 'on';

file2 = horzcat(file,'_Calcium Traces');
saveas(gcf,file2,'jpg');
saveas(gcf,file2,'bmp');


%% Plot Topographic Connectivity

% General
XY = tail(coordata,numOfCols);
XY = [XY.X XY.Y];
xpos = XY(:,1);
xpos = xpos+2;
ypos = XY(:,2);
cellLabel = coordata(:,1);
cellLabel = head(cellLabel,numOfCols);
cellLabel = table2cell(cellLabel);
CoactivityPlot = table2array(SignificantCoactivity);
TopographicConnectivity = figure('Name','Connectivity','NumberTitle','off');

try
    wgPlot(CoactivityPlot,XY,'edgeColorMap',jet,'edgeWidth',0.5,'vertexMarker','.');
end

hold on

Hub3 = [XY Hub3'];
Hub3(Hub3(:,3)==0,:) = [];
Hub3(:,3) = [];

Hub3X = Hub3(:,1);
Hub3Y = Hub3(:,2);

Hub3Plot = scatter(Hub3X,Hub3Y,85,'MarkerFaceColor','#FFFFFF','MarkerEdgeColor','k');

hold on

Hub2 = [XY Hub2'];
Hub2(Hub2(:,3)==0,:) = [];
Hub2(:,3) = [];

Hub2X = Hub2(:,1);
Hub2Y = Hub2(:,2);

Hub2Plot = scatter(Hub2X,Hub2Y,150,'MarkerFaceColor','#A6A6A6','MarkerEdgeColor','k');

hold on

Hub = [XY Hub'];
Hub(Hub(:,3)==0,:) = [];
Hub(:,3) = [];

HubX = Hub(:,1);
HubY = Hub(:,2);

HubPlot = scatter(HubX,HubY,225,'MarkerFaceColor','#000000','MarkerEdgeColor','k');

hold on

WaveOrigin = zeros(1,numOfCols);

if height(waveorigindata)>0
    WOX = zeros(width(waveorigindata),1);
    WOY = zeros(width(waveorigindata),1);
    
    for i = 1:1:width(waveorigindata)
        LR = horzcat('GreenCell',num2str(i));

        WOX(i,1) = XY(waveorigindata.(LR),1);
        WOY(i,1) = XY(waveorigindata.(LR),2);

        WaveOrigin(1,waveorigindata.(LR)) = 1;
    end
    
%     WaveOriginPlot = scatter(WOX,WOY,500,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',2);
end

WaveOrigin = num2cell(WaveOrigin);
WaveOrigin = ['Wave Origins',WaveOrigin];

file3 = horzcat(file,'_Connectivity Map');
saveas(gcf,file3,'jpg');
saveas(gcf,file3,'bmp');

text(xpos,ypos,cellLabel);

file4 = horzcat(file3,' (with labels)');
saveas(gcf,file4,'bmp');


%% Save workspace and Export Data to Excel
newFilename = horzcat('Coactivity_',filename);
save(file);

% Original Raw Data
data = [array2table(t, 'VariableNames',{'Time(s)'}) data];
writetable(data,newFilename,'Sheet',filename);

% Binarised Signal
writetable(newData,newFilename,'Sheet','Binarised Signal');

% Binarised Coactivity
Means = (data.Properties.VariableNames)';
Means(1,:) = [];

Coactivity = [array2table(Means) Coactivity];
writetable(Coactivity,newFilename,'Sheet','Binarised Coactivity');

% Significant Coactivity
SignificantCoactivity = [array2table(Means) SignificantCoactivity];
SignificantCoactivity = [SignificantCoactivity;Hubs;PercentConnect;Nnat];
writetable(SignificantCoactivity,newFilename,'Sheet','Significant Coactivity');% Calcium Traces
xlswritefig(CalciumTraces,newFilename,filename,'B2');

% Topographic Connectivity
xlswritefig(TopographicConnectivity,newFilename,'Significant Coactivity','C5');

% Power Law
xlswritefig(PowerLaw,newFilename,'Significant Coactivity','K5');
