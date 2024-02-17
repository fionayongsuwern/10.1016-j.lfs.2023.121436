%% Data Preparation
% Make sure mean data has headers in the format 'Mean1', 'Mean2'...'Mean30', etc
% Make sure coordinate data has headers 'X' for column with x-coordinates and 'Y' for column with y-coordinates
close all
clear
clc


%% Inputs
rootdir = 'C:\Users\fiona\Desktop\School\Nanyang Technological University\Dissertation\Data\Nnat';
cd(rootdir);

batch = 'Additional';
exp = 3;
ani = 'WT3';
ramp = 0; %either 3G to 11G or just 11G
isl = 5;
cond = 11; %glucose concentration

exp = horzcat('Exp ',num2str(exp));
file = horzcat('Islet ',num2str(isl));
ext = '.xlsx';

if ramp == 0
    ramp = '11G';
    
    cd(fullfile(rootdir,batch,exp,ani,ramp,file));
    
    cond = horzcat(num2str(cond),'G');
    filename = horzcat(file,' ',cond,ext);
    coordfile = horzcat(file,' ',cond,' coords',ext);
    nnatfile = horzcat(file,' ',cond,' nnat cells',ext);
    
elseif ramp == 1
    ramp = '3G to 11G';
    
    cd(fullfile(rootdir,batch,exp,ani,ramp,file));
    
    if cond == 3
        cond = horzcat(num2str(cond),'G');
        file = horzcat(file,' ',cond);
        filename = horzcat(file,ext);
        coordfile = horzcat(file,' coords',ext);
        nnatfile = horzcat(file,' nnat cells',ext);
    elseif cond == 11
        cond = horzcat(num2str(cond),'G');
        file = horzcat(file,' ','3G_',cond);
        filename = horzcat(file,ext);
        coordfile = horzcat(file,' coords',ext);
        nnatfile = horzcat(file,' nnat cells',ext);
    end
end


% Import data
data = readtable(filename);
coordata = readtable(coordfile);
% nnatdata = readtable(nnatfile);
numOfCols = width(data);
numOfRows = height(data);


%% Denoising (smoothing spline)
newData = table;
t = 1:1:numOfRows;
t = t';

SmoothCalciumTraces = figure('Name','Calcium Traces (Smoothed)','NumberTitle','off');
BinarisedCalciumTraces = figure('Name','Calcium Traces (Binarised)','NumberTitle','off');

for i = 1:1:numOfCols
    Mean = horzcat('Mean',num2str(i));
    i
    
    fluorescenceY = data.(Mean);
    fluorescenceYii = fluorescenceY;
    
    fluorescenceYi = fit(t,fluorescenceY,'smoothingspline','SmoothingParam',0.05);
    
    figure(SmoothCalciumTraces);
    
    plot(fluorescenceYi);
    
    legend('off')
    
    baselineMean = feval(fluorescenceYi,1:numOfRows);
    baselineMean = mean(baselineMean);

    for row = 1:1:numOfRows
        if fluorescenceYi(row) > baselineMean
            fluorescenceYii(row) = 1;
        elseif fluorescenceYi(row) <= baselineMean
            fluorescenceYii(row) = 0;
        end
    end
    
    newData = [newData array2table(fluorescenceYii,'VariableNames',{(Mean)})];
    
    hold on
    
    figure(BinarisedCalciumTraces);

    plot(fluorescenceYii);
    
    hold on
end

figure(SmoothCalciumTraces);

figure(BinarisedCalciumTraces);

figure(SmoothCalciumTraces);
ax = gca;
ax.XAxis.Label.String = 'Time';
ax.XAxis.Label.FontSize = 14;
ax.XAxis.Label.FontWeight = 'bold';
ax.YAxis.Label.String = 'Fluorescence';
ax.YAxis.Label.FontWeight = 'bold';
ax.YAxis.Label.FontSize = 14;
ax.YAxis.Exponent = 0;

SmoothCalciumTraces = gcf;

figure(BinarisedCalciumTraces);
ax = gca;
ax.XAxis.Label.String = 'Time';
ax.XAxis.Label.FontSize = 14;
ax.XAxis.Label.FontWeight = 'bold';
ax.YAxis.Label.String = 'Binarised Signal';
ax.YAxis.Label.FontWeight = 'bold';
ax.YAxis.Label.FontSize = 14;
ax.YAxis.Exponent = 0;
ylim([-0.1 1.2])

ExampleCalciumTraces = figure('Name','Calcium Traces','NumberTitle','off');
yyaxis left
plot(fluorescenceY);

hold on

yyaxis left
plot(fluorescenceYi);

hold on

yyaxis right
plot(fluorescenceYii);

legend('off')
ax = gca;
ax.XAxis.Label.String = 'Time';
ax.XAxis.Label.FontSize = 14;
ax.XAxis.Label.FontWeight = 'bold';

yyaxis left
ax.YAxis(1).Label.String = 'Fluorescence';
ax.YAxis(1).Label.FontWeight = 'bold';
ax.YAxis(1).Label.FontSize = 14;
ax.YAxis(1).Exponent = 0;

yyaxis right
ax.YAxis(2).Label.String = 'Binarised Signal';
ax.YAxis(2).Label.FontWeight = 'bold';
ax.YAxis(2).Label.FontSize = 14;
ax.YAxis(2).Limits = [-0.1 1.2];

file1 = horzcat(file,'_Example Processed Calcium Traces');
saveas(gcf,file1,'jpg');
saveas(gcf,file1,'bmp');

% Data set after smoothing
newData = [array2table(t, 'VariableNames',{'Time(s)'}) newData];

close all


% Coactivity

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

save(file,'-v7.3');


% Data shuffling
n = 10000; % number of iterations
for j = 1:1:numOfCols
    M = horzcat('Mean',num2str(j));
    MonteCarloCoactivity.(M) = table;
end

for i = 1:1:n
    MCS = table;
            
    % Generate Dataset
    for column = 1:1:numOfCols
        yName = horzcat('Mean',num2str(column));
        y = newData.(yName);
        random_y = y(randperm(size(y,1)),:);
        MCS = [MCS array2table(random_y,'VariableNames',{yName})];
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
end

save(file,'-v7.3');
disp('data shuffled');


% Comparison to Monte Carlo dataset
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
    A
        
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


SignificantCoactivity = Coactivity;
for i = 1:1:numOfCols
    A = horzcat('Mean',num2str(i));
    i
        
    for j = 1:1:numOfCols
        if SignificantCoactivity.(A)(j) == 1
           SignificantCoactivity.(A)(j) = 0;
        end
    end
end

save(file,'-v7.3');
disp('monte carlo significant coactivity');


% Checking for Power-Law Relationship
threshold = 0.8; % Setting threshold at minimum % coactivity (decimal)
threshold2 = 0.6;
threshold3 = 0.4;
connections = 30; % Setting connectivity % threshold for number of hubs (percent)

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
ylabel('Probability','FontSize',12,'Color','k','FontName','Arial');
yticks([1 10 100]);
yticklabels({'1' '10' '100'});
ax = gca;
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

file2 = horzcat(file,'_Power Law');
saveas(gcf,file2,'jpg');
saveas(gcf,file2,'bmp');


% Plot Calcium Traces
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

file3 = horzcat(file,'_Calcium Traces');
saveas(gcf,file3,'jpg');
saveas(gcf,file3,'bmp');


% Plot Topographic Connectivity

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

% Nnat = zeros(1,numOfCols);
% 
% if height(nnatdata)>0
%     nnatX = zeros(width(nnatdata),1);
%     nnatY = zeros(width(nnatdata),1);
%     
%     for i = 1:1:width(nnatdata)
%         LR = horzcat('LR',num2str(i));
% 
%         nnatX(i,1) = XY(nnatdata.(LR),1);
%         nnatY(i,1) = XY(nnatdata.(LR),2);
% 
%         Nnat(1,nnatdata.(LR)) = 1;
%     end
%     
%     NnatPlot = scatter(nnatX,nnatY,500,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',2);
% end
% 
% Nnat = num2cell(Nnat);
% Nnat = ['Nnat Cells',Nnat];

file4 = horzcat(file,'_Connectivity Map');
saveas(gcf,file4,'jpg');
saveas(gcf,file4,'bmp');

text(xpos,ypos,cellLabel);

file5 = horzcat(file4,' (with labels)');
saveas(gcf,file5,'bmp');


% Pearson's R
N = readmatrix(filename);
% N = N(51:end,:);

% Read sample length, this gives you number of columns in analysis
SampleLength=size(N,2);
numofrows = length(N);

% Correlations (with 5% of sample length smoothing to reduce noise)
for i=1:SampleLength
    for j=1:SampleLength
        X=[];
        Y=[]; 

        X=smooth(N([1:numofrows],i),0.05); %added 5% smooth  
        Y=smooth(N([1:numofrows],j),0.05); %added 5% smooth
        
        [Correlation,Significance]=corrcoef(X,Y);
        
        test(i,j)=Correlation(2,1);
        sign(i,j)=Significance(2,1);
    end
end

% Plot heatmap
PearsonRFig = figure;
imagesc(test(:,:),[0 1]) %these are the correlation values
%drawnow
% title('Heatmap of correlation matrix')
% title({['Heatmap of correlation matrix for ',num2str(file)]},'FontSize',12,'Color','k','FontName','Arial')
xlabel('Cell number','VerticalAlignment','bottom','FontSize',12,'Color','k','FontName','Arial')
xticks([1 numOfCols]);
ylabel('Cell number','VerticalAlignment','top','FontSize',12,'Color','k','FontName','Arial')
yticks([1 numOfCols]);
set(gca,'box','off','XMinorTick','off','YMinorTick','off');

c = colorbar;
c.Ticks = [0 1];
c.Label.String = 'r';
c.Label.VerticalAlignment = 'bottom';
c.Label.FontSize = 12;
c.Label.FontName = 'Arial';

file6 = horzcat(file,'_Correlation Heat Map');
saveas(gcf,file6,'jpg');
saveas(gcf,file6,'bmp');

% This function excludes autocorrelation values (1),and calculates the significance of each corr coeff value.
bb=(mean(nonzeros(triu(test,1))));
cc=mean(nonzeros(triu(sign,1)));
Correlations=test(:,:);
Signif=sign(:,:);

SigPearsonRFig = figure;
imagesc(sign(:,:),[0 0.1]) 
drawnow
colorbar
title('Heatmap of significance matrix')
xlabel('Cell number')
ylabel('Cell number')
title({['Heatmap of significance matrix for ',num2str(file)]});

file7 = horzcat(file,'_Significant Correlation Heat Map');
saveas(gcf,file7,'jpg');
saveas(gcf,file7,'bmp');

AvgPearsonR = ["Pearson's R = "' num2cell(bb)];


% Save workspace and Export Data to Excel
newFilename = horzcat('Coactivity_',ani,'_',filename);

save(file,'-v7.3');
disp('saved');

% Original Raw Data
data = [array2table(t, 'VariableNames',{'Time(s)'}) data];
writetable(data,newFilename,'Sheet',filename);
disp('data');

% Binarised Signal
writetable(newData,newFilename,'Sheet','Binarised Signal');
disp('binarised signal');

% Binarised Coactivity
Means = (data.Properties.VariableNames)';
Means(1,:) = [];

Coactivity = [array2table(Means) Coactivity];
writetable(Coactivity,newFilename,'Sheet','Binarised Coactivity');
disp('binarised coactivity');

% Significant Coactivity
SignificantCoactivity = [array2table(Means) SignificantCoactivity];
SignificantCoactivity = [SignificantCoactivity;Hubs;PercentConnect];
writetable(SignificantCoactivity,newFilename,'Sheet','Significant Coactivity');
disp('significant coactivity');

% Pearson's R
writematrix(AvgPearsonR,newFilename,'Sheet',"Pearson's R",'Range','C2');
disp('pearsons r');

% Pearson's R Coefficients
MeansPR = (data.Properties.VariableNames)';
MeansPR(1,:) = [];

PearsonsR = data.Properties.VariableNames;
PearsonsR(:,1) = [];
PearsonsR = array2table(test,'VariableNames',PearsonsR);

PearsonsR = [array2table(MeansPR) PearsonsR];
writetable(PearsonsR,newFilename,'Sheet',"Pearson's R Coefficients");
disp('pearsons r coefficients');

% Significant Pearson's R
SigPearsonsR = data.Properties.VariableNames;
SigPearsonsR(:,1) = [];
SigPearsonsR = array2table(sign,'VariableNames',SigPearsonsR);

SigPearsonsR = [array2table(MeansPR) SigPearsonsR];
writetable(SigPearsonsR,newFilename,'Sheet',"Significant Pearson's R");
disp('significant pearsons r');

% Calcium Traces
%xlswritefig(CalciumTraces,newFilename,filename,'B2');

% Topographic Connectivity
%xlswritefig(TopographicConnectivity,newFilename,'Significant Coactivity','C5');

% Power Law
%xlswritefig(PowerLaw,newFilename,'Significant Coactivity','K5');

% Pearson's R
%xlswritefig(PearsonRFig,newFilename,"Pearson's R",'C4');

% Significant Pearson's R
%xlswritefig(SigPearsonRFig,newFilename,"Pearson's R",'I4');

close all

cd(rootdir);