%% Data Preparation
close all
clear
clc


%% Import Data
wavetype = '.';

if length(wavetype)>1
    dataname = horzcat('Coactivity*',wavetype,'*.xlsx');
    
    data = dir(dataname);
    num = length(data);
    
    file = horzcat('Power Law for n=',num2str(num),' ',wavetype);
    filename = horzcat(file,'.xlsx');
    
elseif length(wavetype)==1
    dataname = horzcat('Coactivity*','.xlsx');
    
    data = dir(dataname);
    num = length(data);
    
    file = horzcat('Power Law for n=',num2str(num));
    filename = horzcat(file,'.xlsx');
end


%% Import data

% Import Significant Coactivity data
for i = 1:1:num
    Islet = horzcat('Islet',num2str(i));
    SignificantCoactivity.(Islet) = table;
       
    SignificantCoactivity.(Islet) = readtable(data(i).name,'Sheet','Significant Coactivity');
    SignificantCoactivity.(Islet) = table2array(SignificantCoactivity.(Islet)(1:end-3,2:end));
end


%% Calculate Coactivity and Power Law
threshold = 0.8;
connections = 30;

% Calculate Coactivity and Power Law for Pooled Data
Percent = [];
Hubs = [];
numOfHubs = [];
HubsSEM = [];
x = [];
y = [];
X = [];
Y = [];

for i = 1:1:num
    Islet = horzcat('Islet',num2str(i));
    numOfCols = length(SignificantCoactivity.(Islet));
    
    Percent = sum(SignificantCoactivity.(Islet)>threshold)/numOfCols*100;
    
    numOfHubs = sum(Percent>=connections);
    numOfHubs = numOfHubs/numOfCols*100;
    Hubs = [Hubs; numOfHubs];
    
    x = unique(Percent);
    y = zeros(1,length(x));
    
    for j = 1:1:length(x)
        y(j) = (sum(Percent==x(j)))/numOfCols*100;
    end
    
    X = [X x];
    Y = [Y y];
end

numOfHubs = mean(Hubs);
numOfHubs = horzcat('Percentage Hub Cells = ',num2str(numOfHubs),'%');
HubsSEM = std(Hubs)/sqrt(length(Hubs));
HubsSEM = horzcat('SEM = ',num2str(HubsSEM));
Threshold = horzcat('Coactivity Threshold = ',num2str(threshold));
MinCons = horzcat('Min Percentage Connections = ',num2str(connections),'%');

PowerLawFig = figure('Name','Power Law','NumberTitle','off');

try
    [slope, intercept, MSE, R2, S] = logfit(X,Y,'loglog');
    R = horzcat('R² = ',num2str(R2));
    %txt = [R;cellstr(numOfHubs)];
    PowerLaw = ['Power law';R;cellstr(numOfHubs);cellstr(HubsSEM);cellstr(Threshold);cellstr(MinCons)];
catch
    txt = ('Not enough data points to plot a graph.');
end

legend('off');
xlim([0 100]);
ylim([1 100]);
title('Power Law','FontSize',12,'Color','k');
xlabel('Number of Connections','FontSize',12,'Color','k','FontWeight','bold');
ylabel('Probability','FontSize',12,'Color','k','FontWeight','bold');
yticks([1 10 100]);
ax.YAxis.Exponent = 0;
xticks([0 100]);
ax.XAxis.Exponent =0;
set(gca,'box','off');

saveas(gcf,file,'jpg');
saveas(gcf,file,'bmp');

% Calculate Coactivity and Power Law for each Islet
BetaHubs = {'Islet' 'R²' 'Percentage of Hub Cells', 'Hubs'};

for i = 1:1:num
    Islet = horzcat('Islet',num2str(i));
    HubCells.(Islet) = table;
    
    HubCells.(Islet) = readtable(data(i).name,'Sheet','Significant Coactivity');
    HubCells.(Islet) = HubCells.(Islet)(1:end-3,2:end);
    numOfCols = width(HubCells.(Islet));
    
    IsletPercent = sum(table2array(HubCells.(Islet))>=threshold)/(numOfCols-1)*100;
    IsletHub = IsletPercent>=connections;
    numOfIsletHub = sum(IsletPercent>=connections);
    numOfIsletHub = numOfIsletHub/numOfCols*100;
    
    isletX = unique(IsletPercent);
    isletY = zeros(1,length(isletX));
    
    for j = 1:1:length(isletX)
        isletY(j) = (sum(IsletPercent==isletX(j)))/numOfCols*100;
    end
    
    IsletPowerLaw = figure;
    
    try
        [slope, intercept, MSE, R2, s] = logfit(isletX,isletY,'loglog');
    catch
        R2 = NaN;
    end
    
    close(IsletPowerLaw);
    
    HubCellTable = [HubCells.(Islet);num2cell(IsletPercent)];
    HubCells.(Islet) = table;
    IsletHubCells = table;
    
    for j = 1:1:numOfCols
        if IsletHub(1,j) == 1
            IsletHubCells = HubCellTable(end,j);
            
            HubCells.(Islet) = [HubCells.(Islet) IsletHubCells];
        end
    end
        
    BetaHubs = [BetaHubs; data(i).name R2 numOfIsletHub cellstr(strjoin(HubCells.(Islet).Properties.VariableNames))];
end



%% Calculate Connectivity

% Calculate Hub Cells Connactivity
HubsConnectivity = {'Islet' 'Hub Cells' 'Percent Coactivity'};

for i = 1:1:num
    Islet = horzcat('Islet',num2str(i));
    HubCellsConnect = table;
    HubCellsPercentConnect = table;
    
    HubCellsConnect = HubCells.(Islet).Properties.VariableNames';
    
    if isempty(HubCellsConnect)
        HubCellsConnect = cell(1,1);
        HubCellsPercentConnect = cell(1,1);
    else
        HubCellsPercentConnect = table2cell(HubCells.(Islet))';
    end
    
    HubsConnect = cell(length(HubCellsConnect),1);
    HubsConnect{1,1} = data(i).name;
    
    HubsConnectivity = [HubsConnectivity; HubsConnect HubCellsConnect HubCellsPercentConnect];
end

AvgHubsConnect = mean(cell2mat(HubsConnectivity(2:end,3)));
AvgHubsConnectivity = ['Average Connectivity of Hub Cells = ' cellstr(horzcat(num2str(AvgHubsConnect),'%'))];

% Calculate Wave Originator Cells Connectivity
WaveOriginConnectivity = {'Islet' 'Wave Originator Cell' 'Percent Coactivity'};

for i = 1:1:num
    Islet = horzcat('Islet',num2str(i));
    WaveOriginCells.(Islet) = table;
    WaveOriginCellsConnect = table;
    WaveOriginCellsPercentConnect = table;
    
    WOCells.(Islet) = readtable(data(i).name,'Sheet','Significant Coactivity');
    WOCells.(Islet) = WOCells.(Islet)(end-1:end,2:end);
    
    WOTable = table2array(WOCells.(Islet));
    numOfCols = width(WOCells.(Islet));
    
    for j = 1:1:numOfCols
        if WOTable(2,j) == 1
            WaveOrigin = WOCells.(Islet)(end-1,j);
            
            WaveOriginCells.(Islet) = [WaveOriginCells.(Islet) WaveOrigin];
        end
    end
    
    WaveOriginCellsConnect = WaveOriginCells.(Islet).Properties.VariableNames';
    
    if isempty(WaveOriginCellsConnect)
        WaveOriginCellsConnect = cell(1,1);
        WaveOriginCellsPercentConnect = cell(1,1);
    else
        WaveOriginCellsPercentConnect = table2cell(WaveOriginCells.(Islet))';
    end
    
    WaveOriginConnect = cell(length(WaveOriginCellsConnect),1);
    WaveOriginConnect{1,1} = data(i).name;
    
    WaveOriginConnectivity = [WaveOriginConnectivity; WaveOriginConnect WaveOriginCellsConnect WaveOriginCellsPercentConnect];
end

AvgWaveOriginConnect = mean(cell2mat(WaveOriginConnectivity(2:end,3)));
AvgWaveOriginConnectivity = ['Average Connectivity of Wave Originator Cells = ' cellstr(horzcat(num2str(AvgWaveOriginConnect),'%'))];

% Calculate Follower Cells Connectivity
FollowersConnectivity = {'Islet' 'Follower Cell' 'Percent Coactivity'};

for i = 1:1:num
    Islet = horzcat('Islet',num2str(i));
    FollowerCells.(Islet) = table;
    FollowerCellsConnect = table;
    FollowerCellsPercentConnect = table;
    
    FWCells.(Islet) = readtable(data(i).name,'Sheet','Significant Coactivity');
    FWCells.(Islet) = FWCells.(Islet)(end-2:end,2:end);
    
    FWTable = table2array(FWCells.(Islet));
    numOfCols = width(FWCells.(Islet));
    
    for j = 1:1:numOfCols
        if (FWTable(1,j)==0 && FWTable(3,j)==0)
            Followers = FWCells.(Islet)(end-1,j);
            
            FollowerCells.(Islet) = [FollowerCells.(Islet) Followers];
        end
    end
    
    FollowerCellsConnect = FollowerCells.(Islet).Properties.VariableNames';
    
    if isempty(FollowerCellsConnect)
        FollowerCellsConnect = cell(1,1);
        FollowerCellsPercentConnect = cell(1,1);
    else
        FollowerCellsPercentConnect = table2cell(FollowerCells.(Islet))';
    end
    
    FollowersConnect = cell(length(FollowerCellsConnect),1);
    FollowersConnect{1,1} = data(i).name;
    
    FollowersConnectivity = [FollowersConnectivity; FollowersConnect FollowerCellsConnect FollowerCellsPercentConnect];
end

AvgFollowersConnect = mean(cell2mat(FollowersConnectivity(2:end,3)));
AvgFollowersConnectivity = ['Average Connectivity of Follower Cells = ' cellstr(horzcat(num2str(AvgFollowersConnect),'%'))];

% Calculate Non-Hub Cells Connectivity
NHConnectivity = {'Islet' 'Non-Hub Cell' 'Percent Coactivity'};

for i = 1:1:num
    Islet = horzcat('Islet',num2str(i));
    NHCells.(Islet) = table;
    NHCellsConnect = table;
    NHCellsPercentConnect = table;
    
    NHCells.(Islet) = readtable(data(i).name,'Sheet','Significant Coactivity');
    NHCells.(Islet) = NHCells.(Islet)(end-1,2:end);
    
    NHTable = table2array(NHCells.(Islet));
    numOfCols = width(NHCells.(Islet));
    
    for j = 1:1:numOfCols
        if NHTable(1,j)==0
            NH = NHCells.(Islet)(end,j);
            
            NHCells.(Islet) = [NHCells.(Islet) NH];
        end
    end
    
    NHCellsConnect = NHCells.(Islet).Properties.VariableNames';
    
    if isempty(NHCellsConnect)
        NHCellsConnect = cell(1,1);
        NHCellsPercentConnect = cell(1,1);
    else
        NHCellsPercentConnect = table2cell(NHCells.(Islet))';
    end
    
    NHConnect = cell(length(FollowerCellsConnect),1);
    NHConnect{1,1} = data(i).name;
    
    NHConnectivity = [NHConnectivity; NHConnect NHCellsConnect NHCellsPercentConnect];
end

AvgNHConnect = mean(cell2mat(NHConnectivity(2:end,3)));
AvgNHConnectivity = ['Average Connectivity of Non-Hub Cells = ' cellstr(horzcat(num2str(AvgNHConnect),'%'))];


%% Construct bar graphs
barData = [AvgHubsConnect AvgWaveOriginConnect AvgFollowersConnect];

AvgHub = cell2mat(HubsConnectivity(2:end,3));
AvgWaveOrigin = cell2mat(WaveOriginConnectivity(2:end,3));
AvgFollower = cell2mat(FollowersConnectivity(2:end,3));

errLow = [min(cell2mat(HubsConnectivity(2:end,3))) min(cell2mat(WaveOriginConnectivity(2:end,3))) min(cell2mat(FollowersConnectivity(2:end,3)))];
errLow = barData - errLow;
errHigh = [max(cell2mat(HubsConnectivity(2:end,3))) max(cell2mat(WaveOriginConnectivity(2:end,3))) max(cell2mat(FollowersConnectivity(2:end,3)))];
errHigh = errHigh - barData;

BarFig = figure('Name','Average Coactivity','NumberTitle','off');

barX = (1:3);
BarFig1 = bar(barX,barData,'FaceColor','flat');

yticks(0:20:100);
ylabel('Percentage Coactivity')

BarFig1.CData(1,:) = [1 0 0];
BarFig1.CData(2,:) = [0 1 0];
BarFig1.CData(3,:) = [0 0 1];

hold on

AvgHubX = zeros(length(AvgHub),1)+1;
AvgHubPlot = scatter(AvgHubX,AvgHub,10,'k','filled','o');

hold on

AvgWaveOriginX = zeros(length(AvgWaveOrigin),1)+2;
AvgWaveOriginPlot = scatter(AvgWaveOriginX,AvgWaveOrigin,10,'k','filled','o');

hold on

AvgFollowerX = zeros(length(AvgFollower),1)+3;
AvgFollowerPlot = scatter(AvgFollowerX,AvgFollower,10,'k','filled','o');

set(gca,'xticklabel',{'Hubs','Wave Origins','Followers'});

box off

file1 = horzcat(file,'_Average Coactivity Chart (with spread)');
saveas(BarFig,file1,'jpg');
saveas(BarFig,file1,'bmp');

BarFigErr = figure('Name','Average Coactivity','NumberTitle','off');

barX = (1:3);
BarFigErr1 = bar(barX,barData,'FaceColor','flat');

yticks(0:20:100);
ylabel('Percentage Coactivity')

BarFigErr1.CData(1,:) = [1 0 0];
BarFigErr1.CData(2,:) = [0 1 0];
BarFigErr1.CData(3,:) = [0 0 1];

hold on

er = errorbar(barX,barData,errLow,errHigh,'Color','k','LineStyle','none','LineWidth',1,'CapSize',25);

set(gca,'xticklabel',{'Hubs','Wave Origins','Followers'});

box off

file2 = horzcat(file,'_Average Coactivity Chart (with error bar)');
saveas(BarFigErr,file2,'jpg');
saveas(BarFigErr,file2,'bmp');


%% T-Test
TTest = [cellstr('T-Test');'Significance';'Alpha';'T-Stat';'df';'sd'];

[h,p,ci,stats] = ttest2(AvgHub,AvgWaveOrigin,'Alpha',0.001);
TTest1 = [cellstr('Hubs vs Wave Origins');p;0.001;stats.tstat;stats.df;stats.sd];

[h,p,ci,stats] = ttest2(AvgHub,AvgFollower,'Alpha',0.001);
TTest2 = [cellstr('Hubs vs Followers');p;0.001;stats.tstat;stats.df;stats.sd];

[h,p,ci,stats] = ttest2(AvgWaveOrigin,AvgFollower,'Alpha',0.001);
TTest3 = [cellstr('Wave Origins vs Followers');p;0.001;stats.tstat;stats.df;stats.sd];

TTest = [TTest TTest1 TTest2 TTest3];


%% Save workspace and Export Data to Excel
save(file);

writecell(PowerLaw,filename,'Sheet','Power Law','Range','C2:C7');
xlswritefig(PowerLawFig,filename,'Power Law','B9');

writecell(TTest,filename,'Sheet','Power Law','Range','J2:M7');
xlswritefig(BarFig,filename,'Power Law','I9');
xlswritefig(BarFigErr,filename,'Power Law','O9');

writecell(BetaHubs,filename,'Sheet','Beta Hub Cells','Range','A1');

writecell(AvgHubsConnectivity,filename,'Sheet','Hub Cells Connectivity','Range','F5:G5');
writecell(HubsConnectivity,filename,'Sheet','Hub Cells Connectivity','Range','A1');

writecell(AvgWaveOriginConnectivity, filename,'Sheet','Wave Originator Connectivity','Range','F5:G5');
writecell(WaveOriginConnectivity, filename,'Sheet','Wave Originator Connectivity','Range','A1');

writecell(AvgFollowersConnectivity, filename,'Sheet','Follower Cells Connectivity','Range','F5:G5');
writecell(FollowersConnectivity, filename,'Sheet','Follower Cells Connectivity','Range','A1');