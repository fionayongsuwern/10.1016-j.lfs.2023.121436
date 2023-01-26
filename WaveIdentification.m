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


%% Import data
data = readtable(filename);
numOfCols = width(data);
numOfRows = height(data);

oldFolder = cd;
mkdir('Wave Identification');
newFolder = cd + "\" + 'Wave Identification';
cd(newFolder);


%% Plot Calcium Traces

% Plot calcium traces
CalciumTraces = figure('Name','Calcium Traces','NumberTitle','off');

for var = 1:1:numOfCols
    Mean = horzcat('Mean',num2str(var));
    
    figure(CalciumTraces);
    fluorescenceY = data.(Mean);
    plot(fluorescenceY);
    
    hold on
end

xticks(0:500:numOfRows)

ax = gca;
ax.XAxis.Label.String = 'Time';
ax.XAxis.Label.FontSize = 14;
ax.XAxis.Label.FontWeight = 'bold';
ax.YAxis.Label.String = 'Fluorescence';
ax.YAxis.Label.FontWeight = 'bold';
ax.YAxis.Label.FontSize = 14;
ax.YAxis.Exponent = 0;


%% Plot calcium traces after smoothing
t = 1:1:numOfRows;
t = t';

Signal = zeros(0,2);
SmoothCalciumTraces = figure('Name','Calcium Traces (Smoothed)','NumberTitle','off');
BinarisedCalciumTraces = figure('Name','Calcium Traces (Binarised)','NumberTitle','off');

for i = 1:1:numOfCols
    Mean = horzcat('Mean',num2str(i));
    
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
    
    waveStart = strfind(fluorescenceYii',[0 1]);
    waveStart = waveStart';
    waveStart = [waveStart,zeros(length(waveStart),1)+i];
    
    hold on
    
    figure(BinarisedCalciumTraces);

    plot(fluorescenceYii);
    
    hold on
    
    Signal = [Signal; waveStart];
end

figure(SmoothCalciumTraces);
file1 = horzcat(file,'_Smoothed Calcium Traces');
saveas(gcf,file1,'jpg');
saveas(gcf,file1,'bmp');

figure(BinarisedCalciumTraces);
file2 = horzcat(file,'_Binarised Calcium Traces');
saveas(gcf,file2,'jpg');
saveas(gcf,file2,'bmp');

Signal = array2table(Signal,'VariableNames',{'StartTime','Mean'});

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

file3 = horzcat(file,'_Example Processed Calcium Traces');
saveas(gcf,file3,'jpg');
saveas(gcf,file3,'bmp');

Signal = sortrows(Signal,{'StartTime'});

uniqueTime = Signal.StartTime;
uniqueTime = unique(uniqueTime);
Signal1 = {};

for k = 1:1:length(uniqueTime)
    
    l = find(Signal.StartTime==uniqueTime(k));
    l = Signal.Mean(l);
    l = {l'};
    
    Signal1 = [Signal1;l];    
end

Signal1 = [num2cell(uniqueTime) Signal1];
Signal1 = [{'StartTime','Means'};Signal1];


%% Save workspace and Export Data to Excel

newFilename = horzcat('Wave Identification_',file,ext);
newFile = horzcat('Wave Identification_',file);
save(newFile);

writecell(Signal1,newFilename);
xlswritefig(CalciumTraces,newFilename,'Figures','B2');
xlswritefig(SmoothCalciumTraces,newFilename,'Figures','L2');
xlswritefig(BinarisedCalciumTraces,newFilename,'Figures','B25');
xlswritefig(ExampleCalciumTraces,newFilename,'Figures','L25');

cd(oldFolder);