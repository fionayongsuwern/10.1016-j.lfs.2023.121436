%% Data Preparation
% Make sure mean data has headers in the format 'Mean1', 'Mean2'...'Mean30', etc
% Make sure coordinate data has headers 'X' for column with x-coordinates and 'Y' for column with y-coordinates
close all
clear
clc


%% Inputs
file = 'exp1-i4';
ext = '.xlsx';
filename = dir(horzcat('*',file,'*',ext));
filename = filename.name;


%% Import data
data = readtable(filename);
Time = data(:,1);
data.Time_s_ = [];
numOfCols = width(data);
numOfRows = height(data);


%% Plot Animated Calcium Traces

t = table2array(Time);
fluorescenceY = table2array(data);

X = round(max(t),-2);
if X<max(t)
    if rem(X,200)==0
        X = X+200;
    elseif rem(X,200)>0
        X = X+100;
    end
end

YY = max(max(table2array(data)));
Y = round(YY,-2);
if Y<YY
    Y = Y+10000;
end

CalciumTraces = figure('Name','Calcium Traces','NumberTitle','off');

file1 = horzcat(file,'_Calcium Traces.avi');
CalciumAnimation = VideoWriter(file1);
CalciumAnimation.Quality = 100;
CalciumAnimation.FrameRate = 20;

open(CalciumAnimation);

for i = 1:1:numOfRows
    for j = 1:1:numOfCols
        plot(t(i),fluorescenceY(i,j))
        
        hold on
    end
    
    for j = 1:1:numOfCols
        plot(t(1:i),fluorescenceY(1:i,j))
        
        hold on
    end
    
    axis([0 X 0 Y])
    box off
    xlabel('Time','FontWeight','bold','FontSize',14)
    ylabel('Fluorescence','FontWeight','bold','FontSize',14)
    ax = gca;
    ax.YAxis.Exponent = 0;
    
    frame = getframe(gcf);
    writeVideo(CalciumAnimation,frame);
    
%     pause(0.1)
    
    if i~=length(t)
        clf
    end
end

close(CalciumAnimation);

