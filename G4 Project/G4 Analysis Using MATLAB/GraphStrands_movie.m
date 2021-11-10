
%% DATA
AI = readtable("1kf1_Ab_Initio.csv");
C = readtable("1kf1.csv");
%Ab_Initio data
AIData = [AI{: , 1}, AI{: , 2}, AI{: , 3}, AI{: , 4}, AI{: , 5}, AI{: , 6} ];
%All shape parameters1
AllData = [C{: , 1}, C{: , 2}, C{: , 3}, C{: , 4}, C{: , 5}, C{: , 6} ];
%Organize data per GG pair
GGPairs = SeparateStrands(AllData);

%% Global data
global sizepergg params paramnum strandname legends strands v maxp minp lines time;
sizepergg = 9970;
params = {'Shift (dx)' ' Slide (dy)' ' Rise (dz)' ' Tilt (\tau)' ' Roll (\rho)' ' Twist (\Omega)'};
paramnum = 6 - 1;
strandname = {'2-3-4 Strand' '8-9-10 Strand' '14-15-16 Strand' '20-20-22 Strand'};
legends = {'g2g3' 'g3g4' ; 'g8g9' 'g9g10' ; 'g14g15' 'g15g16' ; 'g20g21' 'g21g22' };
strands = [1 3 5 7];
v = VideoWriter('Shift_Slide.avi');
maxp = max(AllData);
minp = min(AllData);
lines = [ line line ; line line ; line line ; line line ];
time = [];

%% Graph Settings and Video Settings
v.Quality = 100;
v.FrameRate = 30;
open(v);

%% Graph the different parameter 2D combinations with 1 plot per strand
datapoints = 10;
figure
for i = 1:paramnum
    for j = i+1: 6
        %Scale figure and do not render real-time
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        %set(gcf,'Visible','off')
        strandsMovie(GGPairs,i,j,datapoints - 1);
        break
        %Save figures individually
        %         file = strcat(file, '.jpg');
        %         exportgraphics(gcf,file, 'Resolution', 300);
        %         close(gcf)
    end
    break
end

% %% Create individual
%
% %Plot for Each Strand (i.e. 2 GG Pairs)
% figure
% tiledlayout(2,2)
% for i = 1:4
%     nexttile
%     strand = strands(i);
%     c = linspace(1,10,length(GGPairs(:,1,strand)));
%     scatter(GGPairs(:,1,strand),GGPairs(:,2,strand),.5)
%     hold on
%     scatter(GGPairs(:,1,strand+1),GGPairs(:,2,strand+1),.5)
%     legend(legends{i,1},legends{i,2})
%     title(strandname{i});
% end

%Deprecated Function: This function created various
function m = SeparateStrands(AllData)
%Start Matrix
GGPairs = zeros(9970,6,8);
j = 1;
size1 = 1;
size2 = 1;
size3 = 1;
size4 = 1;
size5 = 1;
size6 = 1;
size7 = 1;
size8 = 1;
for i = 1:size(AllData)
    switch j
        case 1
            GGPairs(size1,:,j) = AllData(i,:);
            size1 = size1 + 1;
        case 2
            GGPairs(size2,:,j) = AllData(i,:);
            size2 = size2 + 1;
        case 3
            GGPairs(size3,:,j) = AllData(i,:);
            size3 = size3 + 1;
        case 4
            GGPairs(size4,:,j) = AllData(i,:);
            size4 = size4 + 1;
        case 5
            GGPairs(size5,:,j) = AllData(i,:);
            size5 = size5 + 1;
        case 6
            GGPairs(size6,:,j) = AllData(i,:);
            size6 = size6 + 1;
        case 7
            GGPairs(size7,:,j) = AllData(i,:);
            size7 = size7 + 1;
        case 8
            GGPairs(size8,:,j) = AllData(i,:);
            size8 = size8 + 1;
    end
    j = j + 1;
    if (j == 9)
        j = 1;
    end
end
m = GGPairs;
end

function strandsMovie(GGPairs, i, j, datapoints)
global params strandname maxp minp sizepergg v;
%Determine number of data points based on frames

% Write Labels and Annotations
for k = 1:4
    subplot(2,2,k)
    %Labeling and Annotations
    title(strandname{k});
    xlabel(params{i})
    ylabel(params{j})
    axis([minp(i) maxp(i) minp(j) maxp(j)])
end

%Prepare the Vector with missing data. Where each row represents {file num, missing quantity}
missdat = [ 6 9 ; 7 (14 + 9) ; 12 (7 + 14 + 9) ];
missnum = 3;
mdatlast = uint8(1);
nextmiss = uint16(missdat(mdatlast,1)*500 -  missdat(mdatlast,2)); %This is the index of the start of the missing data

%Draw individual frames for each few data points. 
% Draw First Frame individually so as to simplify code in loop
startt = 1;
endd = startt + datapoints;
DrawStrands(GGPairs(startt:endd ,:,:), i, j,'Time (ns): 1',true);
writeVideo(v,getframe(gcf));
frametime = 1;

%Loop throught the rest of the frames
for frame = 2:(sizepergg/datapoints)
    %Data ranges
    startt = endd + 1;
    endd = startt + datapoints;
    frametime = frametime + 1;
    if (nextmiss <= endd && mdatlast <= missnum) % This could be omitted in a better program
        endd = nextmiss;
        timestr = strcat('Time (ns): ', string(frametime));
        DrawStrands(GGPairs(startt:endd ,:,:), i, j,timestr,false);
        frametime = missdat(mdatlast,1)*50;
        
        %Redefine endd for next frame
        mdatlast = mdatlast + 1;
        if (mdatlast <= missnum)
            nextmiss = uint16(missdat(mdatlast,1)*500 -  missdat(mdatlast,2));
        end
    else
        %Define time string
        timestr = strcat('Time (ns): ', string(frametime));
        DrawStrands(GGPairs(startt:endd ,:,:), i, j,timestr,false);
    end
    writeVideo(v,getframe(gcf));
end
close(v)
end

function DrawStrands(GGPairs, i, j, timestr ,first)
% Assign global to local variables for safety
global strands legends lines time;
if first == true
    for k = 1:4
        %Set Up
        subplot(2,2,k)
        strand = strands(k);
        lines(k,1) = line(GGPairs(:,i,strand),GGPairs(:,j,strand),'Marker', '.', 'LineStyle', 'none', 'Color', 'b', 'MarkerSize', 15);
        lines(k,2) = line(GGPairs(:,i,strand+1),GGPairs(:,j,strand+1),'Marker', '.', 'LineStyle', 'none', 'Color', 'r', 'MarkerSize', 15);
        legend(legends{k,1},legends{k,2}, 'AutoUpdate', 'off')
    end
    time = text('String', timestr , 'Position', [-2.4 3.5], 'FontSize', 15);
else
    for k = 1:4
        %Set Up
        subplot(2,2,k)
        
        strand = strands(k);
        set(lines(k,1), 'XData', GGPairs(:,i,strand));
        set(lines(k,1), 'YData', GGPairs(:,j,strand));
        set(lines(k,2), 'XData', GGPairs(:,i,strand+1));
        set(lines(k,2), 'YData', GGPairs(:,j,strand+1));
    end
    set(time, 'String', timestr);
end
end

function DrawAverage(draw,i,j)
end