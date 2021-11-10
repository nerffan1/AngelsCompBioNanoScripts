AI = readtable("1K8P_Ab_Initio.csv");
table = readtable("1K8P.csv");
%Ab_Initio Data
AIData = [AI{: , 1}, AI{: , 2}, AI{: , 3}, AI{: , 4}, AI{: , 5}, AI{: , 6} ];
%All Data
AllData = [table{: , 1}, table{: , 2}, table{: , 3}, table{: , 4}, table{: , 5}, table{: , 6} ];
%Initiate Pool for multithreading
%pool = parpool('threads');

%Define Global Variables
global maxp minp
maxp = max(AllData);
minp = min(AllData);

%% TILED LAYOUT SILHOUETTE
tiledlayout(3,3)
%K mean analysis from MATLAB

for centers = 2:10
idx = kmeans(AIData,centers,'MaxIter',200,'Replicates',250,'Display','off','Options',statset('UseParallel',1));
nexttile
[silh,h] = silhouette(AIData,idx);
xlabel('Silhouette Value')
ylabel('Cluster')
%Calculate score and load to vector
Score = mean(silh); 

title('Silhouette Score:', Score)
end

%% SHIFT - SLIDE - RISE 
[idx,table] = kmeans(AllData,2,'MaxIter',200,'Replicates',250,'Display','off','Options',statset('UseParallel',1));
cluster1 = idx==1;
cluster2 = idx==2;
%New Figure
figure
tiledlayout(1,2)
nexttile
%Plot All data
scatter3(AllData(cluster1,1),AllData(cluster1,2), AllData(cluster1,3),.1,'r.')
hold on
scatter3(AllData(cluster2,1),AllData(cluster2,2),AllData(cluster2,3),.1,'g.')
hold on 
%Plot Ab-Initio Data
scatter3(AIData(:,1),AIData(:,2),AIData(:,3),10,'k.')
hold on
%Plot centroids --- Could I optimize this code?
plot3(table(1,1),table(1,2),table(1,3),'kx','LineWidth',1,'MarkerSize',20)
hold on
plot3(table(2,1),table(2,2),table(2,3),'kx','LineWidth',1,'MarkerSize',20)
title ('dx vs. dy (1K8P)')
xlabel('Shift (dx)')
ylabel('Slide (dy)')
zlabel('Rise (dz)')
%% PITCH - ROLL - TWIST
nexttile
%Plot All data
scatter3(AllData(cluster1,4),AllData(cluster1,5), AllData(cluster1,6),.1,'r.')
hold on
scatter3(AllData(cluster2,4),AllData(cluster2,5),AllData(cluster2,6),.1,'g.')
hold on 
%Plot the Ab-Initio Data
scatter3(AIData(:,4),AIData(:,5),AIData(:,6),10,'k.')
hold on
%Plot Centroids
plot3(table(1,4),table(1,5),table(1,6),'kx','LineWidth',1,'MarkerSize',17)
hold on
plot3(table(2,4),table(2,5),table(2,6),'kx','LineWidth',1,'MarkerSize',17)
title ('\rho vs. \omega (1K8P)')
xlabel('PITCH')
ylabel('Roll (\rho)')
zlabel('Twist (\omega)')

%% Create a count per G-Pair, and count how many times it falls 
figure
pairs = zeros(8,2);
pair = 1; 
for i = 1:size(idx)
    if (idx(i) == 1)
        pairs(pair,1) = pairs(pair,1) + 1;
    else
        pairs(pair,2) = pairs(pair,2) + 1;
    end
    %Change the pair index
    if (pair == 8)
        pair = 0;
    end
    pair = pair + 1;
end
% Create categories (must reorder to make sure order is not changed)
X = categorical({'g3g4','g4g5','g9g10','g10g11','g15g16','g16g17','g21g22','g22g23'});
X = reordercats(X,{'g3g4','g4g5','g9g10','g10g11','g15g16','g16g17','g21g22','g22g23'});
b = bar(X, pairs);
b(1,1).FaceColor = 'r';
b(1,2).FaceColor  = 'g';
title('Clusters per GG-Pair substructure')
xlabel('GG-Pair Structure')
ylabel('Counts')
%% Create Individual plots for each 
%Description: Pass vector data and create various 2D plots for the base
% parameters, and also save them into individual figures
params = {'Shift (dx)' 'Slide (dy)' 'Rise (dz)' 'Tilt (tau)' 'Roll (rho)' 'Twist (omega)'}; 
paramnum = 6 - 1;
%Initiate Graphics objects
close all
figure
set(gcf,'Visible','off')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
clustscat1 =  scatter(AllData(cluster1,1),AllData(cluster1,1),7,'r','.');
hold on 
clustscat2 = scatter(AllData(cluster2,1),AllData(cluster2,1),7, 'g','.');
hold on
cent1 = line(1,1,'Marker', 'x', 'LineStyle', 'none', 'Color', 'b', 'MarkerSize', 20);
cent2 = line(1,1,'Marker', 'x', 'LineStyle', 'none', 'Color', 'b', 'MarkerSize', 20);
hold on
AIscat = line(AIData(:,1),AIData(:,1),'Marker', 'o', 'LineStyle', 'none', 'Color', 'k', 'MarkerSize', 1.5);
axes = set(gcf,'CurrentAxes');
for i = 1:paramnum 
    for j = i: 6
        axes.XLim = [minp(i) maxp(i)];
        axes.YLim = [minp(j) maxp(j)];
        %Plot each cluster separately
        clustscat1.XData = AllData(cluster1,i);
        clustscat1.YData = AllData(cluster1,j);
        clustscat2.XData = AllData(cluster2,i);
        clustscat2.YData = AllData(cluster2,j);
        
        %Plot Centroids
        cent1.XData = table(1,i);
        cent1.YData = table(1,j);
        cent2.XData = table(2,i);
        cent2.YData = table(2,j);
        
        %Plot Ab-Initio Data
        AIscat.XData = AIData(:,i);
        AIscat.YData = AIData(:,j);
          
        %Labeling
        filename = strcat(params{i},' vs. ',params{j});
        title(strcat("1K8P: ", filename));
        xlabel(params{i})
        ylabel(params{j})
  
        %Save figures
        filename = strcat('1K8P_',filename, '.png');
        exportgraphics(gcf,filename, 'Resolution', 300);
    end
end

%% Create a Figure of XYZ using different symbols for each GG pair

delete(pool)