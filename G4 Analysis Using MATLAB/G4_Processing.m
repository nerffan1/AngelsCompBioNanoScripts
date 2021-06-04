%Description: This script is meant to aid the analysis of the G4 structures
% and holds various functions to performs different types of analysis,
% particularly engaged with unsupervised machine learning algorithms.
%Author: Angel-Emilio Villegas Sanchez
%LAST UPDATE: 06/01/21

global KM_pair %This variable is a clustering of individual strands of data
global Features params TI TI_Length minpts pair pairsize Graph_Params pairnum seven_ftrs G4
global Indices AI_Ind menu KM_All WorkDat G4_workspace

G4 = input('Which G4 would you like to analyze:\n-1KF1\n-1K8P\n','s');
if(LoadFiles())
    load(G4_workspace + "_Workspace.mat")
end

%% User Input
menu = input('Would you like a (1) comprehensive analysis, or (2) per-strand analysis');

%Ask for what GG Pair we want to analyze
Pairs_to_analyze = input("What pair(s) would you like to analyze? ");
params
Graph_Params = input("\nWhat parameters would you like to graph? ");
fprintf('You''ve chosen to graph %s and %s\n', params(Graph_Params(1)), params(Graph_Params(2)));
if (menu == 1)
    WorkDat = Features;
    MenuOpt()
else
    for i = 1:length(Pairs_to_analyze)
        pairnum = Pairs_to_analyze(i);
        WorkDat = Features(Indices(:,pairnum),:);
        fprintf('\nTo perform on Pair %i\n',pairnum);
        MenuOpt();
    end
end

%% Additional Functions
%Description: This function creates the menu system
function [] = MenuOpt()
%Ask for options per pair
Input = -2;
while (Input ~= -1)
    fprintf('1)K-Dist Graph\n2)GaugeEpsilons\n3)Gauge Epsilons (No files)\n4)DBSCAN with new TI\n');
    fprintf('5)Dendogram \n6)\n7)K-Means clustering\n');
    fprintf('8)Rolling Average\n9)Time Scatter Plots\n10)Figures for Paper\n11)DBSCAN\n');
    Input = input('Choice: ');
    
    switch Input
        case 1
            PlotK_Dist()
        case 2
            GaugeEpsilon()
        case 3
            GaugeIndivEp()
        case 4
            DBCluster()
        case 5
            Dendrograms()
        case 6
            
        case 7
            Kmeans()
        case 8
            RollingAverage()
        case 9
            TimeScatter()
        case 10
            PaperFigure()
        case 11
            
    end
end
end

%Description: Plot's the smallest Kth Distance for every point n the data.
function [] = PlotK_Dist()
global minpts pair menu WorkDat
figure
ylabel(num2str(minpts) + "th Smallest Distances to points")
xlabel("Number of points ")
kD = pdist2(WorkDat,WorkDat,'euc','Smallest',minpts);
if (menu == 1)
    title("K-Distance Graph: All data")
else
    title("K-Distance Graph (Pair " + num2str(pair) + ")")
end
plot(sort(kD(end,:)));
end

%Description: This function does a DBSCAN on either strands or full data
function [] = DBCluster()
global pair minpts Graph_Params TI AI_Ind menu WorkDat

eps = input('What epsilon would you like?');
labels = dbscan(pair,eps,minpts);
gscatter(WorkDat(:,Graph_Params(1)),WorkDat(:,Graph_Params(2)),labels)
title(params(1) + " vs. " + params(2))
xlabel(params(1))
ylabel(params(2))

if (menu == 1)
    %Plot the TI over time. Find overlap of indices
    DAT_Ind = 1:length(WorkDat);
    DAT_Ind = DAT_Ind(labels == 1); %Get the density connected values
    Conn_Ind = intersect(DAT_Ind,AI_Ind);
    TI_Ind = (AI_Ind == Conn_Ind');
    TI_Ind = any(TI_Ind);
    figure
    plot(1:length(Conn_Ind),TI(TI_Ind))
end
end

% Description: This function tests a range of values between 2 epsilons.
% The epsilon is part the parameter that's changed in the DBSCAN clustering
% method. It creates various files with statistics regarding each
% clustering, and separate files with more detailed information of each
% clustering (corresponding to a unique epsilon).
% NEEDS REVISION TO ACCOUNT FOR comprehensive/per strand
function [] = GaugeEpsilon()
global pair minpts pairsize TI pairnum

%Input
[eps , iter] = epsilonInput();

% Allocate Table Data for major data
Epsilon = linspace(eps(1),eps(2),iter)';
Noise = zeros(iter,1);
Cluster_Count = zeros(iter,1);
Largest_Cluster_Percentage = zeros(iter,1);
Ave_TI_Largest = zeros(iter,1);
Std_TI_Largest = zeros(iter,1);

for i = 1:iter
    
    % Perform DBSCAN using Epsilon(i)
    labels = dbscan(pair,Epsilon(i),minpts);
    labels_u = unique(labels);
    labels_count = length(labels_u);
    
    %Allocate data for table with specific cluster data
    TI_Points = zeros(labels_count,1);
    Ave_TI = zeros(labels_count,1);
    Std_TI = zeros(labels_count,1);
    
    %Find number of clusters and noise
    Cluster_Count(i) = length(labels_u);
    Noise(i) = sum(labels==-1)/pairsize;
    
    %Find Largest Cluster and also save data regarding each cluster
    max = 0;
    for j = 1:Cluster_Count(i)
        
        %Save additional Info
        C = IntersectionTI(labels, labels_u(j));
        TI_Subset = TI(C);
        TI_Length_Subset = length(TI_Subset);
        TI_Points(j) = TI_Length_Subset;
        Ave_TI(j) = mean(TI_Subset);
        Std_TI(j) = std(TI_Subset);
        
        %Max Cluster
        count = sum(labels == labels_u(j));
        if (count > max)
            max = count;
            Ave_TI_Largest(i) = Ave_TI(j);
            Std_TI_Largest(i) = Std_TI(j);
        end
    end
    
    %Save more values
    Largest_Cluster_Percentage(i) = max/pairsize;
    
    %Save table
    Clusters_Info = table(labels_u, TI_Points, Ave_TI, Std_TI);
    filename = strcat("Clusters_Pair_",int2str(pairnum),'_',int2str(minpts),"_",int2str(i));
    save(filename,'Clusters_Info')
end

EpsTable = table(Epsilon, Noise, Cluster_Count,Largest_Cluster_Percentage,Ave_TI_Largest,Std_TI_Largest);
filename = strcat('EpsTable_Pair_',int2str(pairnum),'_',int2str(minpts),'_minpts');
save(filename,'EpsTable')

end
function [] = GaugeIndivEp()
global pair minpts pairsize TI pairnum

%Input
[eps , iter] = epsilonInput();

% Allocate Table Data for major data
Epsilon = linspace(eps(1),eps(2),iter)';
Noise = zeros(iter,1);
Cluster_Count = zeros(iter,1);
Largest_Cluster_Percentage = zeros(iter,1);
Ave_TI_Largest = zeros(iter,1);
Std_TI_Largest = zeros(iter,1);

for i = 1:iter
    
    % Perform DBSCAN using Epsilon(i)
    labels = dbscan(pair,Epsilon(i),minpts);
    labels_u = unique(labels);
    
    %Find number of clusters and noise
    Cluster_Count(i) = length(labels_u);
    Noise(i) = sum(labels==-1)/pairsize;
    
    %Find Largest Cluster and also save data regarding each cluster
    max = 0;
    for j = 1:Cluster_Count(i)
        
        %Save additional Info
        C = IntersectionTI(labels, labels_u(j));
        TI_Subset = TI(C);
        
        %Max Cluster
        count = sum(labels == labels_u(j));
        if (count > max)
            max = count;
            Ave_TI_Largest(i) = mean(TI_Subset);
            Std_TI_Largest(i) = std(TI_Subset);
        end
    end
    
    %Save more values
    Largest_Cluster_Percentage(i) = max/pairsize;
    
end

EpsTable = table(Epsilon, Noise, Cluster_Count,Largest_Cluster_Percentage,Ave_TI_Largest,Std_TI_Largest);
filename = strcat('EpsTable_Pair_',int2str(pairnum),'_',int2str(minpts),'_minpts');
save(filename,'EpsTable')

end

%Description: This function prompts user for input in DBSCAN
function [eps , iter] = epsilonInput()
eps = input("Please tells us the ranges for Epsilon in an n x 2 matrix: ");
iter = input("Please tells us how many iterations to do: ");
end

%Description: This function takes the indices of a clustering performed on
% a superset of data of the TI values. This function finds the indices that
%correspond ,
%MUST UPDATE: Use method in DBCluster() for concise code
function C = IntersectionTI(labels,clusnum)
global pairnum;
%Labels is a logical array subset of Original (OG) data, but we need
%the index numbers. We must convert this to an array of actual indices
DIndices = zeros(1,1000);
DBi = 1;
for i = 1:length(labels)
    if (labels(i) == clusnum)
        DIndices(DBi) = i; %DBIndices belong to a vector of length 9970
        DBi = DBi + 1;
    end
end

%Convert DIndices to OG indices
DIndices = pairnum + (8*(DIndices-1));

%Indices of Ab-Initio in OG Data
Pairab = (49*8 + pairnum):100*8:79760;

%The intersection of the indices corresponding to Ab-Initio data and the
%ones found within chosen cluster
C = intersect(Pairab,DIndices);
% Convert to Ab_Initio indices in order to find the Ab-Initio points in the
% TI data
C = ((C-pairnum-392)/100) + pairnum;

end

%DEPRECATED CODE
function [] = PlotTI()
global seven_ftrs;

% % Plot the TI values within cluster of choice
C = IntersectionTI(labels, 1);
figure
scatter(seven_ftrs(C,1),seven_ftrs(C,7));
title("TI Values in Noise")
ylabel("TI Value")
xlabel("Data Count")

% % Plot OG TI points
% figure
% PairAb = (49*8 + pairnum):100*8:79760;
% PairAb  = ((PairAb - pairnum - 392)/100) + pairnum;
% plot(1:length(PairAb),TI(PairAb));
% title("TI Values (Original Values)")
% ylabel("TI Value")
% xlabel("Data Count")
end

%Description: Dendogram
function [] = Dendrograms()
global pairnum Graph_Params WorkDat
figure
link = linkage(WorkDat, 'average');
dendrogram(link);
cutoff = input('What cutoff would you like? ');
dendrogram(link,'ColorThreshold',cutoff);

%Create cutoff and scatter plot
figure
T = cluster(link,'Cutoff',cutoff,'Criterion','distance');
scatter(WorkDat(:,Graph_Params(1)),WorkDat(:,Graph_Params(2)),20,T,'.');
if (menu==1)
    title("Dendogram for All Data")
else
    title("Dendogram for pair " + num2str(pairnum))
end

end

%Description: This function performs a K-means clustering with various k's
%with silhouette scores to see if there's any
function [] = Kmeans()
global seven_ftrs TI pairnum KM_pair G4 KM_All WorkDat Graph_Params menu AI_Ind Features

%Create time column for data
if (G4 == "1KF1")
    time = [1:2991 3001:3486 3501:5993 6001:10000]' ;
else
    time = (1:10000)' ;
end

inp = input('1) Cluster Plots\n2) silh\n3) Cluster Evaluation\n4) Plot TI Scatter\n5) K-Means TI Histograms\n\nChoice: ');

% 2 clusters since this maximized the Sillhouette scores before.
if (menu==1)
    clusts = input("What k would you like for your clusters (1-16)? ");
    idx = KM_All(:,clusts);
else
    clusts = 2;
    idx = KM_pair(:,pairnum);
end

switch inp
    case 1
        ClusterScat()
    case 2
        silh()
    case 3
        KMeansTimeEval();
    case 4
        TI_Plot_KMeans()
    case 5
        TI_Hist_KMeans()
end

%Description: Simple scatter plot of clusters.
    function [] = ClusterScat()
        figure
        for i = 1:clusts
            clst_ind = (idx == i);
            scatter(WorkDat(clst_ind,Graph_Params(1)),WorkDat(clst_ind,Graph_Params(2)))
            hold on
        end
        if (menu==1)
            title("K-Means clustering (k=2) for all data")
        else
            title("K-Means clustering (k=2) - GG_" + num2str(pairnum))
        end
    end

%Description: Sillhouette score for the K-Means clustering
    function [] = silh()
        figure
        [silh,~] = silhouette(WorkDat,idx);
        xlabel('Silhouette Value')
        ylabel('Cluster')
        %Calculate score and load to vector
        Score = mean(silh);
        title('Silhouette Score:', Score)
    end

%Description: Evaluate the continuity of the k-means clustering
    function time_eval = KMeansTimeEval()
        time_eval = cell(clusts-1,1);
        for centers = 2:clusts
            %Loop through each cluster (i.e. amount of centers) and
            %evaluate the continuity of timestamps by checking the diff()
            %of a sorted set of timestamps (which it already is).
            ind_eval = cell(centers,1);
            for clust = 1:centers
                time = time(idx == clust);         
                figure
                plot(1:length(time),time)
                hold on
                time = diff(time);
                plot(1:length(time),time.*15)
                if (menu == 1)
                    title("K-Means (k= " + num2str(cluster) + ") of all data " + "(cluster " + num2str(clust) + ")")
                else
                    title("K-Means (k=2) of GG_" + num2str(pairnum) + "(cluster " + num2str(clust) + ")")
                end
                unq = unique(time);
                unq_l = length(unq);
                unq = [unq zeros(unq_l,1)];
                
                %Sum the number of times the diff is found
                for k = 1:unq_l
                    unq(k,2) = sum(time == unq(k));
                end
                
                %Then add this to a cell array
                ind_eval(clust) = {unq};
            end
            time_eval(centers - 1) = {ind_eval};
        end
    end

%Description: Plots the TI values  of each clustering
    function [] = TI_Plot_KMeans()
            for i = 1:clusts
                DAT_Ind = 1:length(WorkDat);
                DAT_Ind = DAT_Ind(KM_All(:,clusts) == i); %Get the density connected values
                Conn_Ind = intersect(DAT_Ind,AI_Ind);
                C = (AI_Ind == Conn_Ind');
                C = any(C);
                
                figure
                plot(1:sum(C),seven_ftrs(C,7));
                title("TI Values in Noise (k=" + int2str(i) + ")" + ", Cluster " + int2str(i) )
                ylabel("TI Value")
                xlabel("Data Count")
                figure
                scatter(seven_ftrs(C,3),seven_ftrs(C,4));
                title("K-Means Clustering, (k=" + int2str(i) + ")" + ", Cluster " + int2str(i) )
            end
    end

%Description: Histogram based of TI values for each clustering.
%NOTE: NEEDS REVISION for same reason above
    function [] = TI_Hist_KMeans()
        for centers = 2:clusts
            for i = 1:centers
                C = IntersectionTI(idx, i);
                figure
                histogram(TI(C),10)
                title("Histogram of Cluster " + num2str(i) + " in GG" + num2str(pairnum))
            end
        end
    end
end

%Description: This function graphs rolling averages for subsets of time of
%the data points. The window is user-defined
%REVISION: The Average values of the parameter is the same for every
%'different' pair since I'm not properly getting the TI values of each
%different GG pair
function [] = RollingAverage()
global seven_ftrs pairnum

window = input('\nWhat''s the window of points you''d like in the average? ');
pm = input('/nWhat parameter do you want to graph? ');

%Allocation of resources
i = pairnum:8:800;

MM = movmean(seven_ftrs(i,pm),window);

% for i = 1:length(ranges)-1
%     for j = (ranges(i)+1):ranges(i+1)
%         sum = sum + seven_ftrs(j,pm);
%     end
%     roll_i = roll_i + 1;
%     roll(roll_i) = sum/window;
%     sum = 0;
% end


%Plot the rolling average
figure
plot(1:length(MM),MM)
title("TI Rolling Average of " + num2str(pairnum))
ylabel("TI Average Value")
end

%Description: Scatter plots with color labels denoting time-series
function [] = TimeScatter()

end

%Description: Create Figure for Paper
function [] = PaperFigure()
global KM_pair Features Graph_Params params G4

%Create TiledLayout
figure
subplot(2,3,[1 2 4 5])

%Create a 2D plot of 2 shape parameters

index = [];
if (isfile(G4 + "_Kmeans_All.mat"))
    index = load(G4 + "_Kmeans_All.mat",'index');
    index = index.index;
else
    index = kmeans(Features,2,'MaxIter',300,'Replicates',300,'Display','off','Options',statset('UseParallel',1));
    save(G4 + "_Kmeans_All.mat",'index');
end

%Scatter Plot for Clusters 1/2
c1_x = Features(index == 1,Graph_Params(1));
c1_y = Features(index == 1,Graph_Params(2));
c2_x = Features(index == 2,Graph_Params(1));
c2_y = Features(index == 2,Graph_Params(2));

scatter(c1_x,c1_y,16,[0 0.4470 0.7410],'.')
hold on
scatter(c2_x,c2_y,16,[0.8500 0.3250 0.0980],'.')
xlabel(params(Graph_Params(1)))
ylabel(params(Graph_Params(2)))
title(G4 + ": " + params(Graph_Params(1)) + " vs. " + params(Graph_Params(2)))
a = get(gcf,'Children')
a.FontSize = 16;
hold off

%Add the Silhouette score
subplot(2,3,3)
[silh,h] = silhouette(Features,index);
xlabel('Silhouette Value')
ylabel('Cluster')
Score = mean(silh);
title("Silhouette Score: " + num2str(Score))

%Add the Bar Graph Tile
subplot(2,3,6)
pairs = zeros(8,2);
for r = 1:8
    indices = index(r:8:length(index));
    pairs(r,1) = sum(indices==1);
    pairs(r,2) = sum(indices==2);
end
% Create categories (must reorder to make sure order is not changed)
X = categorical({'GG_1','GG_2','GG_3','GG_4','GG_5','GG_6','GG_7','GG_8'});
X = reordercats(X,{'GG_1','GG_2','GG_3','GG_4','GG_5','GG_6','GG_7','GG_8'});
b = bar(X, pairs);
b(1,1).FaceColor = [0 0.4470 0.7410];
b(1,2).FaceColor  = [0.8500 0.3250 0.0980];
title('Clusters per GG-Pair Structure')
ylabel('Counts')
end

%Description: Load files at beginning of file in order to data ready for
%analysis.
function bool = LoadFiles()
global KM_pair seven_ftrs Features Indices TI TI_Length params minpts pairsize pair G4 AI_Ind KM_All G4_workspace

%Ask user if they want to exchange the original output of Curves in order
%to fix it
dat_exc = input('Would you like to exchange columns of negated data (1,2,4,5)?\n This would fix the negated values (Type 0 or 1): ' );
dat_exc = cast(dat_exc,'logical');
G4_workspace = G4;
if (dat_exc)
    G4_workspace = G4 + "_de"; %Change naming for new variables
end

if (~isfile(G4_workspace + "_Workspace.mat"))
    
    % Load Shape Features
    Features = readmatrix(G4 + ".csv");
    if (dat_exc)
        Features(1:8:length(Features),1:2) = Features(1:8:length(Features),1:2).*-1;
        Features(1:8:length(Features),4:5) = Features(1:8:length(Features),4:5).*-1;
        Features(2:8:length(Features),1:2) = Features(2:8:length(Features),1:2).*-1;
        Features(2:8:length(Features),4:5) = Features(2:8:length(Features),4:5).*-1;
    end
    
    dat_l = length(Features);
    TI = readmatrix('1kf1_TI.dat');
    TI_Length = length(TI);
    Indices = load(G4 + "_Indices.mat");
    Indices = Indices.Indices;
    params = ["Shift (dx)" "Slide (dy)" "Rise (dz)" "Tilt (\tau)" "Roll (\rho)" "Twist (\Omega)" "TI Value"];
    minpts = 12;
    pairsize = length(Features)/8;
    pair = 0;
    AI_Ind = [(49*8 + 1):100*8:dat_l (49*8 + 2):100*8:dat_l (49*8 + 3):100*8:dat_l (49*8 + 4):100*8:dat_l (49*8 + 5):100*8:dat_l (49*8 + 6):100*8:dat_l (49*8 + 7):100*8:dat_l (49*8 + 8):100*8:dat_l];
    AI_Ind = sort(AI_Ind);
    
    %TEMPORARILY CHANGE THIS
    seven_ftrs = load("1KF1_7ftrs.mat");
    seven_ftrs  = seven_ftrs.AI;
    seven_ftrs (1:8:length(seven_ftrs ),1:2) = seven_ftrs (1:8:length(seven_ftrs ),1:2).*-1;
    seven_ftrs (1:8:length(seven_ftrs ),4:5) = seven_ftrs (1:8:length(seven_ftrs ),4:5).*-1;
    seven_ftrs (2:8:length(seven_ftrs ),1:2) = seven_ftrs (2:8:length(seven_ftrs ),1:2).*-1;
    seven_ftrs (2:8:length(seven_ftrs ),4:5) = seven_ftrs (2:8:length(seven_ftrs ),4:5).*-1;
    
    % Load K_means clustering for individual GG pairs
    Kmeans_pair_File = G4 + "_KM_pair.mat";
    if (isfile(Kmeans_pair_File))
        KM_pair = load(Kmeans_pair_File);
        KM_pair = KM_pair.KMEANS;
    else
        KM_pair = zeros(pairsize , 8);
        centers = input('How many clusters would you like to create');
        for p = 1:8
            pair_ind = p:8:length(Features);
            KM_pair(:,p) = kmean
            s(Features(pair_ind,:),centers,'MaxIter',300,'Replicates',300,'Display','off','Options',statset('UseParallel',1));
        end
        save(Kmeans_pair_File,'KM_pair');
    end
    
    % Load K_means clustering for all data
    Kmeans_pair_File = G4 + "_KM.mat";
    if (isfile(Kmeans_pair_File))
        KM_All = load(Kmeans_pair_File);
        KM_All = KM_All.idx;
    else
        centers = input('How many clusters would you like to check for all the data? ');
        KM_All = zeros(length(Features) , centers);
        for p = 1:centers
            KM_All(:,p) = kmeans(Features,p,'MaxIter',500,'Replicates',350,'Display','off','Options',statset('UseParallel',1));
        end
        save(Kmeans_pair_File,'KM_All');
    end
    
    %Save for later
    save(G4_workspace + "_Workspace.mat");
    bool = false;
else
    bool = true;
end
end