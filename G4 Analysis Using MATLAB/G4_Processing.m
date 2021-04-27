%Description: This script is meant to aid the analysis of the G4 structures
% and holds various functions to performs different types of analysis,
% particularly engaged with unsupervised machine learning algorithms.
global KMEANS
global DATA params Indices TI TI_Length minpts pair pairsize Graph_Params pairnum seven_ftrs

if(LoadFiles())
    load('G4_Workspace.mat')
end
%% User InputIndices

%Ask for what GG Pair we want to analyze
Pairs_to_analyze = input("What pair(s) would you like to analyze with dbscan? ");
Graph_Params = input("\nWhat parameters would you like to graph? ");

fprintf('You''ve chosen to graph %s and %s\n', params(Graph_Params(1)), params(Graph_Params(2)));

for i = 1:length(Pairs_to_analyze)
    pairnum = Pairs_to_analyze(i);
    pair = DATA(Indices(:,pairnum),:);
    
    %Ask for options per pair
    Input = -2;
    while (Input ~= -1)
        
        
        fprintf('\nTo perform on Pair %i\n',pairnum);
        fprintf('1)K-Dist Graph\n2)GaugeEpsilons\n3)Gauge Epsilons (No files)\n4) Graph TI and Scatter\n');
        fprintf('5)Dendogram \n6)Scatter plot from Hierarchical Cutoff\n7)K-Means clustering with time\n');
        fprintf('8)Rolling Average\n9)Time Scatter Plots');
        Input = input('Choice: ');
        
        switch Input
            case 1
                if (i == 1)
                    PlotK_Dist(false)
                else
                    PlotK_Dist(true)
                end
                legend
                
            case 2
                GaugeEpsilon()
            case 3
                GaugeIndivEp()
            case 4
                PlotTI()
            case 5
                Dendrograms()
            case 6
                ClusterHier()
            case 7
                Kmeans_Time()
            case 8
                RollingAverage()
            case 9
                TimeScatter()
        end
    end
end

%% Additional Functions

% Description: This function tests a range of values between 2 epsilons.
% The epsilon is part the parameter that's changed in the DBSCAN clustering
% method. It creates various files with statistics regarding each
% clustering, and separate files with more detailed information of each
% clustering (corresponding to a unique epsilon).
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

%Description: Plot's the smallest Kth Distance for every point n the data.
function [] = PlotK_Dist(append)
global minpts pair
if (append)
    hold on
else
    figure
    title("K-Distance Graph")
    ylabel("12th Smallest Distances to points")
    xlabel("Number of points ")
end

kD = pdist2(pair,pair,'euc','Smallest',minpts); % The minpts smallest pairwise distances
plot(sort(kD(end,:)));

end

function [] = PlotTI()
global pair minpts pairnum Graph_Params params TI seven_ftrs;

epsilon = input('Please enter an epsilon value: ');

% Scatter Plot with DBSCAN Labels
figure
labels = dbscan(pair,epsilon,minpts);
gscatter(pair(:,Graph_Params(1)),pair(:,Graph_Params(2)),labels);
title(strcat(params(1)," vs. ",params(2)))
xlabel(params(1))
ylabel(params(2))

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

function [] = Dendrograms()
global DATA pair pairnum
figure
link = linkage(pair, 'average');
dendrogram(link);

cutoff = input('What cutoff would you like? ');
dendrogram(link,'ColorThreshold',cutoff);
%Create cutoff and scatter plot
figure
T = cluster(link,'Cutoff',cutoff,'Criterion','distance');
s = scatter(pair(:,3),pair(:,4),20,T,'.');
s.AlphaData  = .01;
end

function [] = ClusterHier()
global pair

figure
T = clusterdata(pair,'Linkage','median','MaxClust',7);
scatter(pair(:,3),pair(:,4),zeros(9970,1),5,T,'filled')

end

%Description: This function performs a K-means clustering with various k's
%with silhouette scores to see if there's any
function [] = Kmeans_Time()
global pair seven_ftrs TI pairnum KMEANS

%Create time column for data
time_1kf1 = [1:2991 3001:3486 3501:5993 6001:10000]' ;

inp = input('1) Cluster Plots\n2) silh\n3) Cluster Evaluation\n4) Plot TI Scatter\n5) K-Means TI Histograms\n\nChoice: ');

% 2 clusters since this maximized the Sillhouette scores before.
clusts = 2;
idx = KMEANS(:,pairnum);
switch inp
    case 1
        ClusterScat()
    case 2
        silh()
    case 3
        output = KMeansTimeEval();
    case 4
        TI_Plot_KMeans()
    case 5
        TI_Hist_KMeans()
end

%K mean analysis from MATLAB
%pool = parpool('threads');

%Description: Simple scatter plot showing the K-means clusters.
    function output = ClusterScat()
        %Go through each cluster
        figure
        for i = 1:clusts
            clusterpts = (idx == i);
            scatter(pair(clusterpts,1),pair(clusterpts,3))
            hold on
        end
        title("K-Means clustering (k=2) - GG_" + num2str(pairnum))
    end

%Description: Sillhouette score for the K-Means clustering
    function [] = silh()
        figure
        [silh,~] = silhouette(pair,idx);
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
                time = time_1kf1(idx == clust);
                
                figure
                plot(1:length(time),time)
                hold on
                time = diff(time);
                plot(1:length(time),time.*15)
                title("K-Means (k=2) of GG_" + num2str(pairnum) + "(cluster " + num2str(clust) + ")")
                
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
%NOTE: NEEDS REVISION due to function IntersectionTI
    function [] = TI_Plot_KMeans()
        for centers = 2:clusts
            for i = 1:centers
                C = IntersectionTI(idx, i);
                figure
                plot(1:length(C),seven_ftrs(C,7));
                title("TI Values in Noise (k=" + int2str(centers) + ")" + ", Cluster " + int2str(i) )
                ylabel("TI Value")
                xlabel("Data Count")
                figure
                scatter(seven_ftrs(C,3),seven_ftrs(C,4));
                title("K-Means Clustering, (k=" + int2str(centers) + ")" + ", Cluster " + int2str(i) )
            end
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

%Description: Load files at beginning of file in order to data ready for
%analysis.
function bool = LoadFiles()
global KMEANS seven_ftrs DATA Indices TI TI_Length params minpts pairsize pair
clear
if (~isfile('G4_Workspace.mat'))
    % Load additional files
    seven_ftrs = load('7_features.mat');
    seven_ftrs  = seven_ftrs .AI;
    DATA = readmatrix("1kf1.csv");
    Indices = load('Indices.mat');
    Indices = Indices.Indices;
    TI = readmatrix('1kf1_TI.dat');
    TI_Length = length(TI);
    params = ["Shift (dx)" "Slide (dy)" "Rise (dz)" "Tilt (\tau)" "Roll (\rho)" "Twist (\Omega)" "TI Value"];
    minpts = 12;
    pairsize = 9970;
    pair = 0;
    
    % Load K_means clustering.
    Kmeans_File = 'KMEANS_IND.mat';
    if (isfile(Kmeans_File))
        KMEANS = load(Kmeans_File);
        KMEANS = KMEANS.idx;
    else
        KMEANS = zeros(pairsize , 8);
        centers = input('How many clusters would you like to create');
        for p = 1:8
            pair_ind = p:8:79760;
            KMEANS(:,p) = kmeans(DATA(pair_ind,:),centers,'MaxIter',300,'Replicates',300,'Display','off','Options',statset('UseParallel',1));
        end
        save(Kmeans_File,'KMEANS');
    end
    
    %Save for later
    save('G4_Workspace.mat');
    bool = false;
else
    bool = true;
end
end