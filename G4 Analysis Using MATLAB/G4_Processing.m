global DATA params Indices TI TI_Length minpts pair pairsize Graph_Params pairnum
DATA = readmatrix("1kf1.csv");
Indices = load('Indices.mat');
Indices = Indices.Indices; %Load as a matrix and not as a structure
TI = readmatrix('1kf1_TI.dat');
TI_Length = length(TI);
params = ["Shift (dx)" "Slide (dy)" "Rise (dz)" "Tilt (\tau)" "Roll (\rho)" "Twist (\Omega)"];
minpts = 12;
pairsize = 9970;
pair = 0;

%% User Input

%Ask for what GG Pair we want to analyze
Pairs_to_analyze = input("What pair(s) would you like to analyze with dbscan? ");
ranges = input("Please tells us the ranges for Epsilon in an n x 2 matrix: ");
Graph_Params = input("\nWhat parameters would you like to graph? ");

fprintf('You''ve chosen to graph %s and %s\n', params(Graph_Params(1)), params(Graph_Params(2)));

for i = 1:length(Pairs_to_analyze)
    pairnum = Pairs_to_analyze(i);
    pair = DATA(Indices(:,pairnum),:);
    
    fprintf('\nTo perform on Pair %i\n',pairnum);
    
    %Ask for options per pair
    Input = input('What to do?\n1)K-Dist Graph\n2)GaugeEpsilons\n3)Gauge Epsilons (No files)\n4) Graph TI and Scatter\n\n');
    
    switch Input
        %K-Dist Graph(s)
        case 1
            if (i == 1)
                PlotK_Dist(false);
            else
                PlotK_Dist(true);
            end
            legend
            
            %Gauge Various Epsilons
        case 2
            GaugeEpsilon(ranges(i,:));
            
            %Graph Scatter Plot
        case 3
            GaugeIndivEp(ranges(i,:));
            
            % Graph scattering plotting and 
        case 4
            Input = input('Please enter an epsilon value: ');
            PlotTI(Input)
    end
    
end


%% Additional Functions
function [] = GaugeEpsilon(eps)
global pair minpts pairsize TI pairnum

%Allocate memory
iter = 50;

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

function [] = GaugeIndivEp(eps)
global pair minpts pairsize TI pairnum

%Allocate memory
iter = 50;

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

function C = IntersectionTI(labels,clusnum)
% Description: This function takes in labels and gives back the indices for
% intersection of values within a cluster and the TI subset of data

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

%Convert DBIndices to OG indices
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

function [] = PlotK_Dist(append)
global minpts pair pairnum
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

function [] = PlotTI(epsilon)
global pair minpts pairnum Graph_Params params TI;
% Scatter Plot with DBSCAN Labels
figure
labels = dbscan(pair,epsilon,minpts);
gscatter(pair(:,Graph_Params(1)),pair(:,Graph_Params(2)),labels);
title(strcat(params(1)," vs. ",params(2)))
xlabel(params(1))
ylabel(params(2))

% % Plot the TI values within cluster of choice
% C = IntersectionTI(labels, 1);
% figure
% plot(1:length(C),TI(C));
% title("TI Values in Noise")
% ylabel("TI Value")
% xlabel("Data Count")
% 
% % Plot OG TI points
% figure
% PairAb = (49*8 + pairnum):100*8:79760;
% PairAb  = ((PairAb - pairnum - 392)/100) + pairnum;
% plot(1:length(PairAb),TI(PairAb));
% title("TI Values (Original Values)")
% ylabel("TI Value")
% xlabel("Data Count")
end