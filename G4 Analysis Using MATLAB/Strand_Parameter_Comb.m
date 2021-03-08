%Ab_Initio data
%GGPairs = SeparateGGPairs('1kf1.csv');

Message = "Please enter the option you'd like:\n1)Create Indices\n2)GraphPlanes\n";
choice = input(Message);

%% Graph the different parameter 2D combinations with 1 plot per strand
switch (choice)
    case 1
        %Read File
        FILE = input("What's your file name? ",'s');
        PARAMS = readmatrix(FILE);
        length = size(PARAMS);
        %Allocate
        Indices = zeros(length(1)/8,8);
        
        %Fill each Column
        for i = 1:8
            Indices(:,i) = (i:8:length(1))';
        end
        
        %For the next 2, find the ab-initio values for each pair
        
        
    case 2
        
end

%function [] = Graph
params = {'Shift (dx)' 'Slide (dy)' 'Rise (dz)' 'Tilt (tau)' 'Roll (rho)' 'Twist (Omega)'};
paramnum = 6 - 1;
strandname = {'2-3-4 Strand' '8-9-10 Strand' '14-15-16 Strand' '20-21-22 Strand'};
legends = {'g2g3' 'g3g4' ; 'g8g9' 'g9g10' ; 'g14g15' 'g15g16' ; 'g20g21' 'g21g22' };
Planes = [3 5 7 ; 4 6 8 ];
for i = 1:paramnum
    for j = i+1: 6
        figure
        % set(gcf,'Visible','off')
        tiledlayout(1,2)
        for k = 1:2
            %Set Up
            nexttile
            
            %Draw plane
            Pair = Planes(k,:);
            a = scatter(GGPairs(:,i,Pair(1)),GGPairs(:,j,Pair(1)),16,'.');
            a.AlphaData = 0.5;
            hold on
            b = scatter(GGPairs(:,i,Pair(2)),GGPairs(:,j,Pair(2)),16,'.');
            b.AlphaData = 0.5;
            hold on
            c = scatter(GGPairs(:,i,Pair(3)),GGPairs(:,j,Pair(3)),16,'.');
            c.AlphaData = 0.5;
            %Labels
            legend(legends{k,1},legends{k,2})
            title(strandname{k});
            %Labeling
            file = strcat(params{i}," vs. ",params{j});
            title(file)
            xlabel(params{i})
            ylabel(params{j})
        end
    end
end

%Create figure
figure
% set(gcf,'Visible','off')
tiledlayout(2,2)
for k = 1:4
    %Set Up
    nexttile
    strand = strands(k);
    c = linspace(1,10,length(GGPairs(:,1,strand)));
    
    %Graphing info
    l = size(GGPairs(:,i,strand));
    scatter(GGPairs(:,i,strand),GGPairs(:,j,strand),1, c)
    hold on
    scatter(GGPairs(:,i,strand+1),GGPairs(:,j,strand+1),1, c)
    
    %Labels
    legend(legends{k,1},legends{k,2})
    title(strandname{k});
    %Labeling
    file = strcat(params{i},' vs.',params{j});
    title(file)
    xlabel(params{i})
    ylabel(params{j})
    colorbar
end
%Save figures
file = strcat(file, '.jpg');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% exportgraphics(gcf,file, 'Resolution', 300);
% close(gcf)