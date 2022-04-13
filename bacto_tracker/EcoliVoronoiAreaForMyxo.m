function [Myxo_Area, a] = EcoliVoronoiAreaForMyxo(CentroidList_EC, CentroidList_Myxo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code: Sara Rombouts (CBS, Team marcelo Nollmann)
%
% Created: 1/12/2020
%
% Goal of code: function to get the Voronoi Tesselation based on the 
% centroids of Ecoli cells and finally the corresponding Voronoi Area for 
% the Myxo centroid points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 1: Get Voronoi tesselation of E. coli image and its corresponding
% area values
% Correct also for NaN values by taking the average of 3 closest points
% (representing the three closest E. coli cells (or at least their
% centroid))

[v,c] = voronoin(CentroidList_EC);

a = [];
Poly = [];

for i = 1:length(c)
    
    b = [];
    if size(find(c{i}==1),2)==0
        
        x = v(c{i},1);
        y = v(c{i},2);
        A = polyarea(x,y);

        a = vertcat(a,A);

    else
        A = NaN(1);
        a = vertcat(a,A);
    end
end


% note: INSTEAD OF GIVING THE CELLS AN ARBITRARY NUMBER AS
% DENSITY, CALCULATE THE AVERAGE DESITY OF THE 3 NEIREST
% NEIGHBORS AND USE THIS AS THE DENSITY VALUE.

[row, ~] = find(isnan(a));

for j = 1 : length(row)
    J = row(j);
    b = [];
    % Compute the euclidian distances
    dist = sum((CentroidList_EC-CentroidList_EC(J,:)).^2,2);
    dist(row)=0;
    [~,r] = mink(dist,sum(dist==0)+3);
    
    % Calculate the average density (polyarea) of 3 nearest
    % neighbors
    for k = sum(dist==0)+1 : length(r)
        K = r(k);
        x = v(c{K},1);
        y = v(c{K},2);
        A = polyarea(x,y);
        
        b = vertcat(b,A);
        
    end
    
    Avg = mean(b,1);
    a(J,:) = Avg;
end


% STEP 2: Get the cellIndex for each Myxo cell (this cell index represents
% which E. coli cell is located the closest to the Myxo point of interest)

[~, cellIndex] = min(pdist2(CentroidList_EC, CentroidList_Myxo)); 
cellIndex = cellIndex';

% STEP 3: Get the corresponding voronoi area (based on E. coli points) for
% each Myxo point

Myxo_Area = a(cellIndex);

% STEP 4: Save or outputresult


