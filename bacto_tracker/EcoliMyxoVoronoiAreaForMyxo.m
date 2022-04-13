function Values_trackpoints = EcoliMyxoVoronoiAreaForMyxo(VoronoiPoints, CentroidList_EC, CentroidList_Myxo, Image_postprocess)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code: Sara Rombouts (CBS, Team marcelo Nollmann)
%
% Created: 1/12/2020
%
% Goal of code: function to get the Voronoi Tesselation based on the 
% centroids of Ecoli and Myxo cells and finally the corresponding Voronoi 
% Area for the Myxo centroid points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input function: VoronoiPoints (Frame_...) and CentroidList
% (EC_Frame_...)
% output function: Values_trackpoints

% STEP1: retrieve all points (Myxo and E.coli)

% Adjusted  on 25th March, 2021
[C,ia,ib] = intersect(VoronoiPoints,CentroidList_EC,'rows');
CentroidList_EC(ib,:)=[];

Centroids_all = [VoronoiPoints; CentroidList_EC];

% STEP 2: Get the area size of the tesselated cells
[v,c] = voronoin(Centroids_all);

a = [];

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
    dist = sum((Centroids_all-Centroids_all(J,:)).^2,2);
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

% STEP 3: Retrieve from all area values the ones corresponding to the myxo
% points
ind_allpoints = sub2ind([Image_postprocess.ImageSize(1) Image_postprocess.ImageSize(2)],Centroids_all(:,1),Centroids_all(:,2));
ind_trackpoints = sub2ind([Image_postprocess.ImageSize(1) Image_postprocess.ImageSize(2)],CentroidList_Myxo(:,1), CentroidList_Myxo(:,2));

[~,loc]=ismember(ind_trackpoints,ind_allpoints);
Values_trackpoints = a(loc);

% STEP 4: save or output the results