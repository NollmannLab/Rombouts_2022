function [Centroid_list, Total_Voronoi] = Voronoi_tessellation(Image_Voronoi, Image_postprocess)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code: Sara Rombouts (CBS, Team marcelo Nollmann)
%
% Created: 14/09/2020
%
% Goal of code: function to retrieve the centroids of the masks to
% calculate the area of the Voronoi tesselation for each segment
%
% NOTE: This code runs for the TILED AND POST-PROCESSED IMAGES
%
% NOTE: The voronoi tesseltation is slightly adjusted as compared to the
% last code - Here, we do not take the center of mass of the bounding box
% (regionprops - centroids) but we shrink each mask to a single point -
% WHY? We do this to make sure that each point for the voronoi tesselation
% is located inside the mask (e.g.: for curved cells)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 1: Get the centroids of the segments for tracking + BP corrected segments 
%% ==============================================================================

Lx = size(Image_Voronoi,1);
Ly = size(Image_Voronoi,2);


% SEGMENTS TRACKING: Create 8-connected objects from 4-connected objects by removing the
% perimeter of each 4-connected segment
% --------------------------------------------------------------------------------------
Label_postprocess = labelmatrix(Image_postprocess);
Perimeter = bwperim(logical(Label_postprocess),4);
Cell_body = logical(labelmatrix(Image_postprocess))-Perimeter;
Centroid = bwmorph(logical(Cell_body), 'shrink','Inf');
Centroid_test = find(Centroid);

% Check that there is only 1 point retained per segment!
P = Label_postprocess(Centroid_test);
Binc = [1:Image_postprocess.NumObjects];
Counts = hist(P,Binc);
Result = [Binc; Counts]';

% Correction for all segments for which there is more then 1 centroid
% defined
Total_C = [];
for m = 1:size(Result,1)
   if Result(m,2)>1

       P(P==Result(m,1))=0;
       L = Label_postprocess;
       L(L~=Result(m,1))=0;
       L = bwmorph(L, 'shrink', Inf);
       C = find(L);
       Total_C = [Total_C; C];
       
   end
end

Total = Centroid_test(find(P));
Total = [Total; Total_C];

% SEGMENTS WITH BRANCHPONTS AND TORTUOSITY: Get other segments of Myxo in
% the image and get the centroid for those segments
% -----------------------------------------------------------------------

Compl = imcomplement(logical(Image_Voronoi));
DistFromZero = bwdist(Compl);
DistFromZero(DistFromZero==1) = 0;
Final = bwmorph(DistFromZero, 'Thicken', 1);
Final = bwpropfilt(Final, 'area', [100 1000000],4);
Final_conn = bwconncomp(Final,4);
Final_labelmatrix = labelmatrix(Final_conn);

Other_segments = (logical(Final_labelmatrix)~=logical(labelmatrix(Image_postprocess)));

Im_thin = bwmorph(Other_segments, 'thin', Inf);
Im_spur = bwmorph(Im_thin, 'spur', 3);
Im_branch = bwmorph(Im_spur, 'branchpoints');

% Retrieve the PixelIdxList of all branchpoints
BP = find(Im_branch);
[all_row, all_col] = ind2sub([Lx Ly], BP);

% Get all neighboring pixels
Row = [];
Col = [];

for t = 1:size(BP,1)
    r = all_row(t);
    c = all_col(t);
    
    Neighbors = [r-1,c-1;...
        r-1,c;...
        r-1,c+1;...
        r,c-1;...
        r,c;...
        r,c+1;...
        r+1,c-1;...
        r+1,c;...
        r+1,c+1];
    Neighbors(Neighbors<=0)=1;
    Neighbors(Neighbors>Lx)=Lx;
    
    Row = [Row;Neighbors(:,1)];
    Col = [Col;Neighbors(:,2)];
    
end

ind = sub2ind([Lx Ly],Row,Col);
Im_spur(ind) = 0;

%Fill up all circular segments
Im_spur = imfill(Im_spur,4,'holes');

% Get the centroids of resting 
Centroid = bwmorph(logical(Im_spur), 'shrink', Inf);
Total = [Total; find(Centroid)];
Zero = zeros(Lx,Ly);
Zero(Total) = 1;

%% Step 2: Calculate Voronoi
%% =========================

% Calculation of the centroids of the image containing the centroids
% ------------------------------------------------------------------
Centroid_list = table2array(regionprops('table', logical(Zero), 'Centroid'));
% PixelIdxList = regionprops(Im_spur, 'PixelIdxList');

%%

% Calculate the voronoi diagram
% -----------------------------
[v,c] = voronoin(Centroid_list);

% Create image from voronoi diagram

a = [];
Poly = [];

for i = 1:length(c)
    
    b = [];
    if size(find(c{i}==1),2)==0
        
        x = v(c{i},1);
        y = v(c{i},2);
        A = polyarea(x,y);
        a = vertcat(a,A);
        patch(x,y,A);

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
    dist = sum((Centroid_list-Centroid_list(J,:)).^2,2);
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


%% Step 3: Save the densities for all the segments that are included in tracking
%% =============================================================================

Cells = [];
Voronoi = [];
Points = [];
Label_postprocess = double(labelmatrix(Image_postprocess));
for p = 1:size(a,1)
    Cell = Label_postprocess(Centroid_list(p,2),(Centroid_list(p,1)));
    if Cell~=0
        Cells = [Cells; Cell];
        Voronoi = [Voronoi; a(p)];
        Points = [Points;Centroid_list(p,2),Centroid_list(p,1)];
    end
end
Total_Voronoi = [Cells,Voronoi, Points];

[~,idx] = sort(Total_Voronoi(:,1));
Total_Voronoi = Total_Voronoi(idx,:);

    