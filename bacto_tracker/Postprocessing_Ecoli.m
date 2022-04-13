function [Centroid_list, Ec_conn4] = Postprocessing_Ecoli(New_Image, t, Dir_data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code: Sara Rombouts (CBS, Team marcelo Nollmann)
%
% Created: 16/11/2020
%
% Goal of code: function to post-process E. coli images and results in the
% E. coli masks and for each mask the centroid
%
% NOTE: This code runs for the TILED AND POST-PROCESSED IMAGES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: convert output image after tiling/drift correction to binary
% image (setting al nonzero values to 1)
Bin_Ec = logical(New_Image);

% Step 2: Set all pixels which are defined as a Myxo cell body in the
% binary e.coli image to zero (under the assumtion that 1 pixel can only be
% attributed to 1 species)

% Load in the Postprocessed Myxo image
cd(Dir_data)
cd('Analyzed')
cd('Tiling_Drift_PostProcess')

Name = strcat('Frame_', num2str(t, '%03d'));
Image = matfile(Name, 'Writable', true);
Image_Myxo = Image.(Name);
Bin_Myxo = logical(labelmatrix(Image_Myxo));

% Set all Myxo pixels to zero
Bin_Ec = Bin_Ec-Bin_Myxo;
Bin_Ec(Bin_Ec~=1)=0;

% Step 3: Get all 4-connected objects in the binary E.coli image and filter
% them for size (filter out anything smaller then, and including, 10 pixels)

Filter_Ec = bwareaopen(double(Bin_Ec), 11,4);
Ec_conn4 = bwconncomp(Filter_Ec, 4);

% Step 4: Get the centroids of the segments (we do this in the same way as
% we do for the Myxo cells - instead of getting the center of the bounding
% box, we take the center of the maskitself!)

Lx = size(Bin_Ec,1);
Ly = size(Bin_Ec,2);

Compl = imcomplement(logical(labelmatrix(Ec_conn4)));
DistFromZero = bwdist(Compl);
DistFromZero(DistFromZero==1) = 0;
Final = bwmorph(DistFromZero, 'Thicken', 1);
Final = bwpropfilt(Final, 'area', [10 1000000],4);
Final_conn = bwconncomp(Final,4);
Final_labelmatrix = labelmatrix(Final_conn);

Im_thin = bwmorph(logical(Final_labelmatrix), 'thin', Inf);
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
Centroid_list = table2array(regionprops('table', logical(Centroid), 'Centroid'));

% Output of this function = Centroid_list, connected component of
% 'post-processed' E. coli image Ec_conn4


