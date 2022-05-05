%% Code with subfunctions for retrieve parameters of the masks
%% ===========================================================

clear
close all
clc

% Load in the matfile
% -------------------

% Define where the data is stored
Dir_data = uigetdir('/mnt/grey/','Select the folder where the matfiles of Tiled Images are saved. (Tiling_Drift_PostProcess)');
cd(Dir_data)
% Count the number of files
NFiles = size(dir('Frame*.mat'),1);
NFiles_overlap = size(dir('Overlap*.mat'),1);

% Loop over the segments
% ----------------------

parfor k = 1:NFiles
% for k = 1:2
    
    Name = strcat('Frame_', num2str(k, '%03d'));
    Image = matfile(Name, 'Writable', true);
    Image_postprocess = Image.(Name);
    
    Name_overlap = strcat('OverlapRegion_frame_', num2str(k, '%03d'));
    Image_overlap = matfile(Name_overlap, 'Writable', false);
    Image_overlap = Image_overlap.(Name_overlap);
    
    Name_Voronoi = strcat('ImageVoronoi_', num2str(k, '%03d'));
    Image_Voronoi = Image.(Name_Voronoi);
    
    % Flag the cells in overlap regions
    % ---------------------------------
    
    Im_tiled = labelmatrix(Image_postprocess);
    Im_overlap = uint16(logical(labelmatrix(Image_overlap)));
    Cells_overlap = Im_tiled.*Im_overlap;
    Cells_overlap = unique(nonzeros(Cells_overlap));
    
    Cells = zeros(Image_postprocess.NumObjects,1);
    Cells(Cells_overlap) = 1;
    
    Image.OverlapRegionFlag = Cells;
    
    % Flag the cells in border regions (10 pixels border)
    % ---------------------------------------------------

    Name_Area_image = strcat('AreaImage_frame_', num2str(k, '%03d'));
    Area_image = Image.(Name_Area_image);
    
    Area_image = logical(labelmatrix(Area_image));
    Image_border = bwmorph(Area_image, 'shrink',10);
    Image_border = uint16(Area_image-Image_border);
    Border_overlap = Im_tiled.*Image_border;
    Cells_borderoverlap = unique(nonzeros(Border_overlap));
    
    Cells_border = zeros(Image_postprocess.NumObjects,1);;
    Cells_border(Cells_borderoverlap) = 1;
    
    Image.BorderCells_10Pixels = Cells_border;
    
    % Reconstruct the backbone lengths
    % --------------------------------
    
    [Lengths] = backbone_calculation(Image_postprocess);
    
    Image.Backbone = Lengths;
    
    % Retrieve the location of the centroids used for Voronoi tessellation
    % & Voronoi tessellation area size
    % --------------------------------------------------------------------
    

    [Centroid_list, Total_Voronoi] = Voronoi_tessellation(Image_Voronoi, Image_postprocess);
    
    Image.VoronoiPoints = Centroid_list;
    Image.VoronoiTessellation = Total_Voronoi;
    
    disp(strcat('Frame #', num2str(k,'%03d'), ' was analyzed'))
end

