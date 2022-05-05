
%% Main code for Ecoli parameters
%% ==============================

clear all
close all
clc

% Create a new matfile to save E. coli paramters in
Dir_data = uigetdir('/mnt/grey/','Select the folder called Analyzed.');
cd(Dir_data)
cd('Tiling_Drift_PostProcess')

% Count the number of files
NFiles = size(dir('Frame*.mat'),1);

% Loop over the files
parfor k = 1: NFiles
    
    % Get the E.coli files (EC_frame...)
    % Get the Myxo files
    cd(Dir_data)
    cd('Tiling_Drift_PostProcess')
    Name = strcat('Frame_', num2str(k, '%03d'));
    Image = matfile(Name, 'Writable', true);
    Image_postprocess = Image.(Name); % Outputs Frame_001 bwconncomp
    
    Name_EC = strcat('EC_Frame_', num2str(k, '%03d'));
    Image_EC = matfile(Name_EC, 'Writable', true);
    Image_postprocess_EC = Image_EC.(Name_EC); % Outputs EC_Frame_001 bwconncomp
    
    %% 
    % SF 1: Which E. coli segment belong to island/small group/single cell?
    % Input: EC_Frame...
    % Output: Im_EC_final, Image_postprocess_EC.EcoliAssignment -->
    % Both needs to be saved in EC_Frame_...
    % NOTE: Can this info be saved in ConnComp
        
    [Im_EC_final, EcoliAssignment] = EcoliAssignmentIslands(Image_postprocess_EC);
    
    % SF 2: Does a Myxo cell have physical Myxo neigbors?
    % + physical Ecoli neighbors?
    % Are those Ecoli neighbors part of an island/small group/single cells?
    
    [Total_Neighbors, Total_EC_Neighbors, Total_EC_class] = PhysicalNeighborsMyxo(Image_postprocess, Image_postprocess_EC, Im_EC_final);
    
    % SF 3: Voronoi Tesselation based on Ecoli points: Value for Myxo
    % points
    
    CentroidList_EC = Image_EC.CentroidList;
    VoronoiTessellation = Image.VoronoiTessellation;
    CentroidList_Myxo = [VoronoiTessellation(:,4), VoronoiTessellation(:,3)];
    
    [Myxo_Area, a] = EcoliVoronoiAreaForMyxo(CentroidList_EC, CentroidList_Myxo);
    
    % SF 4: Voronoi Tesselation based on Myxo+Ecoli points: Value for Myxo
    % points
    
    VoronoiPoints = Image.VoronoiPoints;
    Values_trackpoints = EcoliMyxoVoronoiAreaForMyxo(VoronoiPoints, CentroidList_EC, CentroidList_Myxo, Image_postprocess);
    
    % SF 5: Voronoi Tesselation based on Ecoli points and 1 Myxo point:
    % Value for 1 Myxo point
    
    % Myxo_Area_OneMyxo = EcoliOneMyxoAreaForMyxo(Image_postprocess, VoronoiTessellation, CentroidList_EC);
    % NOT USED HERE BECAUSE SUPER TIME CONSUMING
    
    % SF 6: Amount of Ecoli biomass in perimeter around Myxo cell
    
    %% Save the variables in a matfile
    % SF1
    Name_ECAssignment = strcat('EC_AssignmentIslands');
    Image_EC.(Name_ECAssignment) = EcoliAssignment;
    
    % SF 2
    Name_matObj = strcat('IntMyxo-EC_Frame_', num2str(k, '%03d'));
    matObject_Int = matfile(Name_matObj, 'writable', true);
    
    Name_Total_Neighbors = strcat('Total_Neighbors');
    Name_Total_EC_Neighbors = strcat('Total_EC_Neighbors');
    Name_Total_EC_class = strcat('Total_EC_class');
    matObject_Int.(Name_Total_Neighbors) = Total_Neighbors;
    matObject_Int.(Name_Total_EC_Neighbors) = Total_EC_Neighbors;
    matObject_Int.(Name_Total_EC_class) = Total_EC_class;
    
    % SF 3
    Name_Myxo_Area = strcat('Myxo_Area_BasedOnECPoints');
    Name_EC_Area = strcat('EC_VoronoiArea');
    matObject_Int.(Name_Myxo_Area) = Myxo_Area;
    Image_EC.(Name_EC_Area) = a;
    
    % SF 4
    Name_MyxoEC_Area = strcat('Myxo_Area_BasedOnECAndMyxoPoints');
    matObject_Int.(Name_MyxoEC_Area) = Values_trackpoints;
    
    disp(strcat('Frame #', num2str(k,'%03d'), ' was saved'))
    
end
