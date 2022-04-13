%% DRIFT CORRECTION AND TILING OF THE 3X3 MOSAIC IMAGE
%% ===================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code: Sara Rombouts (CBS, Team marcelo Nollmann)
%
% 01/09/2020
%
% Goal of code: Main code to handle 3x3 mosaic data acquired by
% RAMM microscope and segmented with Deep Learning network - This code will
% tile the 9 ROIs with a specific architechture into 1 large image, drift
% correct the tiled images over the course of the timelapse and
% post-process the images. Cross correlations for tiling and drift
% correction should be calculated before by python code and will be used
% here as input.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STEP 1: Go to the data
%%

% Select the experiment
% ---------------------
Dir_data = uigetdir('/mnt/grey/','Select the folder where the experiment is saved.');
cd(Dir_data)

% Select the folder where the cross correlations (python - stored in
% txt-files) are stored and retrieve list of all txt-files
% --------------------------------------------------------
Dir_CC = uigetdir(Dir_data, 'Select the folder where the Tiling-CrossCorr are stored.');
% List of cross correlation files
cd(Dir_CC)
List_CC = dir('*.txt');

% Select the folder where we can find the drift correction txt-file and
% load in the data in this file
% -----------------------------
cd(uigetdir(Dir_data, 'Select the folder where the Drift CrossCorr are stored.'))
fileID = fopen('XYshift_DriftCorr_ROI5.txt', 'r');
formatSpec = '%f';
Drift = fscanf(fileID,formatSpec);
% Reformat the drift correction file
oddI = 1:2:size(Drift,1);
evenI = 2:2:size(Drift,1);
Drift = [Drift(oddI),Drift(evenI)];

% Select which images should be stitched and retrieve a list of file names
% ------------------------------------------------------------------------
list = {'BF_normalized','Ecoli_normalized','Ecoli_segmented','Myxo_segmented', 'RGB'};
[indx,tf] = listdlg('ListString',list, 'SelectionMode','single');
Data = list{indx};
% Create list of file names
cd(strcat(Dir_data,'/Segmented_images/ROI_1/',Data));
List_names = dir('*.tif');

% Create matfile to store final bwconnected objects in
% ----------------------------------------------------
cd(Dir_data)
mkdir('Analyzed')
cd('Analyzed')
mkdir('Tiling_Drift_PostProcess')
cd('Tiling_Drift_PostProcess')

% matObject_Tiled = matfile('ConnComp_Tiled_DriftCorr_PostProcess_226-end.mat', 'writable', true);
% matObject_Overlap = matfile('ConnComp_OverlapRegion_226-end.mat', 'writable', true);

%% STEP 2: open loop over all time points
%%
tic
for t = 1:size(List_names,1)
%     for t = 1:5
    
    % Tiling
    % ------
    
    Name = List_names(t).name;
    Name_CC = strcat(Dir_CC,'/',List_CC(t).name);
    
    [Image, Overlap, Area_image] = MosaicTiling(Dir_data, Name, Name_CC, Data);
    
    
    % Drift correction of tilematObjd image
    % -------------------------------
    
    [New_Image, New_Overlap, New_Area_image] = DriftCorr(Image, Overlap, Area_image, Drift,t);
    
    % Post-processing of the Drift-corrected tiled image
    % --------------------------------------------------
    
    switch indx
        
        case 4
            % function for Myxo postprocessing
            [Conn, M] = PostProcess_TiledIm(New_Image);
            
            
            % Retrieve the Overlap region and save in a matfile
            % -------------------------------------------------
            New_Overlap = bwconncomp(logical(New_Overlap),8);
            
            Name_Overlap = strcat('OverlapRegion_frame_', num2str(t, '%03d'));
            matObject_Overlap = matfile(Name_Overlap, 'writable', true);
            matObject_Overlap.(Name_Overlap) = New_Overlap;
            disp(strcat('Overlap frame #', num2str(t,'%03d'), ' was saved'))
            
            % Save the total image area mask
            New_Area_image = bwconncomp(logical(New_Area_image),8);
            
            % Save list of all masks (bwconncomp object) and their locations to
            % reconstruct total image so that you do not
            Name_Tiled = strcat('Frame_', num2str(t, '%03d'));
            Name_TiledVoronoi = strcat('ImageVoronoi_', num2str(t, '%03d'));
            Name_Area_image = strcat('AreaImage_frame_', num2str(t, '%03d'));
            
            matObject_Tiled = matfile(Name_Tiled, 'writable', true);
            matObject_Tiled.(Name_Tiled) = Conn;
            matObject_Tiled.(Name_TiledVoronoi) = M;
            matObject_Tiled.(Name_Area_image) = New_Area_image;
            
            disp(strcat('Frame #', num2str(t,'%03d'), ' was saved'))
            
            %         Lx = 5880;
            %         Ly = 5880;
            %         cd(Dir_data)
            % %         cd('Analyzed')
            %     %     im_name = strcat(num2str(t,'%03d'), '_Myxo_drift_corrected.tif');
            %     im_name = strcat(num2str(t,'%03d'), '_Myxo.tif');
            %
            %         T = Tiff(im_name, 'w');
            %         tagstruct = struct('ImageLength', Lx, ...
            %             'ImageWidth', Ly, ...
            %             'BitsPerSample', 16, ...
            %             'Photometric', Tiff.Photometric.MinIsBlack, ...
            %             'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky, ...
            %             'Compression', Tiff.Compression.None, ...
            %             'SamplesPerPixel',1);
            %         T.setTag(tagstruct);
            %         T.write(uint16(Image));
            %         close(T)
            
            
        case 3
            % function for E. coli postprocessing
            [Centroid_list, Ec_conn4] = Postprocessing_Ecoli(New_Image, t, Dir_data);
            
            
            Name_FinalEc = strcat('EC_Frame_', num2str(t, '%03d'));
            Name_EC_CentroidList = strcat('CentroidList');
            
            matObject_EC = matfile(Name_FinalEc, 'writable', true);
            matObject_EC.(Name_FinalEc) = Ec_conn4;
            matObject_EC.(Name_EC_CentroidList) = Centroid_list;
            
            disp(strcat('EC Frame #', num2str(t,'%03d'), ' was saved'))
            
        otherwise
            disp('No Post-processing is done for this data type.')
    end
end
toc

    