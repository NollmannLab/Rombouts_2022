

%% This program is used to calculate the segmented image of Sara's time-lapse.
%% It works in four steps :
%%    1- Load the trained networks
%%    2- Indicate the folder where the Brightfield images are saved. The
%%    folder should already contain as many folders as there was ROIs imaged
%%    during the experiment. The same architecure will be saved.
%%    3- Point to the folder containing the deconvolved images of E.coli and
%%    calculate for each image the average image (based on the std)
%%    4- Calculate the segmented images and return as output the two binary
%%    images for Myxo and Ecoli.
%%
%% 18-06-2020 : Add the possibility to load a third channel for Myxo fluo. 
%% In that case, the images should already be calculated for the in-focus.
%% ========================================================================

clear
close all
clc

% Ask whether the workspace should be cleaned or we can keep the previous
% data (the trained net already loaded)
% -------------------------------------

Starting_path = '/mnt/PALM_dataserv/DATA/JB/Sara_phd/Deep_Learning/';
Trained_network_folder = uigetdir(Starting_path, 'Select the folder where the trained networks are saved');
cd(Trained_network_folder)

Network = dir('*Trained_network_*.mat');
Nnetwork = size(Network,1);

Net = cell(Nnetwork,1);

for n = 1 : Nnetwork
    load(Network(n).name);
    Net{n} = lgraph;
end

% Select the folder where the brightfield im--cleanAllMasksages are located
% ----------------------------------------------------------

Starting_path = '/mnt/grey/DATA/';
Im_Path_BF = uigetdir(Starting_path, 'Select the BRIGHTFIELD folder');
Folders_BF = dir(strcat(Im_Path_BF, '/ROI_*'));

if ~isempty(Folders_BF)
    
    Nfolders = 0;
    Folders_list = [];
    
    for n = 1 : size(Folders_BF,1)
        if Folders_BF(n).isdir == 1
            Nfolders = Nfolders + 1;
            Folders_list = cat(1, Folders_list, n);
        end
    end
else
    hwarn = warndlg('The folder was expected to contain folders called "ROI_xxx". None were found and this operation is aborted');
    uiwait(hwarn)
    delete(hwarn)
    return
end

% Indicate the folder containing all the deconvolved images of E.coli and 
% from infer the number of frames composing the stacks as well as the image
% dimensions. 
% -----------

Starting_path = '/mnt/grey/DATA/';
Im_Path_Ecoli = uigetdir(Starting_path, 'Select the ECOLI folder where the DECOLVOLVED images are saved');
Tiff_content = dir(strcat(Im_Path_Ecoli, '/*.tif'));

if ~isempty(Tiff_content)
    im_info = imfinfo(strcat(Im_Path_Ecoli, '/', Tiff_content(1).name));
    Nframes = size(im_info,1);
    Lx = im_info(1).Height;
    Ly = im_info(1).Width;
    Ecoli_stack = zeros(Lx,Ly,Nframes);
else
    hwarn = warndlg('The folder was expected to contain tif images. None were found and this operation is aborted');
    uiwait(hwarn)
    delete(hwarn)
    return
end

% Ask whether we want to load a third channel for myxo. If yes indicate the
% folder where the images are saved
% ---------------------------------

Third_channel = questdlg('Do you want do use fluorescence data from the Myxo channel?', 'Third channel?', 'Yes', 'No', 'Yes');

switch Third_channel
    case 'Yes'
        Im_Path_Myxo = uigetdir(Starting_path, 'Select the MYXO folder where the DECOLVOLVED & IN-FOCUS images are saved');
        Tiff_content = dir(strcat(Im_Path_Myxo, '/*.tif'));
        
        if isempty(Tiff_content)
            hwarn = warndlg('The folder was expected to contain tif images. None were found and this operation is aborted');
            uiwait(hwarn)
            delete(hwarn)
            return
        end
end

% Indicate the folder where the segmented data should be saved
% -------------------------------------------------------------

Starting_path = '/mnt/grey/';
Saving_path = uigetdir(Starting_path, 'Select the folder where the data should be saved (make sur you have right permission!)');

Saving_path = strcat(Saving_path, '/Segmented_images/');
mkdir(Saving_path);

% For each image the segmented network is applied
% -----------------------------------------------

classNames = ["Background", "Myxo", "Myxo_contour", "Ecoli", "Ecoli_contour"];
switch Third_channel
    case 'No'
        Lz = 2;
    case 'Yes'
        Lz = 3;
end
im_standardized = zeros(Lx,Ly,Lz);
        
% for nfolder = 1 : Nfolders
for nfolder = 1 : Nfolders
    
    % From the brightfield folder, calculate the number of images to
    % analyze
    % -------
    
    Current_path_BF = strcat(Im_Path_BF, '/', Folders_BF(Folders_list(nfolder)).name);
    Folder_content = dir(strcat(Current_path_BF, '/*.tif'));
    Nim = size(Folder_content,1);
    
    % Create the folders where the segmented images will be saved
    % ----------------------------------------------------------
    
    CurrentSavingPath = strcat(Saving_path, Folders_BF(Folders_list(nfolder)).name);
    mkdir(CurrentSavingPath)
    
    cd(CurrentSavingPath)
    mkdir('BF_normalized')
    mkdir('Ecoli_normalized')
    mkdir('Ecoli_segmented')
    mkdir('Myxo_segmented')
    mkdir('RGB')
    
    switch Third_channel
        case 'Yes'
            mkdir('Myxo_normalized')
    end
    
    for nim = 1 : Nim
        
        % Load the brightfield image
        % --------------------------
        
        cd(Current_path_BF)
        im_name = Folder_content(nim).name;
        im(:,:,1) = imread(im_name);
        
        % Load the deconvolved stack for Ecoli
        % ------------------------------------
        
        cd(Im_Path_Ecoli)
        Template_name = strcat(im_name(1:3), '*', im_name(end-8:end-4), '*.tif');
        Tiff_content = dir(Template_name);
        if size(Tiff_content,1)==1
            im_name_ecoli = Tiff_content(1).name;
            for nframe = 1:Nframes
                Ecoli_stack(:,:,nframe) = imread(im_name_ecoli, 'index', nframe);
            end
            
        else
            fprintf(strcat('The following Ecoli stack was not found : ', ...
                im_name, '\n'));
            continue
        end
        
        % Calculate the average image of the E.coli stack
        % -----------------------------------------------
        
        im(:,:,2) = std(Ecoli_stack,0,3);
        
        % In case load the myxo in-focus image
        % ------------------------------------
        
        switch Third_channel
            case 'Yes'
                
                cd(Im_Path_Myxo)
                Template_name = strcat(im_name(1:3), '*', im_name(end-8:end-4), '*.tif');
                Tiff_content = dir(Template_name);
                if size(Tiff_content,1)==1
                    im_name_myxo = Tiff_content(1).name;
                    im(:,:,3) = imread(im_name_myxo);
                else
                    fprintf(strcat('The following Myxo in-focus image was not found : ', ...
                        im_name, '\n'));
                    continue
                end
        end
                

        % Define the arrays for each class
        % ---------------------------------
        
        Myxo = zeros(Lx, Ly, Nnetwork);
        Ecoli = zeros(Lx, Ly, Nnetwork);
        
        % Normalize the image by using a Gaussian filter
        % ----------------------------------------------
        
        try
            for n_layer = 1 : Lz
                [im_standardized(:,:,n_layer), im(:,:,n_layer)] = Normalization(im(:,:,n_layer), 3);
            end
        catch
            fprintf('Image #%i in ROI #%i was not analyzed due to an error in the normalization process \n',nim, nfolder)
        end
            
        % Save the normalized image
        % -------------------------
        
        cd(strcat(CurrentSavingPath, '/BF_normalized'))
        BF_normalized_name = strcat(num2str(nim, '%.3i'), '_normalized_BF.tif');
        imwrite(im(:,:,1), BF_normalized_name, 'tif', 'Compression', 'none')
        
        cd(strcat(CurrentSavingPath, '/Ecoli_normalized'))
        Ecoli_normalized_name = strcat(num2str(nim, '%.3i'), '_normalized_Ecoli.tif');
        imwrite(im(:,:,2), Ecoli_normalized_name, 'tif', 'Compression', 'none')
        
        switch Third_channel
            case 'Yes'
                cd(strcat(CurrentSavingPath, '/Myxo_normalized'))
                Myxo_normalized_name = strcat(num2str(nim, '%.3i'), '_normalized_Myxo.tif');
                imwrite(im(:,:,3), Myxo_normalized_name, 'tif', 'Compression', 'none')
        end
        
        % Run the CNN on the image. For the new version of the software, each
        % images are normalized before applying the network. Note the images
        % are normalized in such way that the intensity lies between 0 and 255.
        % -------------------------------------------------------------------
        
        for n = 1 : Nnetwork
            
            [Reconstructed_im, ~] = semanticseg(im_standardized, Net{n});
            Myxo(:,:,n) = im2uint8(Reconstructed_im == classNames(2));
            Ecoli(:,:,n) = im2uint8(Reconstructed_im == classNames(4));
        end
               
        av_Myxo = double(255*mean(Myxo/255,3));
        av_Ecoli = double(255*mean(Ecoli/255,3));
        
        RGB = zeros(Lx, Ly, 3);
        RGB(:,:,1) = 0.9952*av_Myxo + 0.138*av_Ecoli;
        RGB(:,:,2) = 0.7893*av_Myxo + 0.6276*av_Ecoli;
        RGB(:,:,3) = 0.2028*av_Myxo + 0.8973*av_Ecoli;
        
        % Display the results
        % -------------------
        
        cd(strcat(CurrentSavingPath, '/Ecoli_segmented'))
        Ecoli_name = strcat(num2str(nim, '%.3i'), '_segmented_Ecoli.tif');
        imwrite(uint8(av_Ecoli), Ecoli_name, 'tif', 'Compression', 'none')
        
        cd(strcat(CurrentSavingPath, '/Myxo_segmented'))
        Myxo_name = strcat(num2str(nim, '%.3i'), '_segmented_Myxo.tif');
        imwrite(uint8(av_Myxo), Myxo_name, 'tif', 'Compression', 'none')
        
        cd(strcat(CurrentSavingPath, '/RGB'))
        RGB_name = strcat(num2str(nim, '%.3i'), '_segmented_RGB.tif');
        imwrite(uint8(RGB), RGB_name, 'tif', 'Compression', 'none')
        
        fprintf('Image #%i in ROI #%i was analyzed \n',nim, nfolder)
    end
end
