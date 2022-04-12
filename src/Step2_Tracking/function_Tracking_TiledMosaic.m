%% TRACKING OF 3x3 MOSAIC IMAGES
%% =============================

close all
clear
clc

%% Part 1: Frame-to-frame AHP tracking
%% -----------------------------------

% Define where the data is stored
Dir_data = uigetdir('/mnt/grey/','Select the folder where the matfiles of Tiled Images are saved.');
cd(Dir_data)
% Count the number of files
NFiles = size(dir('Frame*.mat'),1);

% Define where the pairwise tracks will be saved
cd('..')
mkdir('Tracking')
cd('Tracking')
Dir_save = pwd;
cd(Dir_data)

parfor k = 1 : NFiles-1

    
    % Load in first frame for pairwise tracking by loading in the matfile
    % -------------------------------------------------------------------

    Name = strcat('Frame_', num2str(k, '%03d'));
    Image = matfile(Name, 'Writable', false);
    Image = Image.(Name);
    
    Labelmatrix = labelmatrix(Image);
    ObjectProps = regionprops('table', Labelmatrix, 'Area', 'BoundingBox', 'image', 'perimeter', 'PixelIdxList'); 
    
    % Load in second frame for pairwaise tracking
    % -------------------------------------------
    K = k+1;
    
    Name_K = strcat('Frame_', num2str(K, '%03d'));
    Image_K = matfile(Name_K, 'Writable', false);
    Image_K = Image_K.(Name_K);

    Labelmatrix_K = labelmatrix(Image_K);
    ObjectProps_K = regionprops('table', Labelmatrix_K, 'Area', 'BoundingBox', 'image', 'perimeter', 'PixelIdxList');
    
    for object = 1 : size(ObjectProps,1)
        
        %CREATE FUNCTION FOR THE RETRIEVAL OF PARAMETERS
        % input: object, ObjectProps en LabelMatrixTotal_track
        % output: matrix of Cell IDs of candidates in column 1, cell lengths in
        % column 2, cell area in column 3, overlap in column 4
        % output 2: matrix of Cell ID, cell length and area of cell in frame t
        
        [Cand,Cell] = parameter_retrieval(object, ObjectProps, Labelmatrix_K, ObjectProps_K);
        
        % AHP for determining the Rank of candidates
        % CREATE FUNCTION FOR RANKING THE CANDIDATES BASED ON THE PARAMETERS
        % input : matrix cell and matrix Cand
        % Output: matrix Rank (gives all cell IDs and corresponding Ranks
        
        [Rank] = AHP_approach(Cand,Cell);
        
        % Select the candidate with highest probability based on AHP to be next
        % cell in the track
        
        if isempty(Rank)
            Cell_Track = 0;
            Score = 0;
        elseif size(Rank,1)==1
            Cell_Track = Rank(:,1);
            Score = 1;
        elseif size(Rank,1)>1 && size(find(Rank(:,2)==max(Rank(:,2))),1)==1
            Cell_Track = Rank(find(Rank(:,2)==max(Rank(:,2))),1);
            Score = 2;
        else % size(Rank,1)>1 && size(find(Rank(:,2)==max(Rank(:,2))),1)~=1
            Cell_Track = 0;
            Score = 10;
        end
        
        Track = [object, Cell_Track, Score];
        
        if object==1
            Total_track = Track;
        else
            Total_track = cat(1, Total_track, Track);
        end
    end
        
        Results{k} = Total_track;
        
        %% Part 2: Filtering of the Total_Track
        %% ------------------------------------
        
        Total_filt_track = Filtering_PWtracks(Total_track, ObjectProps, ObjectProps_K);
    
    %% Save the results to a matfile
    Name_PW = strcat(Dir_save, '/PW_Tracks_', num2str(k, '%03d'));
    Fieldname_PW = strcat('PW_Tracks_', num2str(k, '%03d'));
    matObject_PW = matfile(Name_PW, 'Writable', true);
    matObject_PW.(Fieldname_PW) = Total_track;
        
    Name_Filter = strcat(Dir_save, '/PW_Tracks_Filter_', num2str(k, '%03d'));
    Fieldname_Filter = strcat('PW_Tracks_Filter_', num2str(k, '%03d'));
    matObject_Filter = matfile(Name_Filter, 'writable', true);
    matObject_Filter.(Fieldname_Filter) = Total_filt_track;
        
    disp(strcat('Frame #', num2str(k,'%03d'), ' was analyzed and saved'))
end

%% Part 3: Track reconstruction
%% ----------------------------

cd(Dir_save)
Filter_Track = dir('PW_Tracks_Filter*.mat');

% 
% [Full_Track_2, FullScore, FullTrack] = Track_reconstruction(Filter_Track);
[Full_Track_2] = Track_reconstruction(Filter_Track);

Name_FullTrack = strcat(Dir_save, '/FullTrack');
Fieldname_FullTrack = strcat('FullTrack');
FullTrack = matfile(Name_FullTrack, 'Writable', true);
FullTrack.(Fieldname_FullTrack) = Full_Track_2;

disp(strcat('Track reconstruction was analyzed and saved'))
















