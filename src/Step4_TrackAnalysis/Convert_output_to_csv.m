

rootFolder = '/mnt/grey/DATA/rawData_2021/Experiment_50_Sara_WT/003_FastTimeLapse_RAMM_Test/'
analyzedFolder = strcat(rootFolder,'Analyzed')
cd(analyzedFolder)

load('Tracking/FullTrack.mat')
tracksFolder = 'Track_IDs'
mkdir(tracksFolder)
cd(tracksFolder)

number_tracks = size(FullTrack);
number_tracks = number_tracks(1)

for i=1:number_tracks 
    
track_size = size(FullTrack(i).CompleteTrack);
track_length = track_size(2);
track = FullTrack(i).CompleteTrack';

% track format: frame number, mask number
file_name = strcat('Track_ID_',num2str(i),'.csv');
csvwrite(file_name,track);
disp(file_name)
end

%% converts Frame files into lists of centroids for each cellID in a Frame file

imagesFolder = '/Tiling_Drift_PostProcess';
framesFolder = strcat(analyzedFolder,imagesFolder)

outputFolder=  strcat(analyzedFolder,'/masks')
mkdir(outputFolder)

myFiles = dir(fullfile(framesFolder,'Frame_*.mat')); 

for k = 1:length(myFiles)
    baseFileName = myFiles(k).name;
    fullFileName = fullfile(framesFolder, baseFileName);
    disp(fullFileName)
    
    frame = load(fullFileName);

    frame_fields= fields(frame);
    for j=1:length(frame_fields)
        if contains(frame_fields{j},'Frame')
            frame_name = frame_fields{j};
        end
    end
    
    stats=regionprops(frame.(frame_name),'Centroid');
    number_masks = size(stats);
    number_masks = number_masks(1)
    output=zeros(3,number_masks);
    
    for i=1:number_masks 
        output(:,i) = [i, stats(i).Centroid];
    end

    output=output';
    
    % output format: maskID, x-coord, y-coord
    file_name = strcat('Masks_',num2str(k),'.csv');
    outputFullFileName= fullfile(outputFolder, file_name);
    csvwrite(outputFullFileName,output);
    disp(outputFullFileName)

end


