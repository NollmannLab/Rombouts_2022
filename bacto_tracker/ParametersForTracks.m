%% Couple the tracks to the parameters of the segments
%% ===================================================


% Load in the matfile of the tracks
% ---------------------------------

% Define the folder where the Analyzed data is saved
Dir_data = uigetdir('/mnt/grey/','Select the folder where the analyzed data is saved.');

% Define where the data is stored
Dir_tracks = strcat(Dir_data, '/Tracking');
cd(Dir_tracks)

Track = matfile('FullTrack.mat', 'Writable', false);
Track = Track.FullTrack;

% Get the path to the frames with the data
% ----------------------------------------
Dir_parameter = strcat(Dir_data, '/Tiling_Drift_PostProcess');
cd(Dir_parameter)

% Count the number of files
NFiles = size(dir('Frame*.mat'),1);

% Reconstruct the total tracks
% ----------------------------

Total_tracks = zeros(size(Track,1), NFiles);
Track_lengths = [];

for i = 1:size(Track,1)
% for i = 1:100
    Begin = cell2mat(Track(i).Timestamp_StartAndEnd(1));
    End = cell2mat(Track(i).Timestamp_StartAndEnd(2));
    
    T = Track(i).CompleteTrack(2,:);
%     Length = size(Track(i).CompleteTrack(2,:),2);
    
    Total_tracks(i,Begin:End) = T;
%     Track_lengths = [Track_lengths; Length];
end

% Filtering on size
% -----------------
% Size = 35;
% Total_tracks = Total_tracks(Track_lengths>=Size,:);


% Retrieve the data for each frame
% --------------------------------
cd(Dir_parameter)

Total_Density = zeros(size(Total_tracks,1), size(Total_tracks,2));
Total_Area = zeros(size(Total_tracks,1), size(Total_tracks,2));
Total_Backbone = zeros(size(Total_tracks,1), size(Total_tracks,2));
Total_Flag = zeros(size(Total_tracks,1), size(Total_tracks,2));
Total_Border = zeros(size(Total_tracks,1), size(Total_tracks,2));

% NEW 20210111
Total_IntEc = zeros(size(Total_tracks,1), size(Total_tracks,2));
Total_ClassEC = zeros(size(Total_tracks,1), size(Total_tracks,2));
Total_MyxoNeighbors = zeros(size(Total_tracks,1), size(Total_tracks,2));
Total_SizeMyxoNeighbors = zeros(size(Total_tracks,1), size(Total_tracks,2));
Total_DensityEc = zeros(size(Total_tracks,1), size(Total_tracks,2));
Total_DensityEcMyxo = zeros(size(Total_tracks,1), size(Total_tracks,2));


for i = 1:size(Total_tracks,2)
    CellID = nonzeros(Total_tracks(:,i));
    Loc_CellID = find(Total_tracks(:,i));
    
    Name = strcat('Frame_', num2str(i, '%03d'));
    Frame = matfile(Name, 'Writable', true);
    
    % Voronoi densities
    Density = Frame.VoronoiTessellation(:,2);
    D = Density(CellID);
    
    % Flag
    Flag = Frame.OverlapRegionFlag;
    F = Flag(CellID);
    
    % Border Cells
    Border = Frame.BorderCells_10Pixels;
    Bo = Border(CellID);
    
    % Backbone
    Backbone = Frame.Backbone;
    B = Backbone(CellID);
    
    % Row Voronoi point
    Row = Frame.VoronoiTessellation(:,3);
    R = Row(CellID);
    
    % Column Voronoi point
    Col = Frame.VoronoiTessellation(:,4);
    C = Col(CellID);
    
    % Area
    Area = Frame.(Name);
    
    A = [];
    for j = 1:size(CellID,1)
        J = CellID(j);
        Area2 = size(cell2mat(Area.PixelIdxList(J)),1);
        A = [A; Area2];
    end
    
    Total_Density(Loc_CellID,i) = D;
    Total_Area(Loc_CellID,i) = A;
    Total_Backbone(Loc_CellID,i) = B;
    Total_Flag(Loc_CellID,i) = F;
    Total_Border(Loc_CellID,i) = Bo;
    Total_Col(Loc_CellID,i) = C;
    Total_Row(Loc_CellID,i) = R;
    
    % NEW 20210111
    Name_EcParameters = strcat('IntMyxo-EC_Frame_', num2str(i, '%03d'));
    Frame_EcParameters = matfile(Name_EcParameters, 'Writable', true);
    
    % Int Ec
    IntEc = Frame_EcParameters.Total_EC_Neighbors;
    IntEc = IntEc+1;
    
    Final_IntEc = IntEc(CellID);
    
    % ClassEc
    ClassEC = Frame_EcParameters.Total_EC_class;
    ClassEC = cellfun(@sum, ClassEC);
    ClassEC(ClassEC==0)=1;
   
    Final_ClassEC = ClassEC(CellID);
    
    
    % MyxoNeighbors
    MyxoNeighbors_1 = Frame_EcParameters.Total_Neighbors;

    MyxoNeighbors = cellfun(@find, MyxoNeighbors_1, 'UniformOutput', false);
    MyxoNeighbors = imcomplement(cellfun(@isempty, MyxoNeighbors));
    MyxoNeighbors = double(MyxoNeighbors)+1;
    
    Final_MyxoNeighbors = double(MyxoNeighbors(CellID));
    
    % SizeMyxoNeighbors
    SizeMyxoNeighbors = cellfun(@find, MyxoNeighbors_1, 'UniformOutput', false);
    SizeMyxoNeighbors = cell2mat(cellfun(@size, SizeMyxoNeighbors, 'UniformOutput', false));
    SizeMyxoNeighbors = SizeMyxoNeighbors(:,1);
    SizeMyxoNeighbors(SizeMyxoNeighbors==0)=100;
    
    Final_SizeMyxoNeighbors = SizeMyxoNeighbors(CellID);
    
    % DensityEc
    DensityEc = Frame_EcParameters.Myxo_Area_BasedOnECPoints;
    
    Final_DensityEC = DensityEc(CellID);
    
    % DensityEcMyxo
    DensityEcMyxo = Frame_EcParameters.Myxo_Area_BasedOnECAndMyxoPoints;
    
    Final_DensityEcMyxo = DensityEcMyxo(CellID);
    
    Total_IntEc(Loc_CellID,i) = Final_IntEc;
    Total_ClassEC(Loc_CellID,i) = Final_ClassEC;
    Total_MyxoNeighbors(Loc_CellID,i) = Final_MyxoNeighbors;
    Total_SizeMyxoNeighbors(Loc_CellID,i) = Final_SizeMyxoNeighbors;
    Total_DensityEc(Loc_CellID,i) = Final_DensityEC;
    Total_DensityEcMyxo(Loc_CellID,i) = Final_DensityEcMyxo;
    
    disp(strcat('Frame #', num2str(i,'%03d'), ' was processed'))
    
end


Track_Densities = Total_Density;
Track_Densities = log10(Track_Densities);
Track_Densities(Track_Densities==-Inf)=0;
Track_Densities(Track_Densities>=3.5)=6;
Track_Densities(Track_Densities<=2.8 & Track_Densities~=0)=2;
Track_Densities(Track_Densities<=3.5 & Track_Densities>=2.8) = 4;

