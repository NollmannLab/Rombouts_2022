%% Filering of the tracks of final tracks (>10 consecutive frames)
%% ===============================================================

%% Step 1: For what do we want to filter? 
%% --------------------------------------

% Filter for size: yes (=1) or no (=0), and for which size?
Size = 1;
Frames = 35;    % Minimal track length filter

% Filter for tracks where cell travels through an overlap region? If yes (=1) or no (=0)
Overlap = 1;

% Filter for cells at the cell border? If yes (=1) or no (=0)
Border = 1;

%% Step 2: Filtering
%% -----------------

% Option 1: Filtering for size
% ----------------------------

Lengths = [];
for i = 1: size(Total_tracks,1)
   A = size(nonzeros(Total_tracks(i,:)),1);
   Lengths = [Lengths;A];
    
end


Filter_Size = Lengths>=Frames;
Filter_Size = imcomplement(Filter_Size);
Filter_Size = Size*Filter_Size;

% Option 2: Filtering for tracks that travel through the overlap region
% ---------------------------------------------------------------------

Filter_Overlap = sum(Total_Flag,2);
Filter_Overlap = logical(Filter_Overlap);
Filter_Overlap = Overlap*Filter_Overlap;

% Option 3: Filtering out of track where segments are on the image border
% -----------------------------------------------------------------------

Filter_Border = sum(Total_Border,2);
Filter_Border = logical(Filter_Border);
Filter_Border = Border*Filter_Border;

% Tracks to keep are the ones labeled with a zero
% -----------------------------------------------

Total = Filter_Size+Filter_Overlap+Filter_Border;

Total_tracks = Total_tracks(find(Total==0),:);
Total_Density = Total_Density(find(Total==0),:);
Track_Densities = Track_Densities(find(Total==0),:);
Total_Row = Total_Row(find(Total==0),:);
Total_Col = Total_Col(find(Total==0),:);








