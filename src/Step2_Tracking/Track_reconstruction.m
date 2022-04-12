function [Full_Track_2] = Track_reconstruction(Filter_Track)

%% RECONSTRUCTION OF TRACKS
%% ========================

% For every frame we loop over CellIDs and reconstruct the tracks
% We are looping over the structs for every frame in k
% But we want to loop also over all CellIDs and match them to TrackCell
% of previous frames
% We do this by creating a new array in which every row represents a
% track, and every column represents the CellID in the following frame
% (chronological order of CellIDs)
% We use the Sure_Tracks_Frame_k structs, instead of the Raw and
% unfiltered results per frame

% Start reconstructing the full tracks
% Load in variables from matfile
% Reconstruct full tracks - save both the tracks (sequence of Cell ID that
% represent the track of 1 specific cell
% Write full tracks and additional information about the track away in a
% cell (Cell = Full_Track)

Full_Track={};

for j = 1 : size(Filter_Track,1)
    
    Name = strcat('PW_Tracks_Filter_', num2str(j, '%03d'));
    Frame = matfile(Name, 'Writable', false);
    Frame = Frame.(Name);
    TCID = Frame(:,1);
    TC = Frame(:,2);
    Score = Frame(:,3);
    Score(Score==0)=1000;
    
    if j==1
        FullTrack = [TCID,TC];
        FullScore = Score;
        % So we have the original CellID from frame 1 in column 1 and the
        % TrackCell in column 2
    else
        % Loop over every Cell ID in frame k, find the value in the k-th
        % column (because the TrackCell from the previous frame corresponds
        % to the CellID in the next frame
        for l = 1 : size(TCID,1)
            L = TCID(l);
            NextCell = TC(l);
            [row_CellID,~] = find(FullTrack(:,j)==L);       % Find in which row the CellID can be found
            NumTracks = size(FullTrack, 1);                 % We define the amount of rows in the FullTrack matrix to know on which row a new track  should be started
            S = Score(l);
            
            if ~isempty(row_CellID)                         % If the CellID is detected, then we complete the matrix with the TrackCell in the next column
                FullTrack(row_CellID, j+1) = NextCell;
                FullScore(row_CellID, j)= S;
                
            else isempty(row_CellID);                       % This means that a new track is started
                New_row = NumTracks+1;                      % Select the row in which the new track needs to start
                FullTrack(New_row, j) = L;                  % Start new track be copying the CellID in the FullTrack matrix in column corresponding to frame
                FullTrack(New_row, j+1) = NextCell;
                FullScore(New_row, j) = S;
                
            end
        end
    end
    
end

FullScore = [FullScore,zeros(size(FullScore,1),1)];

% Save all tracks from the huge matrix

for p = 1 : size(FullTrack,1)
    
    Track_temp = FullTrack(p,:);
    
    Time = find(Track_temp~=0);
    Start_frame = Time(1,1);
    End_frame = Time(1,end);
    
    if End_frame - Start_frame >10
        Timestamp = {Start_frame, End_frame};
        Scores = FullScore(p,Start_frame:End_frame);
        
        Complete_track = cat(1,Time,Track_temp(Track_temp~=0));
        Total_Track = {[Complete_track], Timestamp, Scores};
        Full_Track = [Full_Track; [Total_Track]];
    end
end

% Save the cell that is combining all the tracks

Col_headers = {'CompleteTrack' 'Timestamp_StartAndEnd' 'Scores'};
Full_Track_2 = cell2struct(Full_Track, Col_headers, 2);