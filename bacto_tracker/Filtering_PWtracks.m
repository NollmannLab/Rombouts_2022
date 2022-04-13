function Total_filt_track = Filtering_PWtracks(Total_track, ObjectProps, ObjectProps_K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code: Sara Rombouts (CBS, Team marcelo Nollmann)
%
% 11/09/2020
%
% Goal of code: function to filter the original pairwaise tracks
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Filtering of the Total_Track

Filtered_Result = Total_track;

UniqTC = unique(Total_track(:,2));
Counts = histc(Total_track(:,2),UniqTC);
Double_tracks = UniqTC(1<Counts);

for i = 2 : size(Double_tracks,1)
    Cell_K = Double_tracks(i);
    
    % Remove all lines where cellID of frame t1 is chosen >1 in column2
    Filtered_Result(Filtered_Result(:,2)==Cell_K,:) = [];
    
    % Define for each CellID from frame t1 that was chosen >1 the
    % cellIDs from frame t0 for which it was chosen >1 (these CellIDs
    % will be the Candidates for the ReverseAHP
    I = find(Total_track(:,2)==Cell_K);
    
    Cand_k = [];
    for c = 1 : size(I,1)
        C = I(c,:);
        Cand_k = vertcat(Cand_k, Total_track(C,1));
    end
    
    % PARAMETERS RETRIEVAL FOR CELLF ROM FRAME t1 AND CANDIDATES FROM
    % FRAME t0
    [Cand,Cell] = reverse_parameter_retrieval(Cell_K, Cand_k, ObjectProps, ObjectProps_K);
    
    % REVERSE AHP APPROACH
    [Rank] = AHP_approach(Cand,Cell);
    
    % Select the candidate with highest probability based on AHP to be next
    % cell in the track
    
    if isempty(Rank)
        Cell_Track = 0;
        Score = 0;
        %             New_Track = [Cell_Track, Cell_K, Score];
    elseif size(Rank,1)==1
        Cell_Track = Rank(:,1);
        Score = 100;
        %             New_Track = [Cell_Track, Cell_K, Score];
    elseif size(Rank,1)>1 && size(find(Rank(:,2)==max(Rank(:,2))),1)==1
        Cell_K_A = Cell_K;
        Cell_Track_A = Rank(find(Rank(:,2)==max(Rank(:,2))),1);
        Score_A = 200;
        %             New_Track = [Cell_Track, Cell_K, Score];
        Cell_Track_B = Rank(find(Rank(:,2)~=max(Rank(:,2))),1);
        Cell_K_B = zeros(size(Cell_Track_B,1),1);
        Score_B = repmat(500, size(Cell_Track_B,1), 1);
        
        Cell_Track = [Cell_Track_A; Cell_Track_B];
        Cell_K = [Cell_K_A;Cell_K_B];
        Score = [Score_A;Score_B];
        
    else
        Cell_Track = Rank(:,1);
        Score = repmat(800, size(Rank,1),1);
        Cell_K= zeros(size(Rank,1),1);
    end
    
    New_Track = [Cell_Track, Cell_K, Score];
    
    if i==2
        Total_new_track = New_Track;
    else
        Total_new_track = cat(1, Total_new_track, New_Track);
    end
    
    
end

Total_filt_track = cat(1, Filtered_Result, Total_new_track);