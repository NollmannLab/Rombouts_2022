function [Cand,Cell] = parameter_retrieval(object, ObjectProps, Labelmatrix_K, ObjectProps_K)

%% Function to select the candidates and retrieve their parameters
%% ---------------------------------------------------------------

% Create output matrixes that will be input for AHP function
% ----------------------------------------------------------

Cell = [];
Cand = [];

% Define the possible candidates in next frame for the track
% ----------------------------------------------------------

% Based on the bounding box
% Retrieve the bonding box of frame k and frame K (k+1)
BB = ObjectProps.BoundingBox(object,:);
BB_k = [BB(1,1)-5 BB(1,2)-5 BB(1,3)+10 BB(1,4)+10];

% Find all objects in BB_k region in next frame K
Crop = imcrop(Labelmatrix_K, BB_k);
Candidates = unique(nonzeros(Crop));

% Labelmatrix_object = LabelMatrix;
% Labelmatrix_object(Labelmatrix_object~=object) = 0;

if isempty(Candidates)
    Cand = [];
else
    Cand_K = [];
    
    
    % Retrieve parameters of the CELL TO BE TRACKED on which we will base the tracking
    % --------------------------------------------------------------------------------
    
    % Retrieve Area under cell mask
    Area = ObjectProps.Area(object);
    
    % Retrieve pixel locations in big Label matrix to define overlap
    Overlap = cell2mat(ObjectProps.PixelIdxList(object));
    
    % Retrieve the backbone of the cell
    ObjBin = cell2mat(ObjectProps.Image(object,:));
    ObjBin = padarray(ObjBin, [1 1], 0, 'both');
    Compl = imcomplement(ObjBin);
    DistFromZero = bwdist(Compl);
    
    Max1 = islocalmax(DistFromZero,1);
    Max2 = islocalmax(DistFromZero,2);
    Overlay = imoverlay(Max1, Max2);
    Overlay = imbinarize(Overlay(:,:,1));
    Overlay = bwmorph(Overlay,'spur',2);
    Length = size(find(Overlay==1),1);         % Actual definition of length of the backbone
    
    % Save parameters in output matrix 'Cell'
    Cell = horzcat(object, Length, Area);
    
    % Retrieve parameters of the CANDIDATES on which we will base the tracking
    % ------------------------------------------------------------------------
    
    
    for c = 1 : size(Candidates,1)
        C = Candidates(c);
        
        % Retrieve Area under candidate mask
        Area_K = ObjectProps_K.Area(C);
        
        % Retrieve pixel locations in big Label matrix to define overlap
        Overlap_K = cell2mat(ObjectProps_K.PixelIdxList(C));
        Overlap_kK = size(intersect(Overlap, Overlap_K),1);
        
        % Retrieve the backbone of the candidate
        ObjBin_K = cell2mat(ObjectProps_K.Image(C,:));
        ObjBin_K = padarray(ObjBin_K, [1 1], 0, 'both');
        Compl_K = imcomplement(ObjBin_K);
        DistFromZero_K = bwdist(Compl_K);
        
        Max1_K = islocalmax(DistFromZero_K,1);
        Max2_K = islocalmax(DistFromZero_K,2);
        Overlay_K = imoverlay(Max1_K, Max2_K);
        Overlay_K = imbinarize(Overlay_K(:,:,1));
        Overlay_K = bwmorph(Overlay_K,'spur',2);
        Length_K = size(find(Overlay_K==1),1);       % Actual definition of length of the backbone
        
        Cand_K = horzcat(C, Length_K, Area_K, Overlap_kK);
        Cand = vertcat(Cand, Cand_K);
        Cand = double(Cand);
    end
    
end

        
        
        
      

            