function [Cand,Cell] = reverse_parameter_retrieval(Cell_K, Cand_k, ObjectProps, ObjectProps_K)


%% Function to retrieve the parameters of the Cell_K and Cand_k
%% ------------------------------------------------------------

% Create output matrixes that will be input for AHP function
% ----------------------------------------------------------
Cell = [];
Cand = [];


% Retrieve parameters of the Cell_K
% ---------------------------------

% Retrieve Area under cell mask
Area_K = ObjectProps_K.Area(Cell_K);

% Retrieve pixel locations in big Label matrix to define overlap
Overlap_K = cell2mat(ObjectProps_K.PixelIdxList(Cell_K));

% Retrieve the backbone of the cell
ObjBin = cell2mat(ObjectProps_K.Image(Cell_K,:));
% NewCol = zeros(size(ObjBin,1),1);
% ObjBin = horzcat(NewCol, ObjBin, NewCol);
% NewRow = zeros(1, size(ObjBin,2));
% ObjBin = vertcat(NewRow, ObjBin, NewRow);
ObjBin = padarray(ObjBin, [1 1], 0, 'both');
Compl = imcomplement(ObjBin);
DistFromZero = bwdist(Compl);

Max1 = islocalmax(DistFromZero,1);
Max2 = islocalmax(DistFromZero,2);
Overlay = imoverlay(Max1, Max2);
Overlay = imbinarize(Overlay(:,:,1));
Overlay = bwmorph(Overlay,'spur',2);
Length_K = size(find(Overlay==1),1);         % Actual definition of length of the backbone

% Save parameters in output matrix 'Cell'
Cell = horzcat(Cell_K, Length_K, Area_K);


% Retrieve parameters of the Cand_k
% ----------------------------------

Candidates = [];

for c = 1 : size(Cand_k,1)
    C = Cand_k(c);
    
    % Retrieve Area under candidate mask
    Area_k = ObjectProps.Area(C);
    
    % Retrieve pixel locations in big Label matrix to define overlap
    Overlap_k = cell2mat(ObjectProps.PixelIdxList(C));
    Overlap_kK = size(intersect(Overlap_k, Overlap_K),1);
    
    % Retrieve the backbone of the candidate
    ObjBin_k = cell2mat(ObjectProps.Image(C,:));
%     NewCol_k = zeros(size(ObjBin_k,1),1);
%     ObjBin_k = horzcat(NewCol_k, ObjBin_k, NewCol_k);
%     NewRow_k = zeros(1, size(ObjBin_k,2));
%     ObjBin_k = vertcat(NewRow_k, ObjBin_k, NewRow_k);
    ObjBin_k = padarray(ObjBin_k, [1 1], 0, 'both');
    Compl_k = imcomplement(ObjBin_k);
    DistFromZero_k = bwdist(Compl_k);
    
    Max1_k = islocalmax(DistFromZero_k,1);
    Max2_k = islocalmax(DistFromZero_k,2);
    Overlay_k = imoverlay(Max1_k, Max2_k);
    Overlay_k = imbinarize(Overlay_k(:,:,1));
    Overlay_k = bwmorph(Overlay_k,'spur',2);
    Length_k= size(find(Overlay_k==1),1);       % Actual definition of length of the backbone
    
    Candidates = horzcat(C, Length_k, Area_k, Overlap_kK);
    
    % Save parameters in output matrix 'Cand'
    Cand = vertcat(Cand, Candidates);
    Cand = double(Cand);
end