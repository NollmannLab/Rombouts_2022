function [Lengths] = backbone_calculation(Image_postprocess)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code: Sara Rombouts (CBS, Team marcelo Nollmann)
%
% Created: 14/09/2020
%
% Goal of code: function to retrieve the backbone of the cells as is done
% in the function 'parameter_retrieval' during tracking
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lengths = [];

% Regionprops to have the isolated image of each segment
% ------------------------------------------------------

Labelmatrix = labelmatrix(Image_postprocess);
ObjectProps = regionprops('table', Labelmatrix, 'image');

% loop over objects
% -----------------

for object = 1 : size(ObjectProps,1)
    
    % Calculation of backbone
    % -----------------------
    
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
    
    Lengths = [Lengths; Length];
end