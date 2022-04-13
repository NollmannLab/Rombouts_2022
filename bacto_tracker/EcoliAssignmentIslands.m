function [Im_EC_final, EcoliAssignment] = EcoliAssignmentIslands(Image_postprocess_EC)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code: Sara Rombouts (CBS, Team marcelo Nollmann)
%
% Created: 1/12/2020
%
% Goal of code: function to define whether the Ecoli segments belong to
% either a island or small group or are single cells
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Im_EC = labelmatrix(Image_postprocess_EC);
Im_EC_bin = logical(Im_EC);
Im_EC_dist = bwdist(Im_EC_bin);

Im_EC2 = Im_EC_dist;
Im_EC2(Im_EC2<=10)=0;
Im_EC2_bin = logical(Im_EC2);
Im_EC2_comp = imcomplement(Im_EC2_bin);

Im_EC2_conn = bwconncomp(Im_EC2_comp,8);
Im_EC2_LM = labelmatrix(Im_EC2_conn);

Im_EC_zero = zeros(Image_postprocess_EC.ImageSize(1), Image_postprocess_EC.ImageSize(2));
% Get for each segment whether it is island, pack of E. coli or single cell
for i = 1:Im_EC2_conn.NumObjects
    Pixels = find(Im_EC2_LM==i);
    B = nonzeros(unique(Im_EC(Pixels)));
    
    if size(B,1) == 0 %E. coli is not found in this inflated island (will probably never happen)
        Im_EC_zero(Pixels) = 1;
        
    elseif size(B,1)==1 % E. coli cell is single cell
        Im_EC_zero(Pixels) = 2;
        
    elseif size(B,1)~=0 & size(B,1)<=20 % E. coli cells form a small group
        Im_EC_zero(Pixels) = 3;
        
    elseif size(B,1)>=21 % E. coli islands form a true island
        Im_EC_zero(Pixels) = 4;
    end
    
end

Im_EC_final = Im_EC_bin.*Im_EC_zero;


% Get for each 4-connected segment its assignment to island/small
% group/single cell
Test = [];
for j = 1:Image_postprocess_EC.NumObjects
   Pixels = cell2mat(Image_postprocess_EC.PixelIdxList(j));
   Assign = unique(Im_EC_final(Pixels));
    
    Test = [Test, Assign];
end

EcoliAssignment = Test;

