function [Total_Neighbors, Total_EC_Neighbors, Total_EC_class] = PhysicalNeighborsMyxo(Image_postprocess, Image_postprocess_EC, Im_EC_final)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code: Sara Rombouts (CBS, Team marcelo Nollmann)
%
% Created: 1/12/2020
%
% Goal of code: function to define the physical neighbors of Myxo cells
% (both Ecoli and Myxo) and to which group (island/small group/single cell)
% the E.coli neighbors belong)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Im_Frame = labelmatrix(Image_postprocess);
Im_EC = labelmatrix(Image_postprocess_EC);

Total_Neighbors = cell(Image_postprocess.NumObjects,1,1);
Total_EC_Neighbors = [];
Total_EC_class = cell(Image_postprocess.NumObjects,1,1);

ObjectProps = regionprops('table', Im_Frame, 'BoundingBox');

for i = 1:Image_postprocess.NumObjects

% Parameter 1: Physical neighbors of Myxo cells
    BB = ObjectProps.BoundingBox(i,:);
    I2 = imcrop(Im_Frame,[BB(:,1)-30, BB(:,2)-30, BB(:,3)+60, BB(:,4)+60]);
    I3 = I2;
    
    I3(I3~=i)=0;
    I3 = bwmorph(I3, 'thicken', 7);
    
    Neighbors = nonzeros(unique(I2(find(I3))));
    Neighbors(Neighbors==i)=[];
    
    if ~isempty(Neighbors)
        Total_Neighbors{i}=Neighbors;
    else
        Total_Neighbors{i}=0;
    end
    
    
% Parameter 2: Physical contact between E. coli and Myxo cells
% Parameter 3: Do Myxo cells contact islands/small groups/single cells?
    
    EC_I2 = imcrop(Im_EC,[BB(:,1)-30, BB(:,2)-30, BB(:,3)+60, BB(:,4)+60]);
    EC_final_crop = imcrop(Im_EC_final,[BB(:,1)-30, BB(:,2)-30, BB(:,3)+60, BB(:,4)+60]);
    
    EC_Neighbors = nonzeros(unique(EC_I2(find(I3))));
    EC_Neighbors_Class = nonzeros(unique(EC_final_crop(find(I3))));
    
    if ~isempty(EC_Neighbors)
        Total_EC_Neighbors=[Total_EC_Neighbors;1];
        Total_EC_class{i} = EC_Neighbors_Class;
    else
        Total_EC_Neighbors=[Total_EC_Neighbors;0];
        Total_EC_class{i} = 0;
    end
    
end