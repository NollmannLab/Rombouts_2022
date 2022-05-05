function Myxo_Area_OneMyxo = EcoliOneMyxoAreaForMyxo(Image_postprocess, VoronoiTessellation, CentroidList_EC)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code: Sara Rombouts (CBS, Team marcelo Nollmann)
%
% Created: 1/12/2020
%
% Goal of code: function to get voronoi value for Myxo point (Myxo cell) of
% interest in E.coli points (Voronoi tesselation based on 1 Myxo point and
% all E.coli points in cropped region)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Im_Frame = labelmatrix(Image_postprocess);
ObjectProps = regionprops('table', Im_Frame, 'BoundingBox');

a = [];
A = [];
for j = 1:size(VoronoiTessellation,1)
    
PointOfInterest = [VoronoiTessellation(j,4), VoronoiTessellation(j,3)];

Centroids = [PointOfInterest; CentroidList_EC];

PixelList = sub2ind([Image_postprocess.ImageSize(1) Image_postprocess.ImageSize(2)], Centroids(:,2), Centroids(:,1));

Myxo_Test_Im = zeros(Image_postprocess.ImageSize(1), Image_postprocess.ImageSize(2));
Myxo_Test_Im(PixelList(1))=1;

Test_Im = zeros(Image_postprocess.ImageSize(1), Image_postprocess.ImageSize(2));
Test_Im(PixelList)=1;

% Get Bounding box of Myxo cell of interest and crop the image (for speed)
BB = ObjectProps.BoundingBox(j,:);
Im_crop_points = imcrop(Test_Im,[BB(:,1)-400, BB(:,2)-400, BB(:,3)+800, BB(:,4)+800]);
Im_crop_Myxo = imcrop(Myxo_Test_Im,[BB(:,1)-400, BB(:,2)-400, BB(:,3)+800, BB(:,4)+800]);

C = table2array(regionprops('table', logical(Im_crop_points), 'Centroid'));
C_Myxo = table2array(regionprops('table', logical(Im_crop_Myxo), 'Centroid'));

Myxo_loc = find(ismember(C,C_Myxo,'rows'));

[v,c] = voronoin(C);
  
    if size(find(c{Myxo_loc}==1),2)==0
        
        x = v(c{Myxo_loc},1);
        y = v(c{Myxo_loc},2);
        A = polyarea(x,y);
       
        a = vertcat(a,A);


    else
        A = NaN(1);
        a = vertcat(a,A);

end


end

Myxo_Area_OneMyxo = a;
