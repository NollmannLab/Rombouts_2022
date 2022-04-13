function [Image, Overlap, Area_image] = MosaicTiling(Dir_data, Name, Name_CC, Data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code: Sara Rombouts (CBS, Team marcelo Nollmann)
%
% Created: 31/08/2020
%
% Corrected: 07/09/2020 - Include tiling based on Myxo_segmented image if
% CC on BF_normalized does not give a result (line 30-38)
%
% Goal of code: function to tile the images (RAMM microscope) to mosaic of
% 3x3. Overlap regions merged by Alpha decompositing.
%
% Adjusted 20210215: Tiling CC correction function implemented
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Cross-correlations file
% Open the txt-file and store in a variable
fileID = fopen(Name_CC, 'r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);

% A is the raw text file containing only 1 column with 24 rows instead of 2 cols and 12
% rows - untangle the cols and rows
oddI = 1:2:size(A,1);
evenI = 2:2:size(A,1);
A = [A(oddI),A(evenI)];

if size(A,1)==24
    % CORRECTION ADDED ON 07/09/2020 (line 30-38)
    % A contains all CC (calculated with Python) based on BF_normalized (row 1
    % to 12) and where the CC based on BF_normalized fails, it cotains the CC
    % based on Myxo_segmetend (row 13 to 24 - contains zeros if CC could be
    % calculated based on BF_normalized)
    % BUT we need to get all CCs for all images (does not matter based on which
    % images it is calculated) - create consensus matric of CC by adding the 2
    % together - new CC matrix contains 12 rows with 12 CCs in x and y
    A = A(1:12,:)+A(13:24,:);
    
end

%% If necessary, Tiling CC correction
A = CorrectionTiling(A);

%% Reconstruct image

% Load in the images
ROI1 = imread(strcat(Dir_data,'/Segmented_images/ROI_1/',Data,'/',Name));
ROI2 = imread(strcat(Dir_data,'/Segmented_images/ROI_2/',Data,'/',Name));
ROI3 = imread(strcat(Dir_data,'/Segmented_images/ROI_3/',Data,'/',Name));
ROI4 = imread(strcat(Dir_data,'/Segmented_images/ROI_4/',Data,'/',Name));
ROI5 = imread(strcat(Dir_data,'/Segmented_images/ROI_5/',Data,'/',Name));
ROI6 = imread(strcat(Dir_data,'/Segmented_images/ROI_6/',Data,'/',Name));
ROI7 = imread(strcat(Dir_data,'/Segmented_images/ROI_7/',Data,'/',Name));
ROI8 = imread(strcat(Dir_data,'/Segmented_images/ROI_8/',Data,'/',Name));
ROI9 = imread(strcat(Dir_data,'/Segmented_images/ROI_9/',Data,'/',Name));


% Create Alpha image
R = ones(2048,2048);
Alpha = zeros(2050);
Alpha(2:2+2047,2:2+2047) = R;
Alpha = bwdist(imcomplement(Alpha));
Alpha = Alpha(2:2+2047,2:2+2047);

% Create an empty array
Image = zeros(5850,5850,9);
Image2 = zeros(5850,5850,9);


% Place first ROI5 in the middle of the new image at position 1901,1901
r_begin_5 = 1901;
r_end_5 = 1901+2047;
c_begin_5 = 1901;
c_end_5 = 1901+2047;

Image(1901:1901+2047,1901:1901+2047,5) = ROI5;

Image2(1901:1901+2047,1901:1901+2047,5) = Alpha;

% Place ROI4 on top of ROI5 - need CC5_4 (CC4_5 = row 3)

Image(r_begin_5-2048-A(3,1):r_end_5-2048-A(3,1),c_begin_5-A(3,2):c_end_5-A(3,2),4) = ROI4;

Image2(r_begin_5-2048-A(3,1):r_end_5-2048-A(3,1),c_begin_5-A(3,2):c_end_5-A(3,2),4) = Alpha;

% Place ROI6 on bottom of ROI5 - need CC5_6 (CC5_6 = row 4)

Image(r_end_5+A(4,1):r_end_5+A(4,1)+2047, c_begin_5+A(4,2):c_end_5+A(4,2),6) = ROI6;

Image2(r_end_5+A(4,1):r_end_5+A(4,1)+2047, c_begin_5+A(4,2):c_end_5+A(4,2),6) = Alpha;

% Place ROI8 at left of ROI5 - need CC5_8 (CC8_5 = row 9)

Image(r_begin_5-A(9,1):r_end_5-A(9,1), c_begin_5-2048-A(9,2):c_end_5-2048-A(9,2),8) = ROI8;

Image2(r_begin_5-A(9,1):r_end_5-A(9,1), c_begin_5-2048-A(9,2):c_end_5-2048-A(9,2),8) = Alpha;

% Place ROI2 at right of ROI5 - need CC5_2 (CC5_2 = 10)

Image(r_begin_5+A(10,1):r_end_5+A(10,1),c_end_5+A(10,2):c_end_5+A(10,2)+2047,2) = ROI2;

Image2(r_begin_5+A(10,1):r_end_5+A(10,1),c_end_5+A(10,2):c_end_5+A(10,2)+2047,2) = Alpha;

% Place ROI9 in the upper-left corner - need CC4_9 and CC 6_9 (CC9_4 = 7 and
% CC9_8 = 1)
Row_pers49 = r_begin_5-2048-A(3,1)-A(7,1);
Row_pers89 = r_begin_5-2048-A(9,1)-A(1,1);
Col_pers49 = c_begin_5-2048-A(3,2)-A(7,2);
Col_pers89 = c_begin_5-2048-A(9,2)-A(1,2);

Image(floor((Row_pers49+Row_pers89)/2):floor((Row_pers49+Row_pers89)/2)+2047, floor((Col_pers49+Col_pers89)/2):floor((Col_pers49+Col_pers89)/2)+2047,9) = ROI9;

Image2(floor((Row_pers49+Row_pers89)/2):floor((Row_pers49+Row_pers89)/2)+2047, floor((Col_pers49+Col_pers89)/2):floor((Col_pers49+Col_pers89)/2)+2047,9) = Alpha;

% Place ROI7 in down-left corner - need CC8_7 and CC6_7 (CC8_7 = 2 and
% CC7_6 = 11)
Row_pers87 = r_end_5-A(9,1)+A(2,1);
Row_pers67 = r_end_5+A(4,1)-A(11,1);
Col_pers87 = c_begin_5-2048-A(9,2)+A(2,2);
Col_pers67 = c_begin_5-2048+A(4,2)-A(11,2);

Image(floor((Row_pers87+Row_pers67)/2):floor((Row_pers87+Row_pers67)/2)+2047, floor((Col_pers87+Col_pers67)/2):floor((Col_pers87+Col_pers67)/2)+2047,7) = ROI7;

Image2(floor((Row_pers87+Row_pers67)/2):floor((Row_pers87+Row_pers67)/2)+2047, floor((Col_pers87+Col_pers67)/2):floor((Col_pers87+Col_pers67)/2)+2047,7) = Alpha;


% Place ROI3 in upper-right corner - need CC4_3 and CC2_3 (CC4_3 = 8 and
% CC3_2 = 5)
Row_pers43 = r_begin_5-2048-A(3,1)+A(8,1);
Row_pers23 = r_begin_5-2048+A(10,1)-A(5,1);
Col_pers43 = c_end_5-A(3,2)+A(8,2);
Col_pers23 = c_end_5+A(10,2)-A(5,2);

Image(floor((Row_pers43+Row_pers23)/2):floor((Row_pers43+Row_pers23)/2)+2047, floor((Col_pers43+Col_pers23)/2):floor((Col_pers43+Col_pers23)/2)+2047,3) = ROI3;

Image2(floor((Row_pers43+Row_pers23)/2):floor((Row_pers43+Row_pers23)/2)+2047, floor((Col_pers43+Col_pers23)/2):floor((Col_pers43+Col_pers23)/2)+2047,3) = Alpha;


% Place ROI1 in down-right corner - need CC2_1 and CC6_1 (CC2_1 = 6 and
% CC6_1 = 12)

Row_pers61 = r_end_5+A(4,1)+A(12,1);
Row_pers21 = r_end_5+A(10,1)+A(6,1);
Col_pers61 = r_end_5+A(4,2)+A(12,2);
Col_pers21 = r_end_5+A(10,2)+A(6,2);

Image(floor((Row_pers61+Row_pers21)/2):floor((Row_pers61+Row_pers21)/2)+2047, floor((Col_pers61+Col_pers21)/2):floor((Col_pers61+Col_pers21)/2)+2047,1) = ROI1;

Image2(floor((Row_pers61+Row_pers21)/2):floor((Row_pers61+Row_pers21)/2)+2047, floor((Col_pers61+Col_pers21)/2):floor((Col_pers61+Col_pers21)/2)+2047,1) = Alpha;

%% Create the final image

% Create the final image of cells
Test = sum(Image2,3);
Image2 = Image2./Test;
Image2(isnan(Image2))=0;
Image = bsxfun(@times, Image, Image2);
Image = sum(Image,3);

Image(Image>=205)=255;
Image(Image>=154 & Image<=205)=204;
Image(Image>=103 & Image<=154)=153;
Image(Image>=52 & Image<=103)=102;
Image(Image~=0 & Image<=52)=51;

% Create image with overlap region
Image3 = Image2;
Image3(Image3~=0)=1;
Image3 = sum(Image3,3);
Image3(Image3==1)=0;
Overlap = logical(Image3);

% Create the total area where the image is
Area_image = logical(Test);
