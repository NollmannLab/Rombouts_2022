function [Conn, M] = PostProcess_TiledIm(New_Image)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code: Sara Rombouts (CBS, Team marcelo Nollmann)
%
% 31/08/2020
%
% Goal of code: function to postprocess the segmented images (RAMM
% microscope - segmetned with DeepLearning Network)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the size of the image
Lx = size(New_Image,1);
Ly = size(New_Image,2);

%% STRINGENCY PIXEL SELECTION

Im_Myxo = New_Image;
Im_Myxo_255 = Im_Myxo;
Im_Myxo_255(Im_Myxo_255~=255) = 0;
Im_Myxo_255 = logical(Im_Myxo_255);

Im_Myxo_255 = bwmorph(Im_Myxo_255, 'hbreak', 1);
Im_Myxo_255 = bwmorph(Im_Myxo_255, 'open', 1);
Im_Myxo_255 = bwpropfilt(Im_Myxo_255, 'area', [1 1000000],4);
Labelmatrix_255 = labelmatrix(bwconncomp(Im_Myxo_255,4));
Props_255 = regionprops('table', Labelmatrix_255, 'Area', 'PixelIdxList');

Im_Myxo_204 = Im_Myxo;
Im_Myxo_204(Im_Myxo_204<=203) = 0;
Im_Myxo_204 = logical(Im_Myxo_204);
Labelmatrix_204 = labelmatrix(bwconncomp(Im_Myxo_204,4));
Props_204 = regionprops('table', Labelmatrix_204, 'Area', 'PixelIdxList');

Im_Myxo_153 = Im_Myxo;
Im_Myxo_153(Im_Myxo_153<=152) = 0;
Im_Myxo_153 = logical(Im_Myxo_153);
Labelmatrix_153 = labelmatrix(bwconncomp(Im_Myxo_153,4));
Props_153 = regionprops('table', Labelmatrix_153, 'Area', 'PixelIdxList');

Im_Myxo_102 = Im_Myxo;
Im_Myxo_102(Im_Myxo_102<=101) = 0;
Im_Myxo_102 = logical(Im_Myxo_102);
Labelmatrix_102 = labelmatrix(bwconncomp(Im_Myxo_102,4));
Props_102 = regionprops('table', Labelmatrix_102, 'Area', 'PixelIdxList');

Im_Myxo_51 = Im_Myxo;
Im_Myxo_51(Im_Myxo_51<=50) = 0;
Im_Myxo_51 = logical(Im_Myxo_51);
Labelmatrix_51 = labelmatrix(bwconncomp(Im_Myxo_51,4));
Props_51 = regionprops('table', Labelmatrix_51, 'Area', 'PixelIdxList');

%% filters by area

Area = [];
for i = 1:size(Props_255(:,1),1)
    Area_255 = table2array(Props_255(i,1));
    
    % T_... gets you the CellID (or line in the Props_-matrix where the area
    % needs to be retrieved from)

    % Area_... gets you the area size of the segment for a given stringency
    % level

    T_204 = unique(nonzeros(Labelmatrix_204(cell2mat(Props_255{i,2}))));
    Area_204 = table2array(Props_204(T_204,1));
    
    T_153 = unique(nonzeros(Labelmatrix_153(cell2mat(Props_255{i,2}))));
    Area_153 = table2array(Props_153(T_153,1));
    
    T_102 = unique(nonzeros(Labelmatrix_102(cell2mat(Props_255{i,2}))));
    Area_102 = table2array(Props_102(T_102,1));
    
    T_51 = unique(nonzeros(Labelmatrix_51(cell2mat(Props_255{i,2}))));
    Area_51 = table2array(Props_51(T_51,1));
    
    Area = [Area_255;Area_204;Area_153;Area_102; Area_51];

    % Save all CellIDs in each stringency level
    Total_2 = [i; T_204; T_153; T_102; T_51];
    
    % Define whether segment is less or more than 5% varying in size when
    % lowering the stringency level
    R=0;
    Total = [];
    for j = 1:size(Area)-1
        R = abs(Area(j+1)-Area(j))/Area(j);
        Total = [Total;R];
        if R > 0.05
            break
        end
    end
    Test{i} = Total;
    Test2{i} = Total_2;
end


New = zeros(Lx*Ly,1);
s= [];
for i = 1:size(Test,2)
    if size(Test{i},1) == 1
        stringency = 255;
        s = [s; stringency];
        Seg = Test2{i}(1);
        PixelList = cell2mat(Props_255{Seg,2});
        New(PixelList) = 1;
    elseif size(Test{i},1) == 2
        stringency = 204;
        s = [s; stringency];
        Seg = Test2{i}(2);
        PixelList = cell2mat(Props_204{Seg,2});
        New(PixelList) = 1;
    elseif size(Test{i},1) == 3
        stringency = 153;
        s = [s; stringency];
        Seg = Test2{i}(3);
        PixelList = cell2mat(Props_153{Seg,2});
        New(PixelList) = 1;
    elseif size(Test{i},1) == 4 && Test{i}(end)>=0.051
        stringency = 102;
        s = [s; stringency];
        Seg = Test2{i}(4);
        PixelList = cell2mat(Props_102{Seg,2});
        New(PixelList) = 1;
    else size(Test{i},1) == 4 & Test{i}(end)<=0.051;
        stringency = 51;
        s = [s; stringency];
        Seg = Test2{i}(5);
        PixelList = cell2mat(Props_51{Seg,2});
        New(PixelList) = 1;
    end
end
M = reshape(New, [Lx Ly]);

%%
% 4-CONNECTIVITY is used BECAUSE SOME SEGMENTS TOUCH EACH when using 8-CONNECTIVITY

Compl = imcomplement(logical(M));
DistFromZero = bwdist(Compl);
DistFromZero(DistFromZero==1) = 0;
Final = bwmorph(DistFromZero, 'Thicken', 1);
Final = bwpropfilt(Final, 'area', [100 1000000],4);
Final_conn = bwconncomp(Final,4);
Final_labelmatrix = labelmatrix(Final_conn);
Im_thin = bwmorph(Final, 'thin', Inf);
Im_spur = bwmorph(Im_thin, 'spur', 3);
Im_branch = bwmorph(Im_spur, 'branchpoints');

New = reshape(Final_labelmatrix, [Lx*Ly 1]);
for n = 1 : Final_conn.NumObjects
    PixelList = Final_conn.PixelIdxList{1,n};
    Num_BP = Im_branch(PixelList);
    Num_BP = nonzeros(Num_BP);
    if size(Num_BP,1)~=0
        New(PixelList)=0;
    end
end
Good_Seg = reshape(New, [Lx Ly]);
Conn_1 = bwconncomp(Good_Seg,4);
Conn_Label = labelmatrix(Conn_1);

%% Turtuosity of the segment to filter out all segments that are either fused or too curved

Images = regionprops('table', Conn_1, 'Image');
Dist = [];
Tortuosity = [];
for k = 1:Conn_1.NumObjects
    Image = cell2mat(Images{k,1});
    Image_thin = bwmorph(Image, 'thin', Inf);
    Image_thin = padarray(Image_thin, [1 1], 0, 'both');
    Image_spur = bwmorph(Image_thin, 'spur', 3);
    Image_Ends = bwmorph(Image_spur, 'endpoints');
    Length = size(find(Image_thin),1);
    [row,col] = find(Image_Ends);

    if ~isempty(row) && size(row,1)~=1
        % define the eucledian distance between those end points

        EucDist = sqrt(((diff(row)^2+diff(col)^2)));
        Dist = [Dist, EucDist];

        % calculate the turtoistity
        Tort = Length/EucDist;
        elseif size(row,1)==1 % Discards all segments for which just 1 endpoint can be found (= case when cell is very short and 'spur' operation reduces the segment to 1 dot)
        Tort = 10;
    else
        Tort = 10;
    end % Discards all segments for which no endpoints can be found (= case where the segment has a round shape)

    Tortuosity = [Tortuosity; Tort];

end

% Define threshold for what we allow
Threshold = 2;
Del_Seg = find(Tortuosity>=Threshold);
for p = 1:size(Del_Seg,1)
    P = Del_Seg(p);
    Conn_Label(Conn_Label==P)=0;
end
Conn = bwconncomp(Conn_Label,4);