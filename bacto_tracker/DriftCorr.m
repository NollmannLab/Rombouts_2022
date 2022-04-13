function [New_Image, New_Overlap, New_Area_image] = DriftCorr(Image, Overlap, Area_image, Drift,t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code: Sara Rombouts (CBS, Team marcelo Nollmann)
%
% 31/08/2020
%
% Goal of code: function to drift correct the tiled images
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find max drift in XYDrift_CrossCorr (python) - this is the amount
% of rows-cols of zeros that will be padded around the image (to achieve circshift)
        
Lim = max(max(abs(Drift)));
Lim = Lim-40;
if Lim <=0
    Lim = 0;
end

% Pad max of absolute values of total drift correlations MINUS 40 pixels (= minimal overlap that will be available already padded around the image) to the array on each side
New_Image = padarray(Image,[Lim Lim],0,'both');
New_Overlap = padarray(Overlap,[Lim Lim],0,'both');
New_Area_image = padarray(Area_image,[Lim Lim],0,'both');

% Circlshift based on the XYDrift_CrossCorr of timestamp calculated for ROI5
New_Image = circshift(New_Image, Drift(t,:));
New_Overlap = circshift(New_Overlap, Drift(t,:));
New_Area_image = circshift(New_Area_image, Drift(t,:));



