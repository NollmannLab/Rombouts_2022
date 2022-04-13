function [im_standardized, im] = Normalization(im, FilterSize)

% The normalization step is actually working in two steps :
% 1- the first step is actually an image processing step in order to remove
% as much noise as possible on the image. This step usually improves the
% contrast on the image and allow for a better homogeneous background.
% 2- the second step is the actual normalization. 
% -----------------------------------------------

[Lx,Ly] = size(im);
I = double(im);
I_filtered = I - imgaussfilt(I,FilterSize);
I2 = I .* I;
I2_filtered = sqrt(imgaussfilt(I2,FilterSize));
c = mean(mean(I2_filtered));
I_processed = I_filtered ./ max(c,I2_filtered);
% I_normalized = I_filtered ./ I2_filtered;

I = reshape(I_processed, [Lx*Ly,1]);
std_I = std(I);
I = I/std_I;
Av_I = mean(I);
I = I - Av_I;

im_standardized = reshape(I, [Lx,Ly]);

% After the normalization procedure, the intensity is redistributed in
% order to match the intensity range of a 16-bit image. This step is not
% used for the network but only for saving the processed original image.
% ----------------------------------------------------------------------

I_processed = reshape(I_processed, [Lx*Ly,1]);
[F_I,I] = ecdf(I_processed);

for n_i = 1 : size(I,1)
    if F_I(n_i) > 0.0005
        Min = I(n_i);
        break
    end
end

for n_i = size(I,1) : -1 : 1
    if F_I(n_i) < 0.9995
        Max = I(n_i);
        break
    end
end

a = (2^16-1)/(Max-Min);
b = 1-a*Min;
im = a*I_processed+b;
im (im<=1) = 1;
im (im>=2^16) = 2^16;

im = reshape(im, [Lx,Ly]);
im = uint16(im);