%Script for Biomineralization hyperspectral imaging
%to read in .raw hyperspectral image files and convert them to reflectance
%
%This file is for using the Specim IQ hyperspectral camera.

clear all;
clc;
close all;

pathrefs = cd;
path = pwd;
file = 'ImageData/2018-09-11_004/capture/2018-09-11_004.raw';
header_file = 'ImageData/2018-09-11_004/capture/2018-09-11_004.hdr';
whiteFile = 'ImageData/2018-09-11_004/capture/2018-09-11_004.raw';
darkForWhiteFile = 'ImageData/2018-09-11_004/capture/DARKREF_2018-09-11_004.raw'; %dark reference taken at same exposure time as white
darkForSampleFile = 'ImageData/2018-09-11_004/capture/DARKREF_2018-09-11_004.raw'; %dark reference taken at same exposure time as sample
[wavelengths, spatial, frames, spectral, tint, settings ] = parseHdrInfo_IQ(path,header_file); %Extract wavelength [nm] and other params.

% Read in image file
img = reshapeImage_IQ(path,file);
white = reshapeImage_IQ(path,whiteFile);
dark = reshapeImage_IQ(path,darkForWhiteFile);
dark_sample = reshapeImage_IQ(path,darkForWhiteFile);

%Plot raw images of sample data at a single wavlength
figure
wavlen = 148; %Define wavelength index to plot raw image
imagesc(img(:,:,wavlen)); %View single wavelength of sample
title(['Raw sample image @ \lambda = ',num2str(wavelengths(wavlen)),' nm']);
xlabel('Pixel number');
ylabel('Frame number');
axis 'equal'

%Plot at a single wavelength and select region to consider as white
%reference white reference
figure
imagesc(white(:,:,wavlen)); %View single wavelength of white
title(['Raw white ref image @ \lambda = ',num2str(wavelengths(wavlen)),' nm']);
xlabel('Pixel number');
ylabel('Frame number');
axis 'equal'

[x_w,y_w] = ginput; %Grab region in which to average white reference
x_w = round(x_w);
y_w = round(y_w);



%% Convert to reflectance
%use average frame of white & dark references
% whiteAvg = zeros(spatial,spectral);
% darkavg = whiteAvg;
% img_dark_sub = zeros(frames,spatial,spectral);
% for i = 1:spatial
%     for j = 1:spectral
%         whiteAvg(i,j) = mean(white(:,i,j));
%         darkavg(i,j) = mean(dark(:,i,j));
%         img_dark_sub(:,i,j) = img(:,i,j)-darkavg(i,j);
%     end
% end

whiteAvg = zeros(spatial,spectral); %initialize average row of white pixels
darkavg = whiteAvg; %initialize average row of dark pixels
img_dark_sub = zeros(frames,spatial,spectral); %initialize image matrix - dark pixels
white_sub = zeros(size(white,1),size(white,2),size(white,3));
whiteavg_sub = whiteAvg;
for i = 1:spatial
    for j = 1:spectral
        darkavg(i,j) = mean(dark(:,i,j));
        whiteAvg(i,j) = mean(mean(white(y_w(1):y_w(2),x_w(1):x_w(2),j)));
        whiteavg_sub(i,j) = whiteAvg(i,j)-darkavg(i,j);
        img_dark_sub(:,i,j) = img(:,i,j)-darkavg(i,j);
    end
end
clear dark
clear white
clear img
figure
imagesc(img_dark_sub(:,:,wavlen))
axis 'equal'
title(['Sample-dark image @ \lambda = ',num2str(wavelengths(wavlen)),' nm']);
xlabel('Pixel number');
ylabel('Frame number');

%Reflectance calculation. Normalize value in image to white reference
%R is normalized for each frame to the same row of white reference pixels
R = zeros(size(img_dark_sub));
for k = 1:frames
    for i = 1:spatial
        for j = 1:spectral
            R(k,i,j) = img_dark_sub(k,i,j)./whiteavg_sub(i,j);
        end
    end
end
%plot reflection 2D image at single wavelength
figure
imagesc(R(:,:,wavlen))
axis 'equal'
title(['Norm. Refl @ \lambda = ',num2str(wavelengths(wavlen)),' nm']);
xlabel('Pixel number');
ylabel('Frame number');

%%
%Click two points to define an area to grab spectra for each pixel
[x,y] = ginput;
x = round(x);
y = round(y);
xy = zeros(abs(x(2)-x(1))*abs(y(1)-y(2)),2);
count = 0;
for i = 1:abs(x(2)-x(1))
    for j = 1:abs(y(2)-y(1))
        count = count+1;
        xy(count,:) = [x(1)+i-1,y(1)+j-1];
    end
end

% %Plot individual spectra for each pixel in selected area
figure
for i = 1:length(xy)
    test = R(xy(i,2),xy(i,1),:);
%     peak = max(test);
%     test_norm = test./peak; %normalize all reflectance maximums to 1
    plot(wavelengths, smooth(test(:),1))
    hold on
end

%Save matrix of pixel position points used to extract reflection data
pixel_matx = reshape(xy(:,1),abs(y(1)-y(2)),abs(x(1)-x(2)));
pixel_maty = reshape(xy(:,2),abs(y(1)-y(2)),abs(x(1)-x(2)));

%Save only R spectra for selected pixel positions
R_mat = R(y(1):1:y(2),x(1):1:x(2),:);
lambda_meas = wavelengths;
%% Spectrum to XYZ
% CMF 390 nm to 830 nm 
% Image wavelegnths to index 1:147
load('Z:\Michel\Masters Project\Codes\Calibration\Data\CIE2DegreeObserver.mat')
lambda_q =  wavelengths(1:147) % 390 nm - 830 nm
% Interpolate the CMF
CMF = interp1(CIE2DegreeObserver(:,1),...
                              CIE2DegreeObserver(:,2:end),...
                              lambda_q);
XYZ = zeros(512,512,3);
for i = 1:512
    for j = 1:512
        spectra = reshape(R(i,j,1:147),[1 147]);
        XYZ(i,j,:) = spectra*CMF;R_mat
    end
end
%% XYZ to CIELAB
tic
lab = zeros(size(XYZ));
for i = 1:512
    for j = 1:512
       lab(i,j,:) = xyz2lab(reshape(XYZ(i,j,:), [1 3]),'WhitePoint','a');
    end
end
toc