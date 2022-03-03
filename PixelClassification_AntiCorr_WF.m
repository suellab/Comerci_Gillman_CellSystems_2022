%% PixelClassification_AntiCorr_WF.m
% This code loads in a time series of 2-color, anticorrelated widefield 
% images and classifies pixels into 2 logical maps. The code:
% 1. Asks the user to select a tiff stack of stabalized images for the 
% phase image (used to find the electrode center) 
% 2. Asks the user to select a tiff stack of stabalized images for each of 
% the two color channels
% 3. Asks the user to input a starting and ending timepoint for analysis
% 4. Asks the user to input a pixel size (in um)
% 5. Asks the user to identify an ROI of background (which is subtracted
% from each image) and an ROI for analysis
% 6. Loads and preprocesses the images
% 7. Asks the user to select a log-ratio quantile to separate the two
% channels from a gallery of quantiles
% 8. Uses that quantile on a region away from the electrode to determine
% the log-ratio threshold to segment the region around the electrode

% Note that the widefield images used in this study had far less
% distinction between the two anti-correlated color channels than the
% confocal images, as is to be expected since the widefield images include
% information from cells above and below the imaging plane. The
% fluorescence intensities also varied quite a bit over time, making
% thresholding uniquely challenging. The strategy employed here assumes the
% area fractions of each color channel in a region far from the electrode
% remain relatively constant over time, and that this region can then be
% used to identify the appropriate threshold for the region near the
% electrode. This worked well for our images over ~12 hours.

% Inputs:
% User will be asked to select:
% 1. The phase tiff stack
% 2. The YFP tiff stack
% 3. The RFP tiff stack
% 4. The initial timepoint
% 5. The final timepoint
% 6. The pixel size (in um)
% 7. The upperleft and lower right corners of a box containing background
% (i.e. outside the biofilm)
% 8. The center of the electrode to be analyzed
% 9. The electrode "name" (i.e. its coordinate pair)
% 10. The appropriate quantile (decimal from 0-1)

% Outputs:
% Saves a .mat file containing:
% 1. YFP_BG/RFP_BG - the YFP and RFP background values
% 2. ROIpts - the center of the electrode being analyzed (used to determine
% the region of the image that is loaded)
% 3. ElecProp - the more precise center of the electrode from the subregion 
% of the input image
% 4. thresh - the thresholds determined from the user selected quantiles
% at each timepoint
% 5. YFP_Map/RFP_Map - a logical map of of pixels in each color channel,
% time is encoded in the 3D
% 6. YFP_Frac/RFP_Frac - the fraction of pixels in each color channel,
% time is encoded in the 2D
% 7. Map - a logical map of signal above background (ie cells/biofilm),
% time is encoded in the 3D
% 8. Y(R)FP_file/_path - the file and path for the YFP and RFP images

% Created 5/15/2021 By Colin J. Comerci - Suel Lab - UCSD
% Modified 2/27/2022 to simplify inputs/outputs and add comments
%% Select Files
close all
clear all
% The following 4 lines allow the user to select images (Phase, YFP, and
% RFP), and input timepoints to be analyzed. It assumes the images are tif
% stacks that have been stabalized
[Phase_file,Phase_path] = uigetfile('*.tif','Select Phase tif stack');
[YFP_file,YFP_path] = uigetfile('*.tif','Select YFP tif stack');
[RFP_file,RFP_path] = uigetfile('*.tif','Select RFP tif stack');
t_i = input('Starting timepoint?');
t_f = input('Ending timepoint?');
pix = input('Pixel size (in um)?');

%%
% The Gaussian filter used to smooth/denoise the raw images
Filt = fspecial('gaussian',round(5*1),1);
% The starting timepoint is used to identify background (i.e. a region 
% outside the biofilm)
cd(RFP_path)
Im_temp_RFP = double(imread(RFP_file,t_i));
cd(YFP_path) 
Im_temp_YFP = double(imread(YFP_file,t_i));
imagesc(Im_temp_RFP+Im_temp_YFP)
% Select top left and bottom right of a region of background
ROIpts_BG = round(ginput(2));
% Background is just the mean
RFP_BG = mean2(Im_temp_RFP(ROIpts_BG(1,2):ROIpts_BG(2,2),ROIpts_BG(1,1):ROIpts_BG(2,1)));
YFP_BG = mean2(Im_temp_YFP(ROIpts_BG(1,2):ROIpts_BG(2,2),ROIpts_BG(1,1):ROIpts_BG(2,1)));
% Select the center of the electrode to be analyzed
ROIpts = round(ginput(1));
ElectrodeName = input('Electrode Number (E_ _)?','s');
cd(Phase_path)
Im_temp_Ph = double(imread(Phase_file,t_i,'PixelRegion',{[ROIpts(1,2)-70 ROIpts(1,2)+70],[ROIpts(1,1)-70 ROIpts(1,1)+70]}));
% The phase image is used to identify the electrode, and then find it's
% centroid. The threshold may need to be changed here to accurately find
% the electrode center, but is not too sensitive
BW = Im_temp_Ph<2500;
ElecProp = regionprops(bwareaopen(BW,50),'Centroid');
%% Load and preprocess images
for i = t_i:t_f
    clear Im_YFP_temp Im_RFP_temp YFP_cont RFP_cont Map_out_temp
    cd(YFP_path) 
        Im_YFP_temp = imfilter(double(imread(YFP_file,i,'PixelRegion',{[ROIpts(1,2)-70 ROIpts(1,2)+70],[ROIpts(1,1)-70 ROIpts(1,1)+70]})),Filt,'symmetric');
        cd(RFP_path)
        Im_RFP_temp = imfilter(double(imread(RFP_file,i,'PixelRegion',{[ROIpts(1,2)-70 ROIpts(1,2)+70],[ROIpts(1,1)-70 ROIpts(1,1)+70]})),Filt,'symmetric');

    % Subtract background from each image and set negative values to zero
    Im_YFP_temp = Im_YFP_temp-YFP_BG;
    Im_YFP_temp(Im_YFP_temp<0) = 0;
    Im_RFP_temp = Im_RFP_temp-RFP_BG;
    Im_RFP_temp(Im_RFP_temp<0) = 0;
    Im_YFP(:,:,i-t_i+1) = Im_YFP_temp;
    Im_RFP(:,:,i-t_i+1) = Im_RFP_temp;
end
%% Find threshold 
% Find maps for regions near and away from the electrode
[X,Y] = meshgrid(1:141,1:141);
dist = pix.*sqrt((X-ElecProp.Centroid(1)).^2+(Y-ElecProp.Centroid(2)).^2);
% Map of the region near the electrode
Map_in = dist>15 & dist<55;
% Map of the region away from the electrode
Map_out = dist>55 & dist<75;
Map_out_temp = (Im_RFP(:,:,1)+Im_YFP(:,:,1))>200;
Map_out_temp = and(bwareaopen(Map_out_temp,4),Map_out);
Map_elec = dist>15;
Im_YFP_temp = Im_YFP(:,:,1);
Im_RFP_temp = Im_RFP(:,:,1);
% Calculates the ratio of the logs of the YFP and RFP intensities, which will be
% used to classify pixels later
ratio_temp = log(Im_YFP_temp(Map_out_temp))./log(Im_RFP_temp(Map_out_temp));
figure
% Show the maps for different quantiles of the log intensity ratio
for i = 1:5
    Quant_temp(i) = quantile(ratio_temp,0.4+i*0.1);
    YFP_Map_temp(:,:,i) = log(Im_YFP_temp)./log(Im_RFP_temp)>=Quant_temp(i);
    RFP_Map_temp(:,:,i) = log(Im_YFP_temp)./log(Im_RFP_temp)<Quant_temp(i);
    subplot(2,5,i)
    imshowpair(Im_YFP_temp./max(max(Im_YFP_temp)),Im_RFP_temp./max(max(Im_RFP_temp)),'ColorChannels','green-magenta')
    title(num2str(0.4+i*0.1))
    subplot(2,5,i+5)
    imshowpair(YFP_Map_temp(:,:,i),RFP_Map_temp(:,:,i),'ColorChannels','green-magenta')
end
% Allows user to select the quantile that best splits the image into two
% pixel classifications (i.e. how many pixels at the initial timepoint are
% in each category)
quant = input('Select quantile=');

%% Create pixel maps
for i = 1:t_f-t_i+1
    clear ratio_temp Im_YFP_temp Im_RFP_temp Map_out_temp Map_elec_temp
    % Create map of pixels with signal
    Map(:,:,i) = (Im_RFP(:,:,i)+Im_YFP(:,:,i))>200;
    Map(:,:,i) = and(bwareaopen(Map(:,:,i),4),Map_in);
    Map_out_temp = (Im_RFP(:,:,i)+Im_YFP(:,:,i))>200;
    Map_elec_temp = and(bwareaopen(Map_out_temp,4),Map_elec);
    Map_out_temp = and(bwareaopen(Map_out_temp,4),Map_out);
    % Create the log-intensity ratio for the pixels away from the electrode
    Im_YFP_temp = Im_YFP(:,:,i);
    Im_RFP_temp = Im_RFP(:,:,i);
    ratio_temp = log(Im_YFP_temp(Map_out_temp))./log(Im_RFP_temp(Map_out_temp));
    % Based on the quantile selected by the user above, the necessary
    % threshold to split the pixels away from the electrode into the same
    % fraction is calculated
    thresh(i) = quantile(ratio_temp,quant);
    % The calculated threshold is used to create pixel classification maps,
    % but now of the pixels near the electrode
    YFP_Map(:,:,i) = and(log(Im_YFP(:,:,i))./log(Im_RFP(:,:,i))>=thresh(i),Map(:,:,i));
    RFP_Map(:,:,i) = and(log(Im_YFP(:,:,i))./log(Im_RFP(:,:,i))<thresh(i),Map(:,:,i));
    YFP_Frac(i) = sum(sum(YFP_Map(:,:,i)))./sum(sum(Map(:,:,i)));
    RFP_Frac(i) = sum(sum(RFP_Map(:,:,i)))./sum(sum(Map(:,:,i)));
    % This is the full pixel classification map
    YFP_Map_full(:,:,i) = and(log(Im_YFP(:,:,i))./log(Im_RFP(:,:,i))>=thresh(i),Map_elec_temp);
    RFP_Map_full(:,:,i) = and(log(Im_YFP(:,:,i))./log(Im_RFP(:,:,i))<thresh(i),Map_elec_temp);
end

%% Save Data
% Save to a file. You may want to change the folder location
cd(Phase_path)
save([Phase_file(1:end-4) '_' ElectrodeName '.mat'],'RFP_BG','YFP_BG','ROIpts',...
    'ElecProp','thresh','YFP_Map','RFP_Map','YFP_Frac','RFP_Frac','Map',...
    'YFP_file','YFP_path','RFP_file','RFP_path','Phase_file','Phase_path')

%% Plots the YFP and RFP fractions as a function of time
figure
plot(YFP_Frac,'g')
hold on
plot(RFP_Frac,'m')
xlabel('Time (in frames)')
ylabel('Fraction of Pixels')
%% Plot Pixel Classification Maps
clear Im_YFP_temp Im_RFP_temp
% Set timepoint
i =36;

Im_YFP_temp = Im_YFP(:,:,i);
Im_YFP_temp(~Map_in) = 0;
Im_RFP_temp = Im_RFP(:,:,i);
Im_RFP_temp(~Map_in) = 0;

figure
subplot(1,2,1)
imshowpair(Im_YFP_temp./max(max(Im_YFP_temp)),Im_RFP_temp./max(max(Im_RFP_temp)),'ColorChannels','green-magenta')
title(['Raw Images at T = ' num2str(i)])
subplot(1,2,2)
imshowpair(YFP_Map(:,:,i),RFP_Map(:,:,i),'ColorChannels','green-magenta')
title('Maps')
%% Plots the histograms of the log-intensity ratio
figure
clear Im_YFP_temp Im_RFP_temp
% set timepoint
i = 19;
Im_YFP_temp = Im_YFP(:,:,i);
Im_RFP_temp = Im_RFP(:,:,i);
histogram(log(Im_YFP_temp(Map_in))./log(Im_RFP_temp(Map_in)),50)