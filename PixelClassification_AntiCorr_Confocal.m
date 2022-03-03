%% PixelClassification_AntiCorr_Confocal.m
% This code loads in a time series of 2-color, anticorrelated images and 
% classifies pixels into 2 logical maps. The code:
% 1. Asks the user to select a folder and the first image for each color
% channel
% 2. Asks the user to identify an ROI of background (which is subtracted
% from each image) and an ROI for analysis
% 3. Loads and preprocesses the images
% 4. Uses histogram normalization to flatten and normalize the signal
% 5. Compares the signal in the two color channels, and classifies each
% pixel above background as either in one color channel or the other

% This code takes a little time to run. It is designed to run on a single
% position at a time, saving the output for further analysis. At the bottom
% are several plots/images to watch and verify the output that can be
% optionally run

% Inputs:
% User will be asked to select:
% 1. The number of timepoints in the timelapse
% 2. The folder to save results
% 3. The first YFP image
% 4. The first RFP image
% 5. The upperleft and lower right corners of a box containing background
% (i.e. outside the biofilm)
% 6. The upperleft and lower right corners of a box that contains the entire
% region to be analyzed

% Outputs:
% Saves a .mat file containing:
% 1. Map - a logical map of signal above background (ie cells/biofilm),
% time is encoded in the 3D
% 2. YFP_Map/RFP_Map - a logical map of of pixels in each color channel,
% time is encoded in the 3D
% 3. YFP_Frac/RFP_Frac - the fraction of pixels in each color channel,
% time is encoded in the 2D
% ROIpts - the opposing corners used to select the ROI for analysis
% Y(R)FP_file/_path - the file and path for the YFP and RFP images
% YFP_BG/RFP_BG - the YFP and RFP background values

% Created 5/15/2021 By Colin J. Comerci - Suel Lab - UCSD
% Modified 2/16/2022 to simplify inputs/outputs and add comments
%%
close all
clear all
% Number of timepoints in image series
num_t = input('Number of timepoints?');
% Select directory to save output
Path_Output = uigetdir(matlabroot,'Select Directory to Save Output');
% Select 1st frane of YFP and RFP image files
[YFP_file,YFP_path] = uigetfile('*.tif','Select 1st YFP Image File');
[RFP_file,RFP_path] = uigetfile('*.tif','Select 1st RFP Image File');


% Timepoint 50 is used to select an ROI for 
% 1. background subtraction - select a region with no biofilm/cells
% 2. the full region for analysis - removing the pillar
% Click the upper left and bottom right corners of a box

% Note: this assumes the file ends with a 3 digit time code and that T=50
% is a good frame to find a region with no biofilm present
cd(YFP_path)
Im_temp_RFP = double(imread([YFP_file(1:end-7) '050.tif']));
cd(RFP_path) 
Im_temp_YFP = double(imread([RFP_file(1:end-7) '050.tif']));
imagesc(Im_temp_RFP+Im_temp_YFP)
ROIpts_BG = round(ginput(2));

% This assumes anything less than the mean + 4STD of a region outside the
% biofilm is background. This is used to identify parts of the image
% outside the biofilm. You can use the image overlay below to check if this
% works for your images
RFP_BG = mean2(Im_temp_RFP(ROIpts_BG(1,2):ROIpts_BG(2,2),ROIpts_BG(1,1):ROIpts_BG(2,1)))...
    + (4*std2(Im_temp_RFP(ROIpts_BG(1,2):ROIpts_BG(2,2),ROIpts_BG(1,1):ROIpts_BG(2,1))));
YFP_BG = mean2(Im_temp_YFP(ROIpts_BG(1,2):ROIpts_BG(2,2),ROIpts_BG(1,1):ROIpts_BG(2,1)))...
    + (4*std2(Im_temp_YFP(ROIpts_BG(1,2):ROIpts_BG(2,2),ROIpts_BG(1,1):ROIpts_BG(2,1))));

ROIpts = round(ginput(2));
% Cycle through the timepoints
for i = 1:num_t
    clear Im_YFP_temp Im_RFP_temp YFP_cont RFP_cont
    cd(YFP_path)
        Im_YFP_temp = imfilter(double(imread([YFP_file(1:end-7) sprintf('%03d',i) '.tif'],'PixelRegion',{[ROIpts(1,2) ROIpts(2,2)],[ROIpts(1,1) ROIpts(2,1)]})),Filt,'symmetric');
      cd(RFP_path)  
        Im_RFP_temp = imfilter(double(imread([RFP_file(1:end-7) sprintf('%03d',i) '.tif'],'PixelRegion',{[ROIpts(1,2) ROIpts(2,2)],[ROIpts(1,1) ROIpts(2,1)]})),Filt,'symmetric');
    % Subtract background from each image and set negative values to zero
    Im_YFP_temp = Im_YFP_temp-YFP_BG;
    Im_YFP_temp(Im_YFP_temp<0) = 0;
    Im_RFP_temp = Im_RFP_temp-RFP_BG;
    Im_RFP_temp(Im_RFP_temp<0) = 0;
    [m,n] = size(Im_YFP_temp);
    % adapthisteq is used to flatten the images (removing the dimming of
    % signal at the biofilm edges) and to normalize their signal levels to
    % ease comparison.
    % Note: 
    % -This normalizes the image by the maximum value in each channel.
    % Other clipping based normalization methods (i.e. normalizing to some
    % quantile) caused problems with the pixel classification. 
    % -Using the uniform distribution behaved similarly to the rayleigh
    % distribution in preliminary tests
    YFP_cont = adapthisteq(Im_YFP_temp./max(max(Im_YFP_temp)),'NumTiles',[round(m/128) round(n/128)],'Distribution','rayleigh');
    RFP_cont = adapthisteq(Im_RFP_temp./max(max(Im_RFP_temp)),'NumTiles',[round(m/128) round(n/128)],'Distribution','rayleigh');
    Im_YFP(:,:,i) = YFP_cont;
    Im_RFP(:,:,i) = RFP_cont;
    % Create a "cell/biofilm map" where any signal above background is
    % considered to be part of cells/biofilm
    Map(:,:,i) = (Im_RFP_temp+Im_YFP_temp)>0;
    Map(:,:,i) = bwareaopen(Map(:,:,i),4);
    % Since the image histograms are normalized, their signal can be
    % compared, with higher signals in YFP or RFP causing the pixel to be
    % added to the appropriate map
    YFP_Map(:,:,i) = and(Im_YFP(:,:,i)./Im_RFP(:,:,i)>=1,Map(:,:,i));
    RFP_Map(:,:,i) = and(Im_YFP(:,:,i)./Im_RFP(:,:,i)<1,Map(:,:,i));
    YFP_Frac(i) = sum(sum(YFP_Map(:,:,i)))./sum(sum(Map(:,:,i)));
    RFP_Frac(i) = sum(sum(RFP_Map(:,:,i)))./sum(sum(Map(:,:,i)));
end
% Save to a file. You may want to change the folder location
cd(Path_Output)
% This saves the data to an output file with a base name derived from the
% YFP image. This worked well for data with a name that has a trailing
% "_C001T001.tif" to denote the color channel and timepoint
save([YFP_file(1:end-13) '.mat'],'Map','YFP_Map','RFP_Map','YFP_Frac',...
    'RFP_Frac','ROIpts','YFP_file','YFP_path','RFP_file','RFP_path','YFP_BG',...
    'RFP_BG')
%% Plots the YFP and RFP fractions as a function of time
figure
plot(YFP_Frac,'g')
hold on
plot(RFP_Frac,'m')
xlabel('Time (in frames)')
ylabel('Fraction of Image')
%% Shows the images and maps at a single timepoint to check result
% this is the timepoint to be plotted. Change as needed
i = 75;
figure
subplot(1,2,1)
imshowpair(Im_YFP(:,:,i),Im_RFP(:,:,i),'ColorChannels','green-magenta')
title(['Raw Images at T = ' num2str(i)])
subplot(1,2,2)
title('Maps')
imshowpair(YFP_Map(:,:,i),RFP_Map(:,:,i),'ColorChannels','green-magenta')