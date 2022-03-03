# Comerci_Gillman_CellSystems_2022
MATLab Code From Comerci, Gillman, et al. Cell Systems. 2022.

Title: Localized electrical stimulation triggers cell-type-specific proliferation in biofilms
Authors: Comerci CJ, Gillman AL, Galera-Laporta L, Gutierrez E, Groisman A, Larkin JW, Garcia-Ojalvo J, Suel GM.

This document provides a brief overview of the published code for the above paper. Note that more detailed documentation is provided in comments inside the individual files.

We are posting the two novel pieces of code developed for this manuscript. They both classify pixels in two color images into two anti-correlated binary maps (i.e. each pixel is determined to be either green or magenta).

1. PixelClassification_AntiCorr_Confocal.m - This code uses histogram normalization to easily classify pixels from two anti-correlated images with high contrast. We used this on confocal data, where the narrow focal plane provided high anti-correlation between our two color channels (i.e. pixels are mostly magenta or green with little signal in the other color channel).

2. PixelClassification_AntiCorr_WF.m - This code thresholds the log-ratio signal levels to classify pixels from two anti-correlated images with lower contrast. This was used on wide field data, where the contrast between the two color channels is much lower (i.e. all pixels have some signal in both channels, and are just slightly brighter in one or the other).

All code can be run and prompts the user for inputs. The inputs are detailed in the comments at the beginning of the code. There are a couple plots at the end of the code that allows the user to monitor how the pixel classification worked for individual timepoints and what the result was as a function of time.
