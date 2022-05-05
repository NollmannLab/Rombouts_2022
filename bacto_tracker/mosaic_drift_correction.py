# CODE FOR CROSS CORRELATION CALCULATION (FOR DRIFT CORRECTION) OF ROI5
# ---------------------------------------------------------------------

"""

@author: Sara Rombouts

Created: 01/09/2020

"""


# The input of this code are the Ecoli images (we use the Ecoli images to calculate the drift over time because they change the least over time
# Ecoli cells do not move over time)
# IMPORTANT NOTE: The ROI that is chosen to calculate the drift on should have (enough) elements in the frame for the cross-correlation to be calculated from

# The output of this code ('shifts' - reflects the XY-translation that needs to be applied to achieve maximal CC over images) 
# is the input for the matlab code for drift correction

# This will be a cross-correlation by constantly changing the reference frame

# STEP 1: Import libraries

from pathlib import Path
import os
import matplotlib.pyplot as plt
from skimage.feature import register_translation
import numpy as np

# STEP 2: Load in the Ecoli data
SavingPath = Path('/mnt/PALM_dataserv/DATA/Sara/DATA/TimeLapseData/PredationAssays/2020_05_20-SampleDilutedEcoliColony/003_FastTimeLapse_RAMM_Test')
StartingPath = Path('/mnt/PALM_dataserv/DATA/Sara/DATA/TimeLapseData/PredationAssays/2020_05_20-SampleDilutedEcoliColony/003_FastTimeLapse_RAMM_Test/Segmented_images/ROI_5')
PathEcoli = StartingPath/'Myxo_segmented'

# We retrieve from the Ecoli directory only the tif-files and sort the list (if we do not sort, then the list is not numbered properly 
# which will result in an incorrect looping through the files)
EcoliFiles = os.listdir(PathEcoli)
list.sort(EcoliFiles)
EcoliFiles = [i for i in EcoliFiles if i.endswith('.tif')]

# STEP 3: Load in Reference image - the first image in the Ecoli sequence will be the reference image for cross correlation

ImgRef = plt.imread(PathEcoli/EcoliFiles[0])

plt.imshow(ImgRef)

# STEP 4: Calculate cross-correlation for each Ecoli image that follows - 
# we loop through sorted list of images and construct the matrix that will indicate the XY offset of the images over time with respect to the 
# reference image (pixel precision)

for file in range(len(EcoliFiles)):

    Img = plt.imread(PathEcoli/EcoliFiles[file])

    shift, error, diffphase = register_translation(ImgRef, Img)
    
    if file==0:
        shifts = shift;
    else:
        shifts = np.vstack((shifts,shift))
        
print(shifts)

# The shifts output will be used as input for the matlab code (DriftCorrection_TestPythonDriftCorr.m) to correct all the images
# Create a filename which will change with each loop
filename = str(SavingPath)+'/XYshift_DriftCorr.txt'
    
# Save the text file in the newly created '/Tiling/' directory
np.savetxt(filename, shifts)
    