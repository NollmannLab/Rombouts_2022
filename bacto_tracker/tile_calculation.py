# CODE FOR CALCULATION OF CROSS-CORRELATION FOR TILING BASED ON BF_NORMALIZED AND MYXO_SEGMETEND IMAGES
# -----------------------------------------------------------------------------------------------------

"""
@author: Sara Rombouts

Created: 01/09/2020

Adjustment of code: MosaicImages.py - Here we save the cross-correlations calculated based on BF_normalized images,
however, when this CC calculation fails (for 2 neighboring images), we calculate only for those images the CC
based on the Myxo_segmented images
We save 24 CCs 
- the first 12 are the ones based on the BF_normalized images (result = 0,0 is calculation fails)
- the last 12 are the ones for which the previous calc fails base don Myxo_segmented (result = 0,0 is calc based on BF_normalized gave result and is not 0,0 when calc is based on Myxo_segmented)
"""

# from pathlib import Path
import os
import matplotlib.pyplot as plt
from skimage.feature import register_translation
import numpy as np
import glob


path_to_data = '/mnt/grey/DATA/rawData_2020/Experiment_9_sara/Segmentation/'

#Create a directory in which the txt-files should be stored
new_dir_name = path_to_data+'Tiling'
os.mkdir(new_dir_name)

Files_ROI1 = glob.glob(path_to_data+'Segmented_images/ROI_1/BF_normalized/*.tif')
Files_ROI2 = glob.glob(path_to_data+'Segmented_images/ROI_2/BF_normalized/*.tif')
Files_ROI3 = glob.glob(path_to_data+'Segmented_images/ROI_3/BF_normalized/*.tif')
Files_ROI4 = glob.glob(path_to_data+'Segmented_images/ROI_4/BF_normalized/*.tif')
Files_ROI5 = glob.glob(path_to_data+'Segmented_images/ROI_5/BF_normalized/*.tif')
Files_ROI6 = glob.glob(path_to_data+'Segmented_images/ROI_6/BF_normalized/*.tif')
Files_ROI7 = glob.glob(path_to_data+'Segmented_images/ROI_7/BF_normalized/*.tif')
Files_ROI8 = glob.glob(path_to_data+'Segmented_images/ROI_8/BF_normalized/*.tif')
Files_ROI9 = glob.glob(path_to_data+'Segmented_images/ROI_9/BF_normalized/*.tif')

FilesMyxo_ROI1 = glob.glob(path_to_data+'Segmented_images/ROI_1/Myxo_segmented/*.tif')
FilesMyxo_ROI2 = glob.glob(path_to_data+'Segmented_images/ROI_2/Myxo_segmented/*.tif')
FilesMyxo_ROI3 = glob.glob(path_to_data+'Segmented_images/ROI_3/Myxo_segmented/*.tif')
FilesMyxo_ROI4 = glob.glob(path_to_data+'Segmented_images/ROI_4/Myxo_segmented/*.tif')
FilesMyxo_ROI5 = glob.glob(path_to_data+'Segmented_images/ROI_5/Myxo_segmented/*.tif')
FilesMyxo_ROI6 = glob.glob(path_to_data+'Segmented_images/ROI_6/Myxo_segmented/*.tif')
FilesMyxo_ROI7 = glob.glob(path_to_data+'Segmented_images/ROI_7/Myxo_segmented/*.tif')
FilesMyxo_ROI8 = glob.glob(path_to_data+'Segmented_images/ROI_8/Myxo_segmented/*.tif')
FilesMyxo_ROI9 = glob.glob(path_to_data+'Segmented_images/ROI_9/Myxo_segmented/*.tif')

list.sort(Files_ROI1)
list.sort(Files_ROI2)
list.sort(Files_ROI3)
list.sort(Files_ROI4)
list.sort(Files_ROI5)
list.sort(Files_ROI6)
list.sort(Files_ROI7)
list.sort(Files_ROI8)
list.sort(Files_ROI9)

list.sort(FilesMyxo_ROI1)
list.sort(FilesMyxo_ROI2)
list.sort(FilesMyxo_ROI3)
list.sort(FilesMyxo_ROI4)
list.sort(FilesMyxo_ROI5)
list.sort(FilesMyxo_ROI6)
list.sort(FilesMyxo_ROI7)
list.sort(FilesMyxo_ROI8)
list.sort(FilesMyxo_ROI9)

for i in range(len(Files_ROI1)):
  
    A1 = plt.imread(Files_ROI1[i])
    A2 = plt.imread(Files_ROI2[i])
    A3 = plt.imread(Files_ROI3[i])
    A4 = plt.imread(Files_ROI4[i])
    A5 = plt.imread(Files_ROI5[i])
    A6 = plt.imread(Files_ROI6[i])
    A7 = plt.imread(Files_ROI7[i])
    A8 = plt.imread(Files_ROI8[i])
    A9 = plt.imread(Files_ROI9[i])
    
    # from top to bottom
        
    shift9_8, error, diffphase = register_translation(A9, A8)
    if np.all(shift9_8==0):
        B9 = plt.imread(FilesMyxo_ROI9[i])
        B8 = plt.imread(FilesMyxo_ROI8[i])
        shiftB9_8, error, diffphase = register_translation(B9, B8)
    else:
        shiftB9_8 = 0,0
    
    
    shift8_7, error, diffphase = register_translation(A8, A7)
    if np.all(shift8_7==0):
        B8 = plt.imread(FilesMyxo_ROI8[i])
        B7 = plt.imread(FilesMyxo_ROI7[i])
        shiftB8_7, error, diffphase = register_translation(B8, B7)
    else:
        shiftB8_7 = 0,0
    
    shift4_5, error, diffphase = register_translation(A4, A5)
    if np.all(shift4_5==0):
        B4 = plt.imread(FilesMyxo_ROI4[i])
        B5 = plt.imread(FilesMyxo_ROI5[i])
        shiftB4_5, error, diffphase = register_translation(B4, B5)
    else:
        shiftB4_5 = 0,0
    
    shift5_6, error, diffphase = register_translation(A5, A6)
    if np.all(shift5_6==0):
        B5 = plt.imread(FilesMyxo_ROI5[i])
        B6 = plt.imread(FilesMyxo_ROI6[i])
        shiftB5_6, error, diffphase = register_translation(B5, B6)
    else:
        shiftB5_6 = 0,0
    
    shift3_2, error, diffphase = register_translation(A3, A2)
    if np.all(shift3_2==0):
        B3 = plt.imread(FilesMyxo_ROI3[i])
        B2 = plt.imread(FilesMyxo_ROI2[i])
        shiftB3_2, error, diffphase = register_translation(B3, B2)
    else:
        shiftB3_2 = 0,0
    
    shift2_1, error, diffphase = register_translation(A2, A1)
    if np.all(shift2_1==0):
        B2 = plt.imread(FilesMyxo_ROI2[i])
        B1 = plt.imread(FilesMyxo_ROI1[i])
        shiftB2_1, error, diffphase = register_translation(B2, B1)
    else:
        shiftB2_1 = 0,0
    
    # from left to right
    shift9_4, error, diffphase = register_translation(A9, A4)
    if np.all(shift9_4==0):
        B9 = plt.imread(FilesMyxo_ROI9[i])
        B4 = plt.imread(FilesMyxo_ROI4[i])
        shiftB9_4, error, diffphase = register_translation(B9, B4)
    else:
        shiftB9_4 = 0,0
    
    shift8_5, error, diffphase = register_translation(A8, A5)
    if np.all(shift8_5==0):
        B8 = plt.imread(FilesMyxo_ROI8[i])
        B5 = plt.imread(FilesMyxo_ROI5[i])
        shiftB8_5, error, diffphase = register_translation(B8, B5)
    else:
        shiftB8_5 = 0,0
    
    shift7_6, error, diffphase = register_translation(A7, A6)
    if np.all(shift2_1==0):
        B7 = plt.imread(FilesMyxo_ROI7[i])
        B6 = plt.imread(FilesMyxo_ROI6[i])
        shiftB7_6, error, diffphase = register_translation(B7, B6)
    else:
        shiftB7_6 = 0,0
    
    shift4_3, error, diffphase = register_translation(A4, A3)
    if np.all(shift4_3==0):
        B4 = plt.imread(FilesMyxo_ROI4[i])
        B3 = plt.imread(FilesMyxo_ROI3[i])
        shiftB4_3, error, diffphase = register_translation(B4, B3)
    else:
        shiftB4_3 = 0,0
    
    shift5_2, error, diffphase = register_translation(A5, A2)
    if np.all(shift5_2==0):
        B5 = plt.imread(FilesMyxo_ROI5[i])
        B2 = plt.imread(FilesMyxo_ROI2[i])
        shiftB5_2, error, diffphase = register_translation(B5, B2)
    else:
        shiftB5_2 = 0,0
    
    shift6_1, error, diffphase = register_translation(A6, A1)
    if np.all(shift6_1==0):
        B6 = plt.imread(FilesMyxo_ROI6[i])
        B1 = plt.imread(FilesMyxo_ROI1[i])
        shiftB6_1, error, diffphase = register_translation(B6, B1)
    else:
        shiftB6_1 = 0,0
    
    
    # Concatenate the arrays
    Total = np.array((shift9_8, shift8_7, shift4_5, shift5_6, shift3_2, shift2_1, shift9_4, shift4_3, shift8_5, shift5_2, shift7_6, shift6_1, shiftB9_8, shiftB8_7, shiftB4_5, shiftB5_6, shiftB3_2, shiftB2_1, shiftB9_4, shiftB4_3, shiftB8_5, shiftB5_2, shiftB7_6, shiftB6_1))
    # Create a filename which will change with each loop
    j = i+1
    name_i = '%03d' %j
    filename = str(new_dir_name)+'/CrossCorrelations_frame'+str(name_i)+'.txt'
    
    # Save the text file in the newly created '/Tiling/' directory
    np.savetxt(filename, Total)
    print(i)