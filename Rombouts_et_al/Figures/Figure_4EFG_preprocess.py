
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 10:16:43 2021

@author: legall
"""


''' 
Loads raw fluorescence images and extract the fluorescence intensity vs time.

Usage:
- modify DATA PATH to match the locations of raw files of the experiment
- 
'''
#%% IMPORTS ###################################################################
import os
import yaml
import sys
import numpy as np
import time
from PIL import Image
import pickle5 as pickle

try:
    from skimage.feature import register_translation as phase_cross_correlation
except:
    from skimage.registration import phase_cross_correlation

from skimage.filters import threshold_triangle

#%% ====================================================================================================================
#                                                   FUNCTIONS
# ======================================================================================================================

def getInfoFromListOfDict(ListOfDict, field, default = False):
        '''
        Function to access a dictionnary value directly from a list of dictionnaries

        Parameters
        ----------
        ListOfDict : List of dictionnaries
        field : string of the dictionnary key to access value
        default : returns the default value within the dictionnary

        Returns
        -------
        val : value of the dictionnary field         

        '''
        if default:
            val = next(list(item.values())[1] for item in ListOfDict if list(item.keys())[0] == field)
        else:
            val = next(list(item.values())[0] for item in ListOfDict if list(item.keys())[0] == field)
        return val

def readYamlDoc(yamlDoc, scriptname):
    '''
    Parameters
    ----------
    yamlDoc : Yaml document with scripts information, including parameters

    Returns
    -------
    inputs : inputs parameters for the script
    outputs : outputs parameters for the script

    '''
    dictionary = yaml.load(yamlDoc, Loader=yaml.FullLoader)

    print("Reading scripts parameters from Yaml document:")
    for script, scriptinfo in dictionary.items():  
        if script == scriptname:
            print('=====  ' + script +' info   =====')
            for item in range(len(scriptinfo)):
                for key, value in scriptinfo[item].items() : #for each script
                    print("--> "+ key + " : " + str(value))
                    if key == "inputs":
                        inputs = value
                    if key == "outputs":
                        outputs = value
    return inputs, outputs


def experimentDataPath(experiment):
    '''
    

    Parameters
    ----------
    experiment : TYPE
        DESCRIPTION.

    Returns
    -------
    folderpath : TYPE
        DESCRIPTION.
    exportpath : TYPE
        DESCRIPTION.

    '''
    if experiment == 'exp8':
        foldername = '/rawData_2020/Experiment_8_sara/Ch_2/'
        exportpath = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data/EC_consumption/'
    elif experiment == 'exp9':
        foldername = '/rawData_2020/Experiment_9_sara/Ch_1/'
        exportpath = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data_with_EC_exp9/EC_consumption/'
    elif experiment == 'exp11':
        foldername = '/rawData_2021/Experiment_11_FTL_A-S+_20210225/004_FastTimeLapse_RAMM_Test/Ch_2/'
        exportpath = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data_Exp11/EC_consumption/'
    elif experiment == 'exp12':
        foldername = '/rawData_2021/Experiment_12_FTL_A+S-_20210226/001_FastTimeLapse_RAMM_Test/Ch_2/'
        exportpath = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data_Exp12/EC_consumption/'
    elif experiment == 'expX_bis':
        foldername = '/mnt/PALM_dataserv/DATA/Sara/DATA/TimeLapseData/PredationAssays/2020_05_20-SampleDilutedEcoliColony/002_FastTimeLapse_RAMM_Test/Ch_2/'
        exportpath = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data_with_EC_expX_bis/EC_consumption/'
    elif experiment == 'exp10_DK':
        foldername = '/mnt/grey/DATA/rawData_2021/Experiment_10_FTL_DelteKill_20210219/004_FastTimeLapse_RAMM_Test/Ch_2/'
        exportpath = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data_Exp10_DK/EC_consumption/'
    elif experiment == 'exp6_wt':
        foldername = '/mnt/grey/DATA/rawData_2020/Experiment_6_sara/tiff_raw_data/Ch_1/'
        exportpath = '/mnt/grey/DATA/AnalyzedData_2021/Experiment_6_sara/EC_consumption/'
    elif experiment == 'exp53':
        foldername = '/mnt/grey/DATA/rawData_2021/Experiment_53_Sara_A+S-/001_FastTimeLapse_RAMM_Test/Ch_2/'
        exportpath = '/mnt/grey/DATA/AnalyzedData_2021/Experiment_53_Sara_A+S-/EC_consumption/'
    elif experiment == 'exp54':
        foldername = '/mnt/grey/DATA/rawData_2021/Experiment_54_Sara_A-S+/002_FastTimeLapse_RAMM_Test/Ch_2/'
        exportpath = '/mnt/grey/DATA/AnalyzedData_2021/Experiment_54_Sara_A-S+/EC_consumption/'
    elif experiment == 'exp55':
        foldername = '/mnt/grey/DATA/rawData_2021/Experiment_55_Sara_A-S+/001_FastTimeLapse_RAMM_Test/Ch_2/'
        exportpath = '/mnt/grey/DATA/AnalyzedData_2021/Experiment_55_Sara_A-S+/EC_consumption/'        
    elif experiment == 'exp50':
        foldername = '/mnt/grey/DATA/rawData_2021/Experiment_50_Sara_WT/003_FastTimeLapse_RAMM_Test/Ch_2/'
        exportpath = '/mnt/grey/DATA/AnalyzedData_2021/Experiment_50_Sara_WT/EC_consumption/'         
    elif experiment == 'exp51':
        foldername = '/mnt/grey/DATA/rawData_2021/Experiment_51_Sara_A+S-/001_FastTimeLapse_RAMM_Test/Ch_2/'
        exportpath = '/mnt/grey/DATA/AnalyzedData_2021/Experiment_51_A+S-/EC_consumption/'    

        
    if experiment == 'expX' or experiment == 'expX_bis' or experiment == 'exp10_DK' or experiment == 'exp51' or experiment == 'exp55':
        folderpath = ''+foldername
    else:
        folderpath = '/mnt/grey/DATA'+foldername
        
    return folderpath, exportpath

def processROI(inputs, folderpath, ROI):
    '''
    

    Parameters
    ----------
    folderpath : TYPE
        DESCRIPTION.
    ROI : TYPE
        DESCRIPTION.

    Returns
    -------
    thresh : TYPE
        DESCRIPTION.
    offset : TYPE
        DESCRIPTION.
    proj_bck : TYPE
        DESCRIPTION.
    proj : TYPE
        DESCRIPTION.
    kymograph : TYPE
        DESCRIPTION.

    '''
    Chan = getInfoFromListOfDict(inputs, 'Chan')
    nb_frames = getInfoFromListOfDict(inputs, 'nb_frames')
    z_focal_plane = getInfoFromListOfDict(inputs, 'z_focal_plane')
    
    # =============================================================================
    #     Open fisrt image to use it as reference (drift)
    # =============================================================================
    pil_img = Image.open(folderpath + '000__FTL_Ch_'+str(Chan)+'_ROI_'+str(ROI)+'.tif')
    pil_img.seek(z_focal_plane)
    EC_000 = np.array(pil_img)
    
    thresh = threshold_triangle(EC_000.astype('float'))
    
    # =============================================================================
    #     Loop through the movies to extract intensity profile
    # =============================================================================
    kymograph = np.zeros((EC_000.shape[0],nb_frames))
    proj_bck = np.zeros((EC_000.shape[0],nb_frames))
    offset = np.zeros((nb_frames,1))
    for i in range(0,nb_frames):
        print('Opening image', i)
        pil_img = Image.open(folderpath + str(i).zfill(3)+ '__FTL_Ch_'+str(Chan)+'_ROI_'+str(ROI)+'.tif')
        if i<50 :
            pil_img.seek(z_focal_plane)
        else:
            pil_img.seek(z_focal_plane)
        EC =  np.array(pil_img) 
        
        # Correct drift using reference image
        shift, error, diffphase = phase_cross_correlation(EC_000, EC)
        EC = np.roll(EC, shift.astype(int), axis=(0, 1))
        
        # Threshold signal to get an estimation of the image background baseline
        temp = EC*(EC<thresh)
        temp = temp.astype('float')
        temp[temp == 0] = 'nan'
        offset[i] = np.nanmean(temp)
        proj_bck[:,i] = np.nanmean(temp.astype('float'), axis=1)

        # Compute projection
        proj = np.nanmean(EC.astype('float'), axis=1)
        kymograph[:,i] = proj
        
    return thresh, offset, proj_bck, proj, kymograph

def dataExport(exportpath, inputs, thresh, offset, proj_bck, proj, kymograph):
    '''
    

    Parameters
    ----------
    exportpath : TYPE
        DESCRIPTION.
    inputs : TYPE
        DESCRIPTION.
    thresh : TYPE
        DESCRIPTION.
    offset : TYPE
        DESCRIPTION.
    proj_bck : TYPE
        DESCRIPTION.
    proj : TYPE
        DESCRIPTION.
    kymograph : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    print('--> Saving file, please wait ...\n')
    Chan = getInfoFromListOfDict(inputs, 'Chan')
    dict_to_export = dict(kymograph=kymograph, proj=proj, offset=offset, proj_bck= proj_bck, thresh=thresh, ROI =ROI, Chan = Chan, experiment=experiment)
    file_to_write = open(exportpath+ experiment+'_Chan' + str(Chan)+'_ROI'+str(ROI)+ '.pkl', "wb")
    pickle.dump(dict_to_export, file_to_write)
    file_to_write.close()
    print('--> Saving done \n')
    
# ======================================================================================================================
#                                                   MAIN
# ======================================================================================================================

if __name__ == '__main__':
 
    scriptname = os.path.basename(sys.argv[0])
    
    # Reads the yaml document with scripts parameters
    yamlDoc = open("ScriptsInfos.yaml", 'r')
    inputs, _ = readYamlDoc(yamlDoc, scriptname)
    
    # Gets parameters to load and export data
    experiment = getInfoFromListOfDict(inputs, 'experiment')
    folderpath, exportpath = experimentDataPath(experiment)
    
    # Process each ROI to get intensity vs time
    ROIs = np.arange(9)+1
    t_start = time.time()
    for ROI in ROIs:
        # ------ Process ROI
        thresh, offset, proj_bck, proj, kymograph = processROI(inputs, folderpath, ROI)
        # ------ Save output data
        dataExport(exportpath, inputs, thresh, offset, proj_bck, proj, kymograph)
    
    print('--> Elapsed time = %s'% round(time.time()-t_start))
    
