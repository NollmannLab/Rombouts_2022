#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 15:59:04 2021

@author: legall
"""


"""
Performs bacterial cluster segmentation.

1- Reads Matlab mat files containing connectivity & pixelsIds lists of segmented cells,
2- applies a dilation with a 10x10 kernel (for pixels of 105 nm) 
3- Labels the resulted dilated image
4- Computing clusters and single cell properties (from dilated and original label data)
5- Exports results in dataframes

Usage:
-  modify PreprocessingInfos.yaml file to specify the desired experiment to process
-  modify PreprocessingInfos.yaml file to specify the kernel size to match your experimental conditions (pixel size)
-  modify DATA PATH in function "experimentDataPath" to match the location of you matlab mat files

"""




#%% IMPORTS ###################################################################
import os
import yaml
import sys
from skimage import measure
import numpy as np
import cv2
import time
import pandas as pd
import h5py


# ======================================================================================================================
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


def experimentDataPath(inputs):
    '''
    

    Parameters
    ----------
    experiment : TYPE
        DESCRIPTION.

    Returns
    -------
    folderpath : TYPE
        DESCRIPTION.
    
    '''
    experiment = getInfoFromListOfDict(inputs, 'experiment')
    if experiment == 'exp8':
        foldername = '/rawData_2020/Experiment_8_sara/Analyzed/Tiling_Drift_PostProcess/'
        exportpath = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data/dataframe/'
    elif experiment == 'exp9':
        foldername = '/rawData_2020/Experiment_9_sara/Segmentation/Analyzed/Tiling_Drift_PostProcess/'
        exportpath = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data_with_EC_exp9/dataframe/'
    elif experiment == 'exp11':
        foldername = '/rawData_2021/Experiment_11_FTL_A-S+_20210225/004_FastTimeLapse_RAMM_Test/Analyzed/Tiling_Drift_PostProcess/'
        exportpath = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data_Exp11/dataframe/'
    elif experiment == 'exp12':
        foldername = '/rawData_2021/Experiment_12_FTL_A+S-_20210226/001_FastTimeLapse_RAMM_Test/Analyzed/Tiling_Drift_PostProcess/'
        exportpath = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data_Exp12/dataframe/'
    elif experiment == 'expX':
        foldername = '/mnt/PALM_dataserv/DATA/Sara/DATA/TimeLapseData/PredationAssays/2020_05_20-SampleDilutedEcoliColony/002_FastTimeLapse_RAMM_Test/Analyzed/Tiling_Drift_PostProcess/'
        exportpath = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data_with_EC_expX/dataframe/'
    elif experiment == 'exp10_DK':
        foldername = '/mnt/grey/DATA/rawData_2021/Experiment_10_FTL_DelteKill_20210219/004_FastTimeLapse_RAMM_Test/Analyzed/Tiling_Drift_PostProcess/'
        exportpath = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data_Exp10_DK/dataframe/'
    elif experiment == 'exp50_wt':
        foldername = '/mnt/grey/DATA/rawData_2021/Experiment_50_Sara_WT/003_FastTimeLapse_RAMM_Test/Analyzed/Tiling_Drift_PostProcess/'
        exportpath = '/mnt/grey/DATA/AnalyzedData_2021/Experiment_50_Sara_WT/'

    if experiment == 'expX' or experiment == 'exp10_DK' or experiment == 'exp50_wt':
        folderpath = ''+foldername
    else:
        folderpath = '/mnt/grey/DATA'+foldername
    
    
    
    return folderpath, exportpath


def load_labelmatrix_matfiles(folderpath,filename):
    '''
    Read Matlab mat files containing connectivity & pixelsIds lists 
    :param folderpath: folder path that contains the mat file
    :return: LB: Label image
    '''
    print('\n *** ' + filename + ' Data loading ***')
    t0 = time.time()
    f = h5py.File(folderpath + filename + '.mat', "r")
    list(f.keys())
    grp = f[filename]    # header to the 
    X=list(grp.items())     # List of fields in 
    
    ImageSize = X[1][1][:]
    PixelIdxList = X[3][1][:]
    
    nrow = ImageSize[0][0].astype('int64')
    ncol = ImageSize[1][0].astype('int64')
    
    t1 = time.time()
    image_out = np.zeros((nrow*ncol,1))
    for it in range(len(PixelIdxList)):
        pixels = f[PixelIdxList[it][0]][:].astype('int64')
        image_out[list(pixels[0][:].flatten())] = it
    t2 = time.time()    
    LB= np.reshape(image_out, (nrow, ncol)).T
    t3 = time.time()
    
    # Timers
    print('Loading time : ', t1-t0)
    print('Conversion time : ', t2-t1)
    print('Reshape time : ', t3-t2)
    
    print(filename + ' loaded ***')
    return LB


def clusterMasks(label_image, Frame_num, inputs):
    # Performing dilation on the mask
    print('Dilating mask...')
    kernel_size = getInfoFromListOfDict(inputs, 'kernel_size')
    kernel=np.ones((kernel_size, kernel_size), np.uint8)     # Creating a 3x3 kernel for BW mask dilation
    BW_dilated=cv2.dilate(np.uint8((label_image > 0.5).astype(np.int_)), kernel)
    print('Mask dilated.')
    
    # Mask labeling
    print('Labeling mask...')
    Clust_labels = measure.label(BW_dilated, background=0)
    
    # Compute clusters properties
    print('Computing clusters props...')
    clusters_props = measure.regionprops_table(Clust_labels, intensity_image=None ,
                                       properties=['label', 
                                                   'centroid', 
                                                   'area', 
                                                   'perimeter',
                                                   'eccentricity',
                                                   'solidity',
                                                   'orientation'], cache = True)
    clusters_df = pd.DataFrame(clusters_props)
    clusters_df = clusters_df.rename(columns = {'label':'cluster_id'})
    clusters_df['Frame_num']= Frame_num
    
    # Compute cells properties
    print('Computing cells props...')
    cells_props = measure.regionprops_table(label_image, Clust_labels, 
                                            properties=['label', 
                                                        'centroid', 
                                                        'area',
                                                        'perimeter', 
                                                        'major_axis_length', 
                                                        'eccentricity',
                                                        'solidity',
                                                        'orientation', 
                                                        'mean_intensity'], cache = True)
    
    cells_df = pd.DataFrame(cells_props)
    cells_df = cells_df.rename(columns = {'label':'cell_id','mean_intensity':'cluster_id'})
    cells_df['cluster_id'] = cells_df['cluster_id'].astype('uint64')
    cells_df['Frame_num']= Frame_num
    
    # Count nb of cells per cluster
    serie_counts = cells_df['cluster_id'].value_counts(ascending=True)
    cells_df['Nb_cells']= cells_df['cluster_id']
    for val in cells_df['cluster_id']:
            cells_df.loc[cells_df['cluster_id']==val, 'Nb_cells'] = serie_counts[val]
    
    return clusters_df, cells_df

def multiscaleSegmentation(inputs, folderpath):
    
    t0 = time.time()
    appended_cells_data = []
    appended_clust_data = []
    nb_frames = getInfoFromListOfDict(inputs, 'nb_frames')
    for i in range(1,nb_frames):
        filename = 'Frame_' + str(i).zfill(3)
        
        # Load Matlab labeled image
        label_image =  load_labelmatrix_matfiles(folderpath, filename).astype('uint64')
        
        # Dilate masks and segment them again
        clusters_df, cells_df = clusterMasks(label_image, i, inputs)
        
        # Gather results
        appended_cells_data.append(cells_df)
        appended_clust_data.append(clusters_df)
    
        del cells_df, clusters_df
        
    cells_df = pd.concat(appended_cells_data, ignore_index=True)
    clusters_df = pd.concat(appended_clust_data, ignore_index=True)
    t1 = time.time()
    print('\n Total computing time : ', t1-t0)
    return clusters_df, cells_df
# ======================================================================================================================
#                                                   MAIN
# ======================================================================================================================

if __name__ == '__main__':
    
    scriptname = os.path.basename(sys.argv[0])
    
    # Reads the yaml document with scripts parameters
    yamlDoc = open("PreprocessingInfos.yaml", 'r')
    inputs, outputs = readYamlDoc(yamlDoc, scriptname)
    
    # Gets parameters to load and export data
    folderpath, exportpath = experimentDataPath(inputs)
    
    # Read Matlab mat files containing connectivity & pixelsIds lists
    clusters_df, cells_df = multiscaleSegmentation(inputs, folderpath)
    
    # Export results
    saveResults = getInfoFromListOfDict(outputs, 'saveResults')
    if saveResults:
        print('Savign')
        clusters_df.to_pickle(exportpath+"clusters_df.pkl")
        cells_df.to_pickle(exportpath+"cells_df.pkl")
    
  