#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 

@author: legall
"""




"""
Performs MSD analysis on single tracks.

Loads tracking analysis matlab files and performs MSD analysis to export 
diffusion and alpha exponent coefficients by fitting the 5 first points
of the MSD for each track.

Usage:
-  modify PreprocessingInfos.yaml file to specify the desired experiment to process
-  modify PreprocessingInfos.yaml file to specify the imsd input parameters to match your experimental conditions (pixel size and time between images)
-  modify DATA PATH in function "experimentDataPath" to match the location of you matlab mat files

"""



#%% IMPORTS ###################################################################
import os
import yaml
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

# Optionally, tweak styles.
mpl.rc('figure',  figsize=(10, 5))
mpl.rc('image', cmap='gray')

import numpy as np
import numpy.matlib
import pandas as pd
import trackpy as tp
from scipy.io import loadmat
import time

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
    # DATA PATH
    experiment = getInfoFromListOfDict(inputs, 'experiment')
    if experiment == 'exp8':
        # foldername = 'Myxo_Sara_Data/' 
        foldername = 'rawData_2020_Experiment_8_sara/'
    elif experiment == 'exp9':
        foldername = 'Myxo_Sara_Data_with_EC_exp9/'
    elif experiment == 'exp11':    
        foldername = 'Myxo_Sara_Data_Exp11/'
    elif experiment == 'exp12':
        foldername = 'Myxo_Sara_Data_Exp12/'
    elif experiment == 'expX':
        foldername = 'Myxo_Sara_Data_with_EC_expX/'
    elif experiment == 'expX_bis':
        foldername = 'Myxo_Sara_Data_with_EC_expX_bis/'
    elif experiment == 'exp10_DK':
        foldername = 'Myxo_Sara_Data_Exp10_DK/'
    elif experiment == 'exp50_wt':
        foldername = 'Experiment_50_Sara_WT/'
    elif experiment == 'exp53':
        foldername = 'Experiment_53_Sara_A+S-/'  
    elif experiment == 'exp55':
        foldername = 'Experiment_55_Sara_A-S+/'
    elif experiment == 'exp72':
        foldername = 'Experiment_72_Sara_AminusSplus/' 
        folderpath = '/mnt/tronador/AnalyzedData_2021/'+foldername 


    cbs = ("miramar","lopevi",'borabora')
    if "atlantis" in os.uname()[1]:
        folderpath = '/home/marcnol/data/myxo/napari/'+foldername
    elif os.uname()[1] in cbs:
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername
        # folderpath = '/mnt/moreno/rawData/2020/'+foldername
    else:
        folderpath = '/home/antoine/Documents/'+foldername
    
    
    return folderpath

def loadmatfiles(folderpath):
    '''
    Loads all matlab .mat files (matrices) within a folder into a single dictionnary. Each matrix is a feature

    :param folderpath: folder path that contains mat files
    :return: my_dict: dictionnary with each field containing one matrix
    '''
    print('\n *** Data loading ***')
    my_dict = {}
    for file in os.listdir(folderpath):
        if file.endswith(".mat"):

            data = loadmat(folderpath+file)
            filename, ext = os.path.splitext(file)
            if filename=='Total_Density':
                my_dict[filename] = np.log10(data.get(filename))
            else:
                my_dict[filename] = data.get(filename)
            print(filename, ' loaded')
    print('--> All files successfully loaded.')

    # Add a field with tracks ID to control results
    ID = np.arange(my_dict.get('Total_tracks').shape[0])  # Create an array with as many tracks ID as there are tracks
    ID = np.random.permutation(ID)  # randomize IDs so that proximal tracks won't show up with similar colors.
    ID = np.matlib.repmat(ID,  my_dict.get('Total_tracks').shape[1], 1) # reshape IDs to match tracks length
    ID = np.swapaxes(ID, 0, 1)  # swap row and col
    my_dict['TrackID'] = ID

    # Add a field with time stamp
    del ID
    ID = np.arange(my_dict.get('Total_tracks').shape[1])  # Create an array with as many tracks ID as there are tracks
    ID = np.matlib.repmat(ID, my_dict.get('Total_tracks').shape[0], 1)  # reshape IDs to match number of tracks
    my_dict['TimeStamp'] = ID

    # Add a field with tracks length
    track_length = np.count_nonzero(my_dict.get('Total_tracks'), axis=1)  # number of cells detected
    track_length = np.matlib.repmat(track_length, my_dict.get('Total_tracks').shape[1], 1)  # reshape IDs to match tracks length
    track_length = np.swapaxes(track_length, 0, 1)  # swap row and col
    my_dict['track_length'] = track_length
    print('--> New fields TrackID and TimeStamp successfully added \n')

    return my_dict

def build_trackpy_input_dict(data_dict, min_track_length = 15):
    '''
    Takes a dictionnary and convert relevant fields into a panda dataframe to feed trackpy.
    :param data_dict:
    :param min_track_length: criteria to select tracks with a minimum nb of frames
    :return: dataframe for trackpy
    '''
    t = {}

    # - Remove empty cells and s
    frame = data_dict['TimeStamp'].flatten()
    rows = data_dict['Total_Row'].flatten()
    cols = data_dict['Total_Col'].flatten()
    area = data_dict['Total_Area'].flatten()
    particle = data_dict['TrackID'].flatten()
    track_length = data_dict['track_length'].flatten()
    
    # Select cells with a min track length
    frame = frame[(track_length > min_track_length)]
    rows = rows[(track_length > min_track_length)]
    cols = cols[(track_length > min_track_length)]
    particle = particle[(track_length > min_track_length)]
    area = area[(track_length > min_track_length)]
    
    # remove cells with null area
    frame = frame[~(area == 0)]
    rows = rows[~(area == 0)]
    cols = cols[~(area == 0)]
    particle = particle[~(area == 0)] 

    t = pd.DataFrame(
        dict(particle=particle, x=cols, y=rows, frame=frame))

    print('--> Conversion done \n')
    return t

def fit_single_msd(im, t_id, n_lag, plot01 = False):
    '''
    Takes a pandas dataframe containing MSDs and fit them with a power law
    :param im: pandas dataframe with indices the lag times and column the MSDs values, each column is a differnt track MSD
    :param t_id: track ID
    :param n_lag: Number of lag time to fit
    :param plot01: plots results
    :return: diffusion coefficient and power law exponent
    '''
    tlag = im.index  
    msds_val = im[t_id].values 
    if plot01 :
        fig, ax = plt.subplots()
        ax.plot(tlag, msds_val, 'k-', alpha=0.1)  # black lines, semitransparent
        ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
               xlabel='lag time $t$')
        
    # Fit the curve with a power law
    pd_msd = pd.DataFrame({'msd': msds_val, 'lagt': tlag}, dtype=np.float64)
    pd_msd = pd_msd.dropna()
    pd_msd = pd_msd.set_index('lagt')['msd']
    n_lag = min(pd_msd.shape[0], n_lag)
    results = tp.utils.fit_powerlaw(pd_msd.iloc[:n_lag], plot=plot01)
    if plot01 :
        ax.set_xscale('linear')
        ax.set_yscale('linear')
        ax.set_xlim([0, pd_msd.index.values.max()])
        text = 'Alpha = %s and D (µm²/s) = %s' % (round(results.values[0, 0], 2), round(results.values[0, 1] / 4, 5))
        anchored_text = AnchoredText(text, loc=2)
        ax.add_artist(anchored_text)

    return results

def fitMSDs(im, inputs):
    D_coeff = []
    alpha = []
    i=1
    t_start = time.time()
    n_lag = getInfoFromListOfDict(inputs, 'n_lag') # nth first points to fit
    for k in im.keys():
        print('Fitting MSD of track ID %s, iteration %s/%s' % (k, i, len(im.keys())))
        results = fit_single_msd(im, t_id=k, n_lag=n_lag, plot01=False)
        D_coeff.append(results.values[0, 1]/4)
        alpha.append(results.values[0, 0])
        i+= 1
    print('Elapsed time, seconds = %s' %round(time.time()-t_start))
    
    return D_coeff, alpha

# ======================================================================================================================
#                                                   MAIN
# ======================================================================================================================

if __name__ == '__main__':
    
    scriptname = os.path.basename(sys.argv[0])
    
    # Reads the yaml document with scripts parameters
    yamlDoc = open("PreprocessingInfos.yaml", 'r')
    inputs, outputs = readYamlDoc(yamlDoc, scriptname)
    
    # Gets parameters to load and export data
    folderpath = experimentDataPath(inputs)  
    
    #  LOADS TRACKS DATA INTO DICTIONNARY 
    data_dict = loadmatfiles(folderpath)

    # #  Build XYT pandas
    min_track_length = getInfoFromListOfDict(inputs, 'min_track_length')
    f = build_trackpy_input_dict(data_dict, min_track_length = min_track_length)

    #  Mean Squared Displacement
    print('Computing MSDs, please wait...')
    t_start = time.time()
    pixsize = getInfoFromListOfDict(inputs, 'pixsize')
    dT = getInfoFromListOfDict(inputs, 'dT')
    max_lagtime = getInfoFromListOfDict(inputs, 'max_lagtime')
    im = tp.imsd(f,pixsize, dT, max_lagtime=max_lagtime) 
    print('--> MSDs computation done. Elapsed time = %s'% round(time.time()-t_start))
    
   
    # Fit all selected tracks
    D_coeff, alpha = fitMSDs(im, inputs)
    
    #%Export results
    D_coeff = np.array(D_coeff)
    D_coeff = np.log10(D_coeff)
    D_coeff_out = np.matlib.repmat(D_coeff,  data_dict.get('Total_tracks').shape[1], 1) 
    D_coeff_out = np.swapaxes(D_coeff_out, 0, 1)  # swap row and col
    
    alpha_out = np.matlib.repmat(alpha,  data_dict.get('Total_tracks').shape[1], 1) 
    alpha_out = np.swapaxes(alpha_out, 0, 1)  # swap row and col

    # - Save results
    saveResults = getInfoFromListOfDict(outputs, 'saveResults')
    if saveResults:
        outfile = folderpath+ 'MSD_D_5frames'
        np.save(outfile, D_coeff_out)
        outfile = folderpath+ 'MSD_n_5frames'
        np.save(outfile, alpha_out)