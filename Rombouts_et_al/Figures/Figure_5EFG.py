#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 11:03:49 2021

@author: legall
"""


''' 
Loads Kymograph fluorescence intensities and fit intensity profiles with exponentials.

Usage:
- modify DATA PATH to match the datashader maps location

'''


#%% IMPORTS ###################################################################
import os
import yaml
import sys
import matplotlib.pyplot as plt
import numpy as np
import pickle5 as pickle
from numpy import inf
from scipy import ndimage 

import scipy.optimize

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
    
    '''
    
    if experiment == 'exp6':
        foldername = 'Experiment_6_sara/' 
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/' 
        foldername = folderpath + foldername +"EC_consumption/"
    elif experiment == 'exp8':
        foldername = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data/EC_consumption/'
    elif experiment == 'exp9':
        foldername = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data_with_EC_exp9/EC_consumption/'
    elif experiment == 'exp11':
        foldername = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data_Exp11/EC_consumption/'
    elif experiment == 'exp12':
        foldername = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data_Exp12/EC_consumption/'
    elif experiment == 'expX_bis':
        foldername = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data_with_EC_expX_bis/EC_consumption/'
    elif experiment == 'exp10_DK':
        foldername = '/mnt/grey/DATA/users/Antoine/Myxo_Sara_Data_Exp10_DK/EC_consumption/'
    elif experiment == 'exp6_wt':
        foldername = '/mnt/grey/DATA/AnalyzedData_2021/Experiment_6_sara/EC_consumption/'
    elif experiment == 'exp53':
        foldername = '/mnt/grey/DATA/AnalyzedData_2021/Experiment_53_Sara_A+S-/EC_consumption/'   
    elif experiment == 'exp54':
        foldername = '/mnt/grey/DATA/AnalyzedData_2021/Experiment_54_Sara_A-S+/EC_consumption/'   
    elif experiment == 'exp55':
        foldername = '/mnt/grey/DATA/AnalyzedData_2021/Experiment_55_Sara_A-S+/EC_consumption/'  
    elif experiment == 'exp50':
        foldername = '/mnt/grey/DATA/AnalyzedData_2021/Experiment_50_Sara_WT/EC_consumption/'  
    
    
    if experiment == 'expX' or experiment == 'exp10_DK' :
        folderpath = ''+foldername
    else:
        folderpath = '/mnt/grey/DATA'+foldername
        
    return folderpath, foldername


def loadKymographIntensities(foldername, inputs):
    Chan = getInfoFromListOfDict(inputs, 'Chan')
    experiment = getInfoFromListOfDict(inputs, 'experiment')
    kymograph = []
    offset = [] 
    print('--> Loading file, please wait ...\n')
    for ROI in ROIs:
        file_to_read = open(foldername+ experiment+'_Chan' + str(Chan)+'_ROI'+str(ROI)+ '.pkl',  'rb')
        data = pickle.load(file_to_read)
        kymograph.append(data['kymograph'])
        offset.append(data['offset'])
    kymograph = np.sum(np.array(kymograph), axis=0)
    offset = np.sum(np.array(offset), axis=0)
    print('--> File loaded. \n')
    
    return kymograph, offset

def processKymograph(kymograph, offset):
    # =============================================================================
    #      Get image offset, smooth it, and correct kymograph
    # =============================================================================
    
    offset_filt = ndimage.median_filter(offset, size=31)
    kymograph_corr = kymograph.astype('float') - offset_filt.transpose()

    if i==0:
        c = 'b'
    elif i==1:
        c = 'orange'
    else:
        c = 'g'

    # =============================================================================
    #      Fit Gaussian to background
    # =============================================================================
    def gaus(x,a,x0,sigma,off):
        return a*np.exp(-(x-x0)**2/(2*sigma**2)) + off

    y = kymograph_corr[:,nb_frames-1]
    y = ndimage.median_filter(y, size=501)
    x = range(y.shape[0])
    
    popt = np.array([3.1815218 , 1031.1085082 ,  665.78067835, 0])
    kymograph_corr2 = kymograph_corr / (gaus(x,*popt)[:, None]+0)

    # =============================================================================
    #     Normalize kymograph with profile at t=0
    # =============================================================================
    temp = kymograph_corr2.mean(axis=1)>kymograph_corr2.mean()
    kymograph_corr_mask = np.repeat(temp[:,np.newaxis], nb_frames, axis=1)
    kymograph_corr_masked = kymograph_corr2*kymograph_corr_mask
    kymograph_corr_masked = kymograph_corr_masked*(kymograph_corr_masked>0)
    kymograph_norm = kymograph_corr_masked/kymograph_corr_masked[:,10][:,None]
    
    KProj = kymograph_norm
    
    # Do a weighted average ignoring NaN and Inf values
    KProj[KProj == -inf] = np.NaN
    KProj[KProj == inf] = np.NaN
    masked_data = np.ma.masked_array(KProj, np.isnan(KProj))
    weights=kymograph_corr[:,0]
    weights = weights*(weights>weights.max()/2)
    average = np.ma.average(masked_data, axis=0, weights=weights)
    KProj = average.filled(np.nan)

    fig = plt.figure(101, figsize=(4,4))
    if (i == 0) or (i==2):
        lines, = plt.plot(KProj, ".",  markersize=2)
        
        plt.ylim(0,1.1)
        plt.xlim(10,700)
        plt.ylabel('EC density, norm', fontsize=20)
        plt.xlabel('Time, frame number', fontsize=20)
        plt.legend()
        plt.title(experiment, fontsize=24)
        plt.tight_layout()
        
        x0,x1 = plt.gca().get_xlim()
        y0,y1 = plt.gca().get_ylim()
        plt.gca().set_aspect(abs(x1-x0)/abs(y1-y0))
        
        if fit_exp:
            xs = np.arange(KProj.shape[0])
            params, cv = scipy.optimize.curve_fit(monoExp, xs, KProj, bounds = ([0.9, 0, 0], [1.1, 10**10, 0.2]))
            m, t, b = params
            # plot the results
            plt.plot(xs, monoExp(xs, m, t, b), c=lines.get_color(), lw=1, label= labels[i] +' : a=%5.2f, b=%5.0f, c=%5.2f' % tuple(params))
                            
            # inspect the parameters
            print(f"Y = {m} * e^(-{t} * x) + {b}")
            print(f"Tau = {t } frames")
            plt.legend( prop={'size': 6})
    

def monoExp(x, A, F, O):
    return A * np.exp(-x / F) + O

# =============================================================================
#                                   MAIN
# =============================================================================
plt.close('all')
if __name__ == '__main__':
    
    scriptname = os.path.basename(sys.argv[0])
    
    # Reads the yaml document with scripts parameters
    yamlDoc = open("ScriptsInfos.yaml", 'r')
    inputs, _ = readYamlDoc(yamlDoc, scriptname)
    
    # Gets parameters to load and export data
    experiment = getInfoFromListOfDict(inputs, 'experiment')
    _, foldername = experimentDataPath(experiment)


    fit_exp = getInfoFromListOfDict(inputs, 'fit_exp') # Fit exponential decay
    labels = ['Bottom', 'Middle', 'Top']
    
    for i in np.arange(3):
        # =========================================================================
        #     Parameters for ROI and Channel
        # =========================================================================
        ROIs = getInfoFromListOfDict(inputs, 'ROIs')[i]
        Chan = getInfoFromListOfDict(inputs, 'Chan')
        nb_frames = getInfoFromListOfDict(inputs, 'nb_frames')
        
        # =========================================================================
        #      Load data
        # =========================================================================
        kymograph, offset = loadKymographIntensities(foldername, inputs)
        
        # =========================================================================
        #      Process data
        # =========================================================================
        processKymograph(kymograph, offset)
        
       
