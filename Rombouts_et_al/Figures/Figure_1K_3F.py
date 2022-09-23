#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 09:53:12 2022

@author: legall
"""



''' 
Loads dataframe files with all features for all experiments and computes the mean
histograms of speed for all tracks.

Usage:
- modify DATA PATH to match the dataframe files of the experiments
- provide the desired experiment names in the lists corresponding to each 
condition (exp_WT, exp_As, exp_aS)
- make sure you provide the correct Exp_dT (time between images for each exp) 
and pix_size values (same value for all expierments)
- make sure to select the proper parameters in "Filter data range" section 
depending on whether you want to plot histograms for scouts/swarms/everybody...
'''


#%% IMPORTS ###################################################################
import os
import yaml
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle5 as pickle
import  matplotlib as mpl
import time

#%% FIFURE PROPS ###################################################

# mpl.rcParams.update(mpl.rcParamsDefault)
#### FIGURE
mpl.rcParams["figure.figsize"] = (10, 6)
mpl.rcParams["savefig.bbox"] = "tight"

#### FONT
mpl.rcParams["font.size"] = 24
mpl.rcParams['font.family'] = 'Arial'

#### AXES
mpl.rcParams["axes.grid"] = False
mpl.rcParams["axes.linewidth"] = 2

#### PLOTS
mpl.rcParams['lines.linewidth'] = 4

#### TICKS
mpl.rcParams["xtick.top"] = False
mpl.rcParams["xtick.major.width"] = 2
mpl.rcParams["xtick.major.size"] = 10
mpl.rcParams["xtick.minor.width"] = 2
mpl.rcParams["xtick.minor.size"] = 5
mpl.rcParams["ytick.right"] = False
mpl.rcParams["ytick.major.width"] = 2
mpl.rcParams["ytick.major.size"] = 10
mpl.rcParams["ytick.minor.width"] = 2
mpl.rcParams["ytick.minor.size"] = 5

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

        
def load_dataframe_exp(experiment = 'exp8'):
    # LOADS TRACKS DATA
    
    data = []
    if experiment == 'expX':
        foldername = 'Myxo_Sara_Data_with_EC_expX/'
        folderpath = '/mnt/grey/DATA/users/Antoine/'+foldername 
        fh = folderpath+"dataframe/all_features.pkl"
    elif experiment == 'expX_bis':
        foldername = 'Myxo_Sara_Data_with_EC_expX_bis/'
        folderpath = '/mnt/grey/DATA/users/Antoine/'+foldername 
        fh = folderpath+"dataframe/all_features_win_size_5.pkl"    
    elif experiment == 'exp6':
        foldername = 'Experiment_6_sara/' 
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername 
        fh = folderpath+"dataframe/all_features_WithNAN_win_size_5.pkl"
    elif experiment == 'exp8':
        foldername = 'Myxo_Sara_Data/' 
        folderpath = '/mnt/grey/DATA/users/Antoine/'+foldername 
        fh = folderpath+"dataframe/all_features_win_size_5.pkl"
    elif experiment == 'exp9':
        foldername = 'Myxo_Sara_Data_with_EC_exp9/'
        folderpath = '/mnt/grey/DATA/users/Antoine/'+foldername 
        fh = folderpath+"dataframe/all_features.pkl"
    elif experiment == 'exp10_DK':
        foldername = 'Myxo_Sara_Data_Exp10_DK/'
        folderpath = '/mnt/grey/DATA/users/Antoine/'+foldername 
        fh = folderpath+"dataframe/all_features_win_size_5.pkl"
    elif experiment == 'exp11':
        foldername = 'Myxo_Sara_Data_Exp11/'
        folderpath = '/mnt/grey/DATA/users/Antoine/'+foldername 
        fh = folderpath+"dataframe/all_features_win_size_5.pkl"
    elif experiment == 'exp12':
        foldername = 'Myxo_Sara_Data_Exp12/'
        folderpath = '/mnt/grey/DATA/users/Antoine/'+foldername 
        fh = folderpath+"dataframe/all_features.pkl"
    elif experiment == 'exp50':
        foldername = 'Experiment_50_Sara_WT/'
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername 
        fh = folderpath+"dataframe/all_features_WithNAN_win_size_5.pkl"
    elif experiment == 'exp51':
        foldername = 'Experiment_51_Sara_A+S-/'
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername 
        fh = folderpath+"dataframe/all_features_WithNAN_win_size_5.pkl"
    elif experiment == 'exp53':
        foldername = 'Experiment_53_Sara_A+S-/'  
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername 
        fh = folderpath+"dataframe/all_features_WithNAN_win_size_5.pkl"
    elif experiment == 'exp54':
        foldername = 'Experiment_54_Sara_A-S+/'  
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername 
        fh = folderpath+"dataframe/all_features_WithNAN_win_size_5.pkl"             
    elif experiment == 'exp55':
        foldername = 'Experiment_55_Sara_A-S+/'   
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername 
        fh = folderpath+"dataframe/all_features_WithNAN_win_size_5.pkl"                     
    elif experiment == 'exp66':
        foldername = 'Experiment_66_Anna_WT_SideOfECColony/' 
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername 
        fh = folderpath+"dataframe/all_features_WithNAN_win_size_5.pkl"
    elif experiment == 'exp70':
        foldername = 'Experiment_70_Sara_WT_cytosolicGFP/' 
        folderpath = '/mnt/tronador/AnalyzedData_2021/'+foldername 
        fh = folderpath+"dataframe/all_features_WithNAN_win_size_5.pkl"
    elif experiment == 'exp71':
        foldername = 'Experiment_71_A+S-/' 
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername 
        fh = folderpath+"dataframe/all_features_WithNAN_win_size_5.pkl"    
    elif experiment == 'exp72':
        foldername = 'Experiment_72_Sara_AminusSplus/' 
        folderpath = '/mnt/tronador/AnalyzedData_2021/'+foldername 
        fh = folderpath+"dataframe/all_features_WithNAN_win_size_5.pkl"
    else:
        print()
    
    
    with open(fh, "rb") as fh:
      data.append(pickle.load(fh))
    
    data = pd.concat(data, axis=1)     
    print('\n *** Tracks data for', experiment ,'successfully loaded ***\n')
    return data


def formatDataframe(data):
    # -- Compute max travelled distances
    # Get each track XYT coordiantes
    df_test = data.loc[:,['Total_Col', 'Total_Row','TimeStamp', 'TrackID']]
    # Get first and last position of each track
    df_out = df_test.groupby('TrackID').agg(['first', 'last'])#.stack()
    # Compute euclidian distances        
    df = np.sqrt( (df_out[('Total_Col', 'last')]-df_out[('Total_Col', 'first')])**2 + (df_out[('Total_Row', 'last')]-df_out[('Total_Row', 'first')])**2)
    df = df.reset_index(name='distance')
    df['distance_norm']  = df['distance'] / (df_out[('TimeStamp', 'last')] - df_out[('TimeStamp', 'first')])
    # Put the results back into the input dataframe
    data = data.merge(df[['TrackID', 'distance','distance_norm']],
                                  left_on=['TrackID'], 
                                  right_on=['TrackID'], 
                                  how='left')
    # -- Modify fields
    data['clust_area_x_mean'] = data.groupby(['TrackID'])['clust_area_x'].transform('mean') 
    data['instant_speed_mean2'] = data.groupby(['TrackID'])['instant_speed'].transform('mean') 
    fields_to_log = ['clust_area_x','clust_area_y','clust_area_x_mean','Total_Density_mean']
    data[fields_to_log] = np.log10(data[fields_to_log])
    fields_to_10 = ['instant_speed','instant_speed_from', 'instant_speed_to','instant_speed_mean', 'instant_speed_min', 'instant_speed_max', 'instant_speed_std']
    data[fields_to_10] = 10**(data[fields_to_10])
    return data

def filterDataRange(data, inputs):
    # -- Filter data range
    density_min = getInfoFromListOfDict(inputs, 'density_min')
    density_max = getInfoFromListOfDict(inputs, 'density_max')
    Nb_cells_min = getInfoFromListOfDict(inputs, 'Nb_cells_min')
    Nb_cells_max = getInfoFromListOfDict(inputs, 'Nb_cells_max')
    speed_min = getInfoFromListOfDict(inputs, 'speed_min')
    speed_max = getInfoFromListOfDict(inputs, 'speed_max')
    distance_norm_min = getInfoFromListOfDict(inputs, 'distance_norm_min')
    distance_norm_max = getInfoFromListOfDict(inputs, 'distance_norm_max')
    Mask_width = getInfoFromListOfDict(inputs, 'Mask_width')
    t1 = getInfoFromListOfDict(inputs, 't1')
    t2 = getInfoFromListOfDict(inputs, 't2')
    
    
    ind = data.index[ (data['Total_Density'] >= density_min) & (data['Total_Density'] <= density_max)
                     & (data['Nb_cells'] >= Nb_cells_min) & (data['Nb_cells'] <= Nb_cells_max)
                     & (data['instant_speed'] >= speed_min) & (data['instant_speed'] <= speed_max)
                     & (data['distance_norm'] >= distance_norm_min) & (data['distance_norm'] <= distance_norm_max)
                     & (data['track_length'] >= 0) & (data['track_length'] <= 70000)
                     & (data['instant_speed'] >= 0) & (data['instant_speed'] <= 70000)
                     & (data['Total_Col'] >= Mask_width) & (data['Total_Col'] <= 5880- Mask_width)
                     & (data['Total_Row'] >= Mask_width) & (data['Total_Row'] <= 5880- Mask_width)
                     & (data['TimeStamp'] >= t1) & (data['TimeStamp'] <= t2)].tolist()
    return ind

def loadAllDataframes(inputs):
    exp_WT = getInfoFromListOfDict(inputs, 'exp_WT')
    exp_aS = getInfoFromListOfDict(inputs, 'exp_aS')
    exp_As = getInfoFromListOfDict(inputs, 'exp_As')
    
    strains = [exp_WT,exp_aS,exp_As]
    
    data = [] 
    TrackID_max = 0
    for strain in range(len(strains)):
        for experiment in strains[strain]:
        
            # LOADS TRACKS DATA
            data_temp = load_dataframe_exp(experiment)
            data_temp = formatDataframe(data_temp)
            data_temp['TrackID'] = data_temp['TrackID'] + TrackID_max
            data_temp['Experiment'] = experiment
            TrackID_max = data_temp['TrackID'].max()
            data.append(data_temp)  
    
    data = pd.concat(data, axis=0, ignore_index=True) 
    return data

def extractSpeedFromData(data, inputs):
    exp_WT = getInfoFromListOfDict(inputs, 'exp_WT')
    exp_aS = getInfoFromListOfDict(inputs, 'exp_aS')
    exp_As = getInfoFromListOfDict(inputs, 'exp_As')
    
    strains = [exp_WT, exp_aS, exp_As]
    
    df_speed = pd.DataFrame()
    for strain in range(len(strains)):
        for experiment in strains[strain]:
            print("extracting speeds from ", experiment)
            data_strain = data[data['Experiment'] == experiment]
            ind = filterDataRange(data_strain, inputs)
            df_speed=pd.concat([df_speed,   data_strain.loc[ind, 'instant_speed'].reset_index(drop=True)],axis=1,ignore_index=True)
    # Give column names based on strains list of names
    df_speed.columns = [item for sublist in strains for item in sublist]
    return df_speed

def computeSpeedDist(df_speed, inputs):
    exp_WT = getInfoFromListOfDict(inputs, 'exp_WT')
    exp_aS = getInfoFromListOfDict(inputs, 'exp_aS')
    exp_As = getInfoFromListOfDict(inputs, 'exp_As')
    strains_names = ['exp_WT', 'exp_aS', 'exp_As']
    
    strains = [exp_WT,exp_aS,exp_As]
    Exp_dT = getInfoFromListOfDict(inputs, 'Exp_dT')    
    pix_size = getInfoFromListOfDict(inputs, 'pix_size')            # Pixel size in µm
    
    # -- Compute distributions        
    df_dist = pd.DataFrame()
    bins = np.linspace(0, 10, 250)
    for strain in range(len(strains)):
        fig, ax = plt.subplots(figsize=(8,8))
        for experiment in strains[strain]:
            df_dist[experiment], _, _ = plt.hist(df_speed[experiment]*pix_size/Exp_dT[experiment]*60, bins=bins,histtype=u'step', density=True)
        
        plt.ylim([1e-05, 1e1])
        plt.xlabel("Instant speed, µm/min", fontsize = 36)
        plt.gca().set_xlim(0,5)
        plt.ylabel("Normalized density", fontsize = 36)
        plt.gca().set_yscale('log')
        plt.title(strains_names[strain])
        fig.tight_layout()
        
    fig, ax = plt.subplots(figsize=(8,8))
    for strain_names in strains_names:
        Mean = df_dist[locals()[strain_names]].mean(axis=1)
        Std = df_dist[locals()[strain_names]].std(axis=1)
        ax.plot(bins[:-1], Mean, label=strain_names)
        ax.fill_between(bins[:-1], Mean - Std, Mean + Std, alpha=0.2)
    plt.ylim([1e-05, 1e1])
    plt.xlabel("Instant speed, µm/min", fontsize = 24)
    plt.gca().set_xlim(0,5)
    plt.ylabel("Normalized density", fontsize = 24)
    plt.gca().set_yscale('log')
    plt.legend()
    fig.tight_layout()
    
# ======================================================================================================================
#                                                   MAIN
# ======================================================================================================================

if __name__ == '__main__':

    scriptname = os.path.basename(sys.argv[0])
    
    # Reads the yaml document with scripts parameters
    yamlDoc = open("ScriptsInfos.yaml", 'r')
    inputs, _ = readYamlDoc(yamlDoc, scriptname)
    
    # Loads all experiments data
    data = loadAllDataframes(inputs)
    
    # Filter data range and extract speed
    df_speed = extractSpeedFromData(data, inputs)
    
    # Computes speed distributions
    computeSpeedDist(df_speed, inputs)

    

    
    
      
    
    
    