#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 15:22:51 2021

@author: legall
"""


''' 
Loads dataframe file with all features for a single experiment and computes
2D histograms of the Voronoi area (wrongly termed density here) vs cluster size.

Usage:
- modify DATA PATH to match the dataframe files of the experiments in function "load_dataframe_exp"
- provide the desized experiment name in parameter "experiment" of ScriptsInfos.yaml
'''

#%% IMPORTS ###################################################################
import os
import yaml
import sys
from matplotlib import cm
from matplotlib.colors import ListedColormap
from astropy.visualization import simple_norm
import matplotlib.colors as colors
import numpy as np
import matplotlib.pyplot as plt

import datashader as ds

import pandas as pd

import pickle5 as pickle

import  matplotlib as mpl


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
   
def load_dataframe_exp(experiment):
    '''
    Loads the dataframe analysis data (single file).
    
    :param experiment: string with the name of the experment to load.
    :return: data: dataframe analyusis data
    :return: folderpath: location of the dataframe data.
    '''
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
    print(' *** Tracks data successfully loaded ***\n')
    return data, folderpath


def squarify(fig):
    '''
    Makes the input figure aspect ratio square.

    Parameters
    ----------
    fig : Figure handle.

    Returns
    -------
    None.

    '''
    w, h = fig.get_size_inches()
    if w > h:
        t = fig.subplotpars.top
        b = fig.subplotpars.bottom
        axs = h*(t-b)
        l = (1.-axs/w)/2
        fig.subplots_adjust(left=l, right=1-l)
    else:
        t = fig.subplotpars.right
        b = fig.subplotpars.left
        axs = w*(t-b)
        l = (1.-axs/h)/2
        fig.subplots_adjust(bottom=l, top=1-l)
        
        
def Hist2D_twinaxis(inputs):
    '''
    Computes 2D histograms using Datashader.

    Parameters
    ----------
    inputs : inputs parameters from ScriptsInfos.yaml

    Returns
    -------
    None.

    '''
    #%% LOADS TRACKS DATA
    experiment = getInfoFromListOfDict(inputs, 'experiment')
    data, folderpath = load_dataframe_exp(experiment)
    
    print(data.columns.tolist())

    #%% -- Compute max travelled distances
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
    #%% -- Modify fields
    data['clust_area_x_mean'] = data.groupby(['TrackID'])['clust_area_x'].transform('mean') 
    
    fields_to_log = ['clust_area_x','clust_area_y','clust_area_x_mean','Total_Density_mean']
    data[fields_to_log] = np.log10(data[fields_to_log])
    fields_to_10 = ['instant_speed','instant_speed_from', 'instant_speed_to','instant_speed_mean', 'instant_speed_min', 'instant_speed_max', 'instant_speed_std']
    data[fields_to_10] = 10**(data[fields_to_10])
    
    #%% -- Filter data range
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
   
    #%% ------ FORMAT FOR DATASHADER DISPLAY
    # Get fields
    X_field = getInfoFromListOfDict(inputs, 'X_field')
    Y_field = getInfoFromListOfDict(inputs, 'Y_field')
    C_field = 'instant_speed'
    
    df = data.loc[ind, [X_field, Y_field, C_field]]


    #%%  ------ BUILD DATASHADER AGGREGATOR
    # Make custom colormap
    cmp = cm.get_cmap('Greens', 512)
    newcolors = cmp(np.linspace(0,1, 512))
    newcolors[0, :] = np.array([1, 1, 1, 1])
    newcmp = ListedColormap(newcolors)
        
        
    
    plot_width = getInfoFromListOfDict(inputs, 'plot_width')
    plot_height = getInfoFromListOfDict(inputs, 'plot_height')
    
    # # Default plot ranges:
    x_range = getInfoFromListOfDict(inputs, 'x_range')
    y_range = getInfoFromListOfDict(inputs, 'y_range')
    
    
    cvs = ds.Canvas(x_range=x_range, y_range=y_range, plot_width=plot_width, plot_height=plot_height)
    
    agg = cvs.points(df,X_field, Y_field)
    
    mat = agg.data.astype('float64')
    
    #%%  Display results
    
    fig, ax = plt.subplots()
    
    norm = simple_norm(mat, 'log')
    im = ax.imshow(mat, origin='lower', 
                   norm=colors.LogNorm(vmin=mat[mat>0].min(), vmax=mat.max()), 
                   cmap=newcmp, extent=[x_range[0], x_range[1],y_range[0], y_range[1]],
                   interpolation='none')
    plt.ylabel("Local density (log)")
    plt.xlabel("Cluster size (log, au)")
    plt.title(experiment, y=1.24)
    
    
    ax.set(xlim = [x_range[0], x_range[1]], ylim = [y_range[0], y_range[1]])
    
    
    fig.colorbar(im)
    plt.show()
    
    ax2 = ax.twiny()
    xmin = x_range[0]-np.log10(550)
    xmax = x_range[1]-np.log10(550)
    ax2.set(xlim = [xmin, xmax])
    ax2.set_xlabel('Number of cells per cluster (log)')
    
    squarify(fig)
    fig.canvas.mpl_connect("resize_event", lambda evt: squarify(fig))
    
# ======================================================================================================================
#                                                   MAIN
# ======================================================================================================================

if __name__ == '__main__':

    
    scriptname = os.path.basename(sys.argv[0])
    
    # Reads the yaml document with scripts parameters
    yamlDoc = open("ScriptsInfos.yaml", 'r')
    inputs, _ = readYamlDoc(yamlDoc, scriptname)

    # Computes 2D histogram   
    Hist2D_twinaxis(inputs)
    
