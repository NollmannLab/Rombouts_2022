#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 12:18:59 2022

@author: legall
"""



''' 
Loads dataframe files with all features for all experiments and computes the mean
of normalized 2D histograms for each condition.
Displays the difference of one confition versus another.

Usage:
- modify DATA PATH to match the dataframe files of the experiments
- provide the desired experiment names in the lists corresponding to each condition (exp_WT, exp_As, exp_aS)
- Make sure to set the parameter "distance_norm_min" to either 2 or 0 to filter
cells that did not move or to not filter them out, resp.
'''

#%% IMPORTS ###################################################################
import os
import yaml
import sys
from matplotlib import cm
from astropy.visualization import simple_norm
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np
import matplotlib.pyplot as plt

import datashader as ds

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

def normalize2max(im):
    ''' normalize to max '''
    im = im-np.min(im)
    return im/np.max(im)

def normalize2sum(im):
    ''' normalize to sum '''
    return im/np.sum(im)

def load_dataframe_exp(experiment):
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
    elif experiment == 'ex70':
        foldername = 'XXX/' 
        folderpath = '/mnt/tronador/XXX/'+foldername 
        fh = folderpath+"dataframe/all_features_WithNAN_win_size_5.pkl"
    elif experiment == 'exp71':
        foldername = 'Experiment_71_A+S-/' 
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername 
        fh = folderpath+"dataframe/all_features_WithNAN_win_size_5.pkl"    
    elif experiment == 'exp72':
        foldername = 'XX/' 
        folderpath = '/mnt/tronador/XXX/'+foldername 
        fh = folderpath+"dataframe/all_features_WithNAN_win_size_5.pkl"
    else:
        print()
    
    
    with open(fh, "rb") as fh:
      data.append(pickle.load(fh))
    
    data = pd.concat(data, axis=1)     
    print('\n *** Dataframe data for', experiment ,'successfully loaded ***\n')
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
    
def computesHist2D(inputs):
    '''
    Loads dataframe files for lists of experiments and computes 2D histograms for each experiment

    Parameters
    ----------
    inputs : inputs parameters from ScriptsInfos.yaml

    Returns
    -------
    arr : Normalized 2D Histograms for each histogram.
    weights : Weights for each experiment histogram

    '''
    exp_WT = getInfoFromListOfDict(inputs, 'exp_WT')
    exp_aS = getInfoFromListOfDict(inputs, 'exp_aS')
    exp_As = getInfoFromListOfDict(inputs, 'exp_As')
    
    strains = [exp_WT,exp_aS,exp_As]
    strains_names = ['exp_WT', 'exp_aS', 'exp_As']
    
    
    arr = {}        #Stores normalized histograms
    weights = {}    #Stores each histogram sum to do weighted average
    
    t0 = time.time()
    for strain in range(len(strains)):
        for experiment in strains[strain]:
            #%% LOADS TRACKS DATA
            print('Loading experiment '+ experiment)
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
            plot_width = getInfoFromListOfDict(inputs, 'plot_width')
            plot_height = getInfoFromListOfDict(inputs, 'plot_height')
            
            # # Default plot ranges:
            x_range = getInfoFromListOfDict(inputs, 'x_range')
            y_range = getInfoFromListOfDict(inputs, 'y_range')
    
            
            
            cvs = ds.Canvas(x_range=x_range, y_range=y_range, plot_width=plot_width, plot_height=plot_height)
            agg = cvs.points(df,X_field, Y_field)
            
            mat = agg.data.astype('float64')
            weights[experiment] = np.sum(mat)
            mat = normalize2sum(mat)
    
            arr[experiment] = mat

    print('Elapsed time: %s' % (time.time() - t0)) 
    return arr, weights
    
def Hist2D_Mean(arr, weights, inputs):  
    exp_WT = getInfoFromListOfDict(inputs, 'exp_WT')
    exp_aS = getInfoFromListOfDict(inputs, 'exp_aS')
    exp_As = getInfoFromListOfDict(inputs, 'exp_As')
    
    strains = [exp_WT,exp_aS,exp_As]
    strains_names = ['exp_WT', 'exp_aS', 'exp_As']
    #%%  ------ DATASHADER AGGREGATOR PARAMETERS
    plot_width = getInfoFromListOfDict(inputs, 'plot_width')
    plot_height = getInfoFromListOfDict(inputs, 'plot_height')
    
    # # Default plot ranges:
    x_range = getInfoFromListOfDict(inputs, 'x_range')
    y_range = getInfoFromListOfDict(inputs, 'y_range')
    #%% Compute mean 2D histograms
    weighted_mean = True
    arr_mean = np.zeros((len(strains_names), plot_width,plot_height)) 
    tot_weights = {}
    if weighted_mean:
        for count, strain in enumerate(strains_names):
            tot_weights[strain] = 0
            for experiment in strains[count]:
                tot_weights[strain] = tot_weights[strain]+weights[experiment]
    
    
    for count, strain in enumerate(strains_names):
        for experiment in strains[count]:
            if weighted_mean:
                arr_mean[count] = arr_mean[count] + arr[experiment]*weights[experiment]/tot_weights[strain]
            else:
                arr_mean[count] = arr_mean[count] + arr[experiment]/len(strains[count])
    
    return arr_mean

def disp_Hist2D_Mean(arr_mean, inputs, i_strain=0):
    '''
    

    Parameters
    ----------
    arr_mean : TYPE
        DESCRIPTION.
    inputs : TYPE
        DESCRIPTION.
    i_strain : element number in list ['exp_WT', 'exp_aS', 'exp_As']  

    Returns
    -------
    None.

    '''
    #%%  Display single histograms
    strains_names = ['exp_WT', 'exp_aS', 'exp_As']  
    # Make custom colormap
    cmp = cm.get_cmap('Greens', 512)
    newcolors = cmp(np.linspace(0,1, 512))
    newcolors[0, :] = np.array([1, 1, 1, 1])
    newcmp = ListedColormap(newcolors)
            
    # Histogram to plot
    mat = arr_mean[i_strain]
            
    # Default plot ranges:
    x_range = getInfoFromListOfDict(inputs, 'x_range')
    y_range = getInfoFromListOfDict(inputs, 'y_range')
        
    fig, ax = plt.subplots()
    
    norm = simple_norm(mat, 'log')
    im = ax.imshow(mat, origin='lower', 
                    norm=colors.LogNorm(vmin=mat[mat>0].min(), vmax=mat.max()), 
                    cmap=newcmp, extent=[x_range[0], x_range[1],y_range[0], y_range[1]],
                    interpolation='none')
    plt.ylabel("Local density (log)")
    plt.xlabel("Cluster size (log, au)")
    plt.title(strains_names[i_strain], y=1.24)
    
    
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
    
    
def HistComp(arr_mean, inputs, i_ref=0, i_diff=1, op = '-'):
    '''
    

    Parameters
    ----------
    arr_mean : TYPE
        DESCRIPTION.
    inputs : TYPE
        DESCRIPTION.
    i_ref : reference histogram. Default is WT, i_ref=0
    i_diff : hist to make difference with. The default is 1 (exp_aS).
    op : operator to use (either - or /)
    Returns
    -------
    None.

    '''
    strains_names = ['exp_WT', 'exp_aS', 'exp_As']    
    # =============================================================================
    #                           Histogram differences
    # =============================================================================
    #%%  Custom diverging colormaps
    
        
    colors_list = ["aqua","blue","darkblue", "black", "black", "darkred","red","yellow"]
    nodes = [0.0, 0.3, 0.4, 0.4999, 0.5001, 0.6, 0.7, 1.0]
    cmap2 = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colors_list)))
    
    colors_list = ["blue","darkblue", "black", "black", "darkred","red"]
    nodes = [0., 0.4, 0.4999, 0.5001, 0.6, 1]
    cmap = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colors_list)))
    
    #%% Plot histograms ratios with custom colormpas
    
    print('Displaying ', strains_names[i_diff], op, '', strains_names[i_ref], ' Histograms')
    
    if op == '/':
        mat = arr_mean[i_diff]/arr_mean[i_ref]
        mat[mat==np.inf]= np.nanmax(mat[mat != np.inf])
        # If colormap centered on 1
        linthresh=10**-5
        norm=colors.SymLogNorm(linthresh=linthresh, linscale=linthresh, vmin=np.nanmin(mat[mat>0]), vmax=np.nanmax(mat), base=10)
        cmap = cmap
    else:
        mat = arr_mean[i_diff]-arr_mean[i_ref]
        mat[mat==0]=np.nan
        mat[mat==np.inf]=np.nan
        # If colormap centered on 0
        linthresh=0.000001#0.000001
        norm=colors.SymLogNorm(linthresh=linthresh, linscale=linthresh*1000, vmin=np.nanmin(mat), vmax=-np.nanmin(mat), base=10)
        cmap = cmap2
        
    # Default plot ranges:
    x_range = getInfoFromListOfDict(inputs, 'x_range')
    y_range = getInfoFromListOfDict(inputs, 'y_range')
            
    fig, ax = plt.subplots()
    im = ax.imshow(mat, origin='lower', norm=norm, cmap=cmap, 
                        extent=[x_range[0], x_range[1],y_range[0], y_range[1]], 
                        interpolation='none')
    plt.ylabel("Local density (log)")
    plt.xlabel("Cluster size (log, au)")
    plt.title('{0} {1} {2}'.format(strains_names[i_diff], op, strains_names[i_ref]), y=1.24)
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
    
        
        # Computes MAIN
        arr, weights = computesHist2D(inputs)
        
        # Mean histograms for ['exp_WT', 'exp_aS', 'exp_As']   
        arr_mean = Hist2D_Mean(arr, weights, inputs) 
        
        # Display one single histogram
        disp_Hist2D_Mean(arr_mean, inputs, i_strain=0)
        
        # Computes difference of2 histograms
        HistComp(arr_mean, inputs, i_ref=0, i_diff=1, op = '-')
    
    
    
        