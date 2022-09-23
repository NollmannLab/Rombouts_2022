#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 16:19:20 2022

@author: legall
"""

''' 
Loads dataframe files with all features for all experiments and computes the 
transitions between classes.

Usage:
- modify DATA PATH to match the dataframe files of the experiments
- provide the desired experiment names in the lists corresponding to each 
condition (exp_WT, exp_As, exp_aS)
- Modify the variable "strain" in Main to pick the strain of interest to display.
'''

#%% IMPORTS ###################################################################
import os
import yaml
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle5 as pickle
import time

import datashader as ds

from matplotlib import cm
from matplotlib.colors import ListedColormap
from astropy.visualization import simple_norm
from matplotlib.collections import LineCollection
import  matplotlib as mpl


#%%

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
# mpl.rcParams["axes.set_box_aspect"] = 1
#### PLOTS
mpl.rcParams['lines.linewidth'] = 4

#### TICKS
# mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["xtick.top"] = False
mpl.rcParams["xtick.major.width"] = 2
mpl.rcParams["xtick.major.size"] = 10
# mpl.rcParams["xtick.minor.visible"] = True
mpl.rcParams["xtick.minor.width"] = 2
mpl.rcParams["xtick.minor.size"] = 5

# mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["ytick.right"] = False
mpl.rcParams["ytick.major.width"] = 2
mpl.rcParams["ytick.major.size"] = 10
# mpl.rcParams["ytick.minor.visible"] = True
mpl.rcParams["ytick.minor.width"] = 2
mpl.rcParams["ytick.minor.size"] = 5

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
    return data

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
            data_temp, _ = load_dataframe_exp(experiment)
            data_temp=data_temp[['Total_Col', 'Total_Row','TimeStamp', 'TrackID', 'Total_Density', 'Nb_cells']]
            data_temp['Experiment'] = experiment
            data_temp['TrackID'] = data_temp['TrackID'] + TrackID_max
            TrackID_max = data_temp['TrackID'].max()
            data.append(data_temp)  
    
    data = pd.concat(data, axis=0, ignore_index=True) 
    return data

def filterDataRange(data, inputs):
    # -- Filter data range
    classes = ['herd', 'loners','scouts']
    class_sel = [0,1,2]
    data['Class'] = np.nan
    for i in range(len(class_sel)):
        if classes[class_sel[i]] == 'scouts' :
            density_min = 4.5
            density_max = 10**6
            Nb_cells_min = 0
            Nb_cells_max = 20
        elif classes[class_sel[i]] == 'loners' :
            density_min = 0
            density_max = 4.5
            Nb_cells_min = 0
            Nb_cells_max = 2
        elif classes[class_sel[i]] == 'herd' :
            density_min = 0
            density_max = 4.5
            Nb_cells_min = 0
            Nb_cells_max = 10**6
        elif classes[class_sel[i]] == 'rafts' :
            density_min = 0
            density_max = 10**6
            Nb_cells_min = 2
            Nb_cells_max = 600
        elif classes[class_sel[i]] == 'swarms' :
            density_min = 0
            density_max = 10**6
            Nb_cells_min = 600
            Nb_cells_max = 10**6
        elif classes[class_sel[i]] == 'all' :
            density_min = 0
            density_max = 10**6
            Nb_cells_min = 0
            Nb_cells_max = 10**6

        distance_norm_min = 0
        distance_norm_max = 10**6
        Mask_width = 400
        t1 = 0
        t2 = 10**6
        ind = data.index[ (data['Total_Density'] >= density_min) & (data['Total_Density'] <= density_max)
                          & (data['Nb_cells'] >= Nb_cells_min) & (data['Nb_cells'] <= Nb_cells_max)
                          & (data['distance_norm'] >= distance_norm_min) & (data['distance_norm'] <= distance_norm_max)
                          & (data['Total_Col'] >= Mask_width) & (data['Total_Col'] <= 5880- Mask_width)
                          & (data['Total_Row'] >= Mask_width) & (data['Total_Row'] <= 5880- Mask_width)
                          & (data['TimeStamp'] >= t1) & (data['TimeStamp'] <= t2)].tolist()
        
        ind_endTracks = data.index[ (data['Total_Density'] == -1) & (data['Nb_cells'] == -1)].tolist()
        ind = np.sort(ind+ind_endTracks)
        data.loc[ind, 'Class'] = class_sel[i]
    return data


def mostfrequent(array):
    if True in np.isreal(array):
        u, c = np.unique(array[~np.isnan(array) & (array>-1)], return_counts=True)
        if c.size==0:
            output = -2
        else:
            output = u[c.argmax()]
    else:
        output = -2
    return output

def computeClasses(data, inputs):
    data['Class_row'] = np.nan
    data['Class_col'] = np.nan
    
    # Do a rolling window to smooth class transitions that are too transient
    print('\n *** Applying rolling filter on classes, please wait... ***')
    t = time.time()
    temp = data.groupby('TrackID')['Class'].rolling(10, center=True,min_periods=1).apply(mostfrequent).reset_index()
    print('Elapsed: %s' % (time.time() - t))   
    
    
    data.loc[(temp['Class'] >= 0) & (temp['Class']< 0.5), 'Class_row'] = 0
    data.loc[(temp['Class'] >= 0) & (temp['Class']< 0.5), 'Class_col'] = 0
    
    data.loc[(temp['Class'] >= 0.5) & (temp['Class']< 2), 'Class_row'] = 2
    data.loc[(temp['Class'] >= 0.5) & (temp['Class']< 2), 'Class_col'] = 0
    
    data.loc[(temp['Class'] >= 2) & (temp['Class']< 3.5), 'Class_row'] = 1
    data.loc[(temp['Class'] >= 2) & (temp['Class']< 3.5), 'Class_col'] = np.sqrt(3)

    #%% -- Sort dataframe
    # Sort dataframe rows by TrackID then Time
    data = data.sort_values(["TrackID", "TimeStamp"], ascending = (True, True))  
    
    # Find ends of tracks (change of track ID) and add NAN to be recognized as individual tracks with Datashader
    data["diff"] = data["TrackID"].diff(periods=-1)
    RowsToCopy= data.loc[(data["diff"] != 0) , :].copy()
    
     # Set values to np.nan in fake row except in ['TimeStamp','TrackID']
    cols = [col for col in RowsToCopy.columns if col not in ['TimeStamp','TrackID']]
    RowsToCopy[cols] = np.nan
    RowsToCopy['TimeStamp'] = RowsToCopy['TimeStamp']+1
    RowsToCopy['Class_row'] = np.nan
    RowsToCopy['Class_col'] = np.nan
    RowsToCopy['Total_Density'] = -1
    RowsToCopy['Nb_cells'] = -1
    RowsToCopy['Class'] = -1
    data = pd.concat([data, RowsToCopy])

    # Sort again dataframe rows by TrackID then Time
    data = data.sort_values(["TrackID", "TimeStamp"], ascending = (True, True)) 
    
    return data, temp

def buildMap(data):
    #%% ------ DATASHADER
    data["Diff"] = data["Class"].diff(periods=-1)
    data["Diff"] = np.sign(data["Diff"])
    
    # Get fields
    X_field = 'Class_row'
    Y_field = 'Class_col'
    C_field = 'Diff'
    
    df = data[[X_field, Y_field, C_field,"Class"]]
    df = df[df['Class'].notna()]
    noise = np.random.normal(0, .1, df[X_field].shape)
    df[X_field] =  (df[X_field]+5)*3 + noise*2
    noise = np.random.normal(0, .1, df[X_field].shape)
    df[Y_field] =   (df[Y_field]+5)*3 + noise*2
    
    # # Default plot ranges:
    x_range = (df[X_field].min() , df[X_field].max())
    y_range = (df[Y_field].min() , df[Y_field].max())
    
    plot_width = 5000#round((x_range[1]+1)/5)
    plot_height = 5000#round((y_range[1]+1)/5)
    
    cvs = ds.Canvas(x_range=x_range, y_range=y_range, plot_width=round(plot_width), plot_height=round(plot_height))
    
    agg = cvs.line(df,X_field, Y_field,  ds.count())
    
    return agg, x_range, y_range


def dispMap(agg, x_range, y_range):
    # ------ Edit colormap
    cmp = cm.get_cmap('jet', 10**6)
    newcolors = cmp(np.linspace(0,1, 10**6))
    newcolors[0, :] = np.array([1, 1, 1, 1])
    newcmp = ListedColormap(newcolors)  
    
    fig, ax = plt.subplots()
    norm = simple_norm(agg.values, 'linear')
    plt.imshow(agg.values, cmap = newcmp, extent=[x_range[0], x_range[1],y_range[0], y_range[1]], norm = norm, origin='lower', interpolation="none") 
    ax.axis('square')
    ax.set_yticks([], minor=False)
    ax.set_xticks([], minor=False)

    ax.set_ylim(13.5, 21.5)
    ax.set_xlim(13.5, 22.5)
    
def transitionHist(data):
    fig, ax = plt.subplots()
    data[data['Class']>-1]['Class'].hist()
    
    #%% SHOW FREQUENCY OF TRANSITIONING BETWEN CLASSES
    data["Diff"] = data["Class"].diff(periods=-1)
    data['Diff'][data['Diff'] > 0] = 1
    data['Diff'][data['Diff'] < 0] = 1
    data['Diff'][data['Diff'].isnull()] = 0 # Class zero has Diff = Nan when no changes
    
    group = data.groupby(['TrackID'])['Diff'].value_counts().unstack()
    group[group.isnull()] = 0
    
    counts = group[1]
    
    _, ax = plt.subplots()
    bins = np.linspace(0, 50, 50)
    plt.hist(counts[np.isfinite(counts)].values, bins = bins)
    ax.set_yscale('log')
    x0,x1 = ax.get_xlim()
    y0,y1 = ax.get_ylim()
    ax.set_aspect(abs(x1-x0)/abs(y1-y0)*1e5)
    ax.set_xlabel('Number of transitions per track')
    ax.set_ylabel('Counts')
    # ax.set_title(strains_names[strain], y=1.04)
    
def dispMapSingleTrack(data):
    #%% PLOT SINGLE TRACKS COLORCODED WITH TIME
    
    fig, ax = plt.subplots()
    i=0
    #%%
    i=i+1
    ind_toplot = data.index[ (data['TrackID'] == i) ].tolist()
    while (len(ind_toplot)<100 or len(data.loc[ind_toplot, 'Class'].unique())<3):
        i=i+1
        ind_toplot = data.index[ (data['TrackID'] == i) ].tolist()
    
    x = data.loc[ind_toplot, 'Class_row'] 
    x = (x+5)*3 + np.random.normal(0, .1, x.shape)
    y = data.loc[ind_toplot, 'Class_col']
    y = (y+5)*3 + np.random.normal(0, .1, y.shape)
    c = data.loc[ind_toplot, 'TimeStamp']
    
    
    
    points = np.array([x[~np.isnan(x)], y[~np.isnan(x)]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    if 'lc' in locals():
        try:
            cb.remove()
        except ValueError:
            print("Oops!  cb could not be removed...")
        lc.remove()
        del lc

    lc = LineCollection(segments, cmap="jet", linewidth=3)
    # set color to date values
    lc.set_clim(0, 700)
    lc.set_array(c[~np.isnan(x)])
    ax.add_collection(lc)
    plt.title('TrackID = ' + str(i))
    
    plt.axis('square')
    plt.ylim(13.5, 21.5)
    plt.xlim(13.5, 22.5)
    cb = fig.colorbar(lc)
    print(data[(data['TrackID'] == i)]['Class'].unique())
    
def computeFluxes(data, temp):
    df_flux = data[['Class']]
    df_flux['new_Class'] = np.nan
    df_flux.loc[(temp['Class'] >= 0) & (temp['Class']< 0.5), 'new_Class'] = 0
    df_flux.loc[(temp['Class'] >= 0.5) & (temp['Class']< 2), 'new_Class'] = 1
    df_flux.loc[(temp['Class'] >= 2) & (temp['Class']< 3.5), 'new_Class'] = 3
    df_flux.loc[temp['Class'] >= 3.5, 'new_Class'] = 7
    df_flux["Flux"] = df_flux["new_Class"].diff(periods=-1)
    
    fig, ax = plt.subplots()
    ax = df_flux[df_flux["Flux"] != 0]["Flux"].hist(bins = 50)
    rects = ax.patches

    # Make some labels.
    labels = [f"label{i}" for i in range(len(rects))]
    
    for rect, label in zip(rects, labels):
        height = rect.get_height()
        if height> 0:
            ax.text(
                rect.get_x() + rect.get_width() / 2, height + 5, height, ha="center", va="bottom"
            )
    
    # #
    # fig, axs = plt.subplots(1,4)
    # i=0
    # plt.subplots_adjust(wspace = .01)
    # ind_toplot = [1277,1687,5238,24196] #exp8
    
    # for i in range(len(ind_toplot)):
    #     ind = data.index[ (data['TrackID'] == ind_toplot[i]) ].tolist()
    
    #     x = data.loc[ind, 'Class_row'] 
    #     x = (x+5)*3 + np.random.normal(0, .1, x.shape)
    #     y = data.loc[ind, 'Class_col']
    #     y = (y+5)*3 + np.random.normal(0, .1, y.shape)
    #     c = data.loc[ind, 'TimeStamp']
        
        
    #     points = np.array([x[~np.isnan(x)], y[~np.isnan(x)]]).T.reshape(-1, 1, 2)
    #     segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
      
    #     lc = LineCollection(segments, cmap="jet", linewidth=3)
        
    #     # set color to date values
        
    #     lc.set_array(c[~np.isnan(x)])
    #     axs[i].add_collection(lc)
    #     lc.set_clim(0, 700)
    #     axs[i].axis('square')
    #     axs[i].set_yticks([16.5], minor=False)
    #     axs[i].set_xticks([16.5], minor=False)
    #     axs[i].yaxis.grid(True, which='major')
    #     axs[i].xaxis.grid(True, which='major')
    #     axs[i].set_ylim(13.5, 21.5)
    #     axs[i].set_xlim(13.5, 22.5)
    #     for tick in axs[i].xaxis.get_major_ticks():
    #         tick.tick1line.set_visible(False)
    #         tick.tick2line.set_visible(False)
    #         tick.label1.set_visible(False)
    #         tick.label2.set_visible(False)
    #     for tick in axs[i].yaxis.get_major_ticks():
    #         tick.tick1line.set_visible(False)
    #         tick.tick2line.set_visible(False)
    #         tick.label1.set_visible(False)
    #         tick.label2.set_visible(False)
    #     axs[i].set_xticklabels([])
    #     axs[i].set_yticklabels([])
#%% ====================================================================================================================
#                                                   MAIN
# ======================================================================================================================

if __name__ == '__main__':
    scriptname = os.path.basename(sys.argv[0])
    
    # Reads the yaml document with scripts parameters
    yamlDoc = open("ScriptsInfos.yaml", 'r')
    inputs, _ = readYamlDoc(yamlDoc, scriptname)
    
    # Loads all experiments data
    data = loadAllDataframes(inputs)
    print('Data loaded')
    data = formatDataframe(data)
    print('Data formated')
    data = filterDataRange(data, inputs)
    print('Data filtered')
    data, temp = computeClasses(data, inputs)
    
    # Get lists of experiments grouped by strains
    exp_WT = getInfoFromListOfDict(inputs, 'exp_WT')
    exp_aS = getInfoFromListOfDict(inputs, 'exp_aS')
    exp_As = getInfoFromListOfDict(inputs, 'exp_As')
    
    # Pick the group of experimetns to display
    strain = exp_aS     #  <-- Modify this variable to pick the strain of interest.
    data_to_plot = data[data['Experiment'].isin(strain)]
    
    print('Classes computed')
    agg, x_range, y_range = buildMap(data_to_plot)
    print('Map computed')
    dispMap(agg, x_range, y_range)
    print('Map displayed')
    transitionHist(data_to_plot)
    print('Transitions histogram displayed')
    computeFluxes(data_to_plot, temp)
    print('Fluxes computed')
    
   
    
    
    
    