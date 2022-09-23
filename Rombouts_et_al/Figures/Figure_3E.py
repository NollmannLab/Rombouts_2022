#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 10:16:36 2022

@author: legall
"""



''' 
Loads datashader maps generated with Figure3_generate_tracks_datashader_maps_v20220207.py
and computes similarity measures.

Usage:
- modify DATA PATH to match the datashader maps location

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
from skimage.metrics import structural_similarity
from skimage.measure import label, regionprops_table

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

def load_datashader_maps(experiment, plot_width):
    '''
    

    Parameters
    ----------
    experiment : TYPE
        DESCRIPTION.
    plot_width : TYPE
        DESCRIPTION.

    Returns
    -------
    data : TYPE
        DESCRIPTION.
    folderpath : TYPE
        DESCRIPTION.

    '''
    # LOADS TRACKS DATA
    
    data = []
    if experiment == 'expX':
        foldername = 'Myxo_Sara_Data_with_EC_expX/'
        folderpath = '/mnt/grey/DATA/users/Antoine/'+foldername 
    elif experiment == 'expX_bis':
        foldername = 'Myxo_Sara_Data_with_EC_expX_bis/'
        folderpath = '/mnt/grey/DATA/users/Antoine/'+foldername 
    elif experiment == 'exp6':
        foldername = 'Experiment_6_sara/' 
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername 
    elif experiment == 'exp8':
        foldername = 'Myxo_Sara_Data/' 
        folderpath = '/mnt/grey/DATA/users/Antoine/'+foldername 
    elif experiment == 'exp9':
        foldername = 'Myxo_Sara_Data_with_EC_exp9/'
        folderpath = '/mnt/grey/DATA/users/Antoine/'+foldername 
    elif experiment == 'exp10_DK':
        foldername = 'Myxo_Sara_Data_Exp10_DK/'
        folderpath = '/mnt/grey/DATA/users/Antoine/'+foldername 
    elif experiment == 'exp11':
        foldername = 'Myxo_Sara_Data_Exp11/'
        folderpath = '/mnt/grey/DATA/users/Antoine/'+foldername 
    elif experiment == 'exp12':
        foldername = 'Myxo_Sara_Data_Exp12/'
        folderpath = '/mnt/grey/DATA/users/Antoine/'+foldername 
    elif experiment == 'exp50':
        foldername = 'Experiment_50_Sara_WT/'
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername 
    elif experiment == 'exp51':
        foldername = 'Experiment_51_Sara_A+S-/'
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername 
    elif experiment == 'exp53':
        foldername = 'Experiment_53_Sara_A+S-/'  
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername 
    elif experiment == 'exp54':
        foldername = 'Experiment_54_Sara_A-S+/'  
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername 
    elif experiment == 'exp55':
        foldername = 'Experiment_55_Sara_A-S+/'   
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername 
    elif experiment == 'exp66':
        foldername = 'Experiment_66_Anna_WT_SideOfECColony/' 
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername 
    elif experiment == 'exp70':
        foldername = 'Experiment_70_Sara_WT_cytosolicGFP/' 
        folderpath = '/mnt/tronador/AnalyzedData_2021/'+foldername 
    elif experiment == 'exp71':
        foldername = 'Experiment_71_A+S-/' 
        folderpath = '/mnt/grey/DATA/AnalyzedData_2021/'+foldername 
    elif experiment == 'exp72':
        foldername = 'Experiment_72_Sara_AminusSplus/' 
        folderpath = '/mnt/tronador/AnalyzedData_2021/'+foldername 
        
    else:
        print('Not yet implemented.')
    
    fh = folderpath+'datashader/'+ experiment+'_tracks_classes' +'_'+str(plot_width)+ '_v3.pkl'
    with open(fh, "rb") as fh:
      data.append(pickle.load(fh))
    
    print(' *** Tracks data successfully loaded ***\n')
    return data, folderpath

def loadDatashaders(inputs):
    exp_WT = getInfoFromListOfDict(inputs, 'exp_WT')
    exp_aS = getInfoFromListOfDict(inputs, 'exp_aS')
    exp_As = getInfoFromListOfDict(inputs, 'exp_As')
    
    strains = [exp_WT,exp_aS,exp_As]
        
    plot_width = 5880

    experiments = [item for sublist in strains for item in sublist]
    data = {}
    
    for experiment in experiments:
        data[experiment] , folderpath = load_datashader_maps(experiment, plot_width)
    print(' ***     All tracks data successfully loaded     ***\n')

    print(' *** Datashader data successfully loaded ***\n')
    return data


def compute_similarity_map(track_dict, classes_to_compare=['scouts','swarm']):
    '''
    Takes a data dictionnary of classes maps and compute the similarity map
    track_dict = data[experiment][0] # Dictionnary with tracks arrays per class
    
    
    Parameters
    ----------
    track_dict : data[experiment][0] # Dictionnary with tracks arrays per class
    
    classes_to_compare : List of 2 classes to compare 
        DESCRIPTION. The default is ['scouts','herd'].

    Returns
    -------
    result : TYPE
        DESCRIPTION.
    imA : TYPE
        DESCRIPTION.
    imB : TYPE
        DESCRIPTION.

    '''
    
    classes = list(track_dict.keys())         # List of population classes
    class_sel = [i for i in range(len(classes)) if classes[i] in classes_to_compare]
    
    
    # Get tracks images and normalize them
    im1 = track_dict[classes[class_sel[0]]]
    im1 = (im1 - im1.min()) / (im1.max() - im1.min())
    im2 = track_dict[classes[class_sel[1]]]
    im2 = (im2 - im2.min()) / (im2.max() - im2.min())
    
    # threshold
    im1 = im1>0
    im2 = im2>0             
    
    # Remove borders
    pad = round(im1.shape[0]*0.04)
    mask = np.ones((im1.shape[0]-2*pad,im1.shape[1]-2*pad))
    mask = np.pad(array=mask, pad_width=pad, mode='constant', constant_values=0)
    im1 = im1*mask
    im2 = im2*mask
    
    im1 = (im1*512).astype('uint16')
    im2 = (im2*512).astype('uint16')
    
    try:
        data_range=im2.max()-im2[im2>0].min()
    except ValueError:  #raised if `y` is empty.
        data_range=im2.max()-0

    (score, diff) = structural_similarity(im1, im2, win_size=3, full=True, data_range=data_range)
    # (score, diff) = structural_similarity(im1, im2, win_size=3, full=True, data_range=im2.max()-im2[im2>0].min())
    
    print("Image similarity", score)
    
    result = abs(diff*(im1>im1.min()))
    result = (result * 2**16).astype("uint16")
 
    return result, im1, im2

def compute_similarity_map_all(data, inputs):
    exp_WT = getInfoFromListOfDict(inputs, 'exp_WT')
    exp_aS = getInfoFromListOfDict(inputs, 'exp_aS')
    exp_As = getInfoFromListOfDict(inputs, 'exp_As')
    strains = [exp_WT,exp_aS,exp_As]
    experiments = [item for sublist in strains for item in sublist]
    result = {}
    for experiment in experiments:
        result[experiment],_,_ = compute_similarity_map(data[experiment][0], classes_to_compare=['scouts','swarms'])
        print('--> ', experiment, ' similarity map computed')
    return result
    
    
def quantifySimilarity_to_df(result, inputs):
    exp_WT = getInfoFromListOfDict(inputs, 'exp_WT')
    exp_aS = getInfoFromListOfDict(inputs, 'exp_aS')
    exp_As = getInfoFromListOfDict(inputs, 'exp_As')
    strains = [exp_WT,exp_aS,exp_As]
    experiments = [item for sublist in strains for item in sublist]
    
    # ------ SIMILARITY QUANTIFICATION
    field = 'area' #'label','centroid', 'area','perimeter', 'major_axis_length', 'eccentricity','solidity','orientation', 'mean_intensity'
    style = 'bar0'
    scaling = 'log'
    norm = False
    nb_bins = 50
    
    xdata = []
    ydata = []
    for experiment in experiments:
        label_result = label(result[experiment]>0, connectivity=result[experiment].ndim, background=0)
        props_result = regionprops_table(label_result,result[experiment], properties=[field], cache = True)
        df_props_result = pd.DataFrame(props_result)
        # log-scaled bins
        if experiment==experiments[0]:
            if scaling == 'log':
                bins = np.logspace(0, np.log10(df_props_result[field].max()), nb_bins)
            else:
                bins = np.linspace(df_props_result[field].min(), df_props_result[field].max(), nb_bins)
            widths2 = (bins[1:] - bins[:-1])
            bins2 = bins
        else:
            bins = bins2
        
        # # Calculate histogram
        hist = np.histogram(df_props_result[field], bins=bins)
        # normalize by bin width
        if norm:
            hist_norm = hist[0] / widths2
            hist_norm = hist_norm / np.sum(hist_norm)
        else:
            hist_norm = hist[0]
        
        xdata_temp = bins[:-1]
        ydata_temp = hist_norm
        
        xdata.append(xdata_temp)
        ydata.append(ydata_temp)
        
    df = pd.DataFrame({'X': xdata, 'Y': ydata}, index = experiments)
    return df, widths2


def disp_SimQuantif_results(df, widths2, inputs, field = 'area', style = 'bar0', scaling = 'log', norm = False, nb_bins = 50):
    exp_WT = getInfoFromListOfDict(inputs, 'exp_WT')
    exp_aS = getInfoFromListOfDict(inputs, 'exp_aS')
    exp_As = getInfoFromListOfDict(inputs, 'exp_As')
    strains_names = ['exp_WT', 'exp_aS', 'exp_As']
    strains = [exp_WT,exp_aS,exp_As]
    # ------ SHOW SIMILARITY QUANTIFICATION RESULTS
    ax2 = plt.subplot(1,1,1)
    for strain in range(len(strains)):
                
        bins = df['X'].mean()
        Mean = df.loc[strains[strain],'Y'].mean()
        Std = np.stack(df.loc[strains[strain],'Y']).std(axis=0)
        
        bins = bins[np.nonzero(Mean)]
        Std = Std[np.nonzero(Mean)]
        Mean = Mean[np.nonzero(Mean)]
        
        
        if style == 'bar':
            ax2.bar(bins, Mean, widths2, linewidth=0,alpha=0.5)
        else:
            ax2.plot(bins, Mean, label=strains_names[strain])
            ax2.fill_between(bins, Mean - Std, Mean + Std, alpha=0.2)
            
        ax2.set_yscale('linear')
        ax2.set_xscale(scaling)
        ax2.set_title(field, fontsize = 20)
        ax2.legend()
        ax2.set_xlabel('Correlated track size', fontsize=20)
        ax2.set_ylabel('Counts', fontsize=20)
        ax2.grid()
        if scaling == 'log':
            x0,x1 = plt.gca().get_xlim()
            y0,y1 = plt.gca().get_ylim()
            y0 = 0
            x0 = 10**0
            plt.gca().set_aspect(abs(np.log(x1)-np.log(x0))/abs(y1-y0)*0.5)
        else:
            x0,x1 = plt.gca().get_xlim()
            y0,y1 = plt.gca().get_ylim()
            plt.gca().set_aspect(abs(x1-x0)/abs(y1-y0))
            
def disp_SimQuantif_SingleResult(df, widths2, inputs, istrain = 1, field = 'area', style = 'bar0', scaling = 'log', norm = False, nb_bins = 50):
    exp_WT = getInfoFromListOfDict(inputs, 'exp_WT')
    exp_aS = getInfoFromListOfDict(inputs, 'exp_aS')
    exp_As = getInfoFromListOfDict(inputs, 'exp_As')
    
    strains = [exp_WT,exp_aS,exp_As]
    ax3 = plt.subplot(1,1,1)
    
    for experiment in strains[istrain]:
                
        bins = df.loc[experiment,'X']
        Mean = df.loc[experiment,'Y']
        
        
        bins = bins[np.nonzero(Mean)]
        Mean = Mean[np.nonzero(Mean)]
        
        
        if style == 'bar':
            ax3.bar(bins, Mean, widths2, linewidth=0,alpha=0.5)
        else:
            ax3.plot(bins, Mean, label=experiment)
            
        ax3.set_yscale('linear')
        ax3.set_xscale(scaling)
        ax3.set_title(field, fontsize = 20)
        ax3.legend()
        ax3.set_xlabel('Correlated track size', fontsize=20)
        ax3.set_ylabel('Counts', fontsize=20)
        ax3.grid()
        if scaling == 'log':
            x0,x1 = plt.gca().get_xlim()
            y0,y1 = plt.gca().get_ylim()
            y0 = 0
            x0 = 10**0
        else:
            x0,x1 = plt.gca().get_xlim()
            y0,y1 = plt.gca().get_ylim()
# ======================================================================================================================
#                                                   MAIN
# ======================================================================================================================

if __name__ == '__main__':
    
    scriptname = os.path.basename(sys.argv[0])
    
    # Reads the yaml document with scripts parameters
    yamlDoc = open("ScriptsInfos.yaml", 'r')
    inputs, _ = readYamlDoc(yamlDoc, scriptname)
    
    # load datashaders for several experiments to compare
    data = loadDatashaders(inputs)
    
    # Computes similarity index map for each experiment
    result = compute_similarity_map_all(data, inputs)
    
    # Quantify similarity and format results into dataframe
    df, widths2 = quantifySimilarity_to_df(result, inputs)

    # Display similarity histogram for one strain
    istrain = 0
    disp_SimQuantif_SingleResult(df, widths2, inputs, istrain = istrain, field = 'area', style = 'bar0', scaling = 'log', norm = False, nb_bins = 50)    
      
    # Display the mean similarity histogram for each strain
    disp_SimQuantif_results(df, widths2, inputs, field = 'area', style = 'bar0', scaling = 'log', norm = False, nb_bins = 50)
    
                    
    plt.gcf().tight_layout()
    plt.gca().set_xlim(0,1000)
    plt.gca().set_ylim(0,1000)
    plt.gca().set_aspect(1.0/plt.gca().get_data_ratio(), adjustable='box')
        
                
                
       


