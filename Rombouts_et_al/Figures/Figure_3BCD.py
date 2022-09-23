#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 16:42:07 2022

@author: legall
"""


''' 
Loads datashader maps generated for one set of experiment with Figure3_generate_tracks_datashader_maps_v20220207.py
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
from matplotlib import cm
from matplotlib.colors import ListedColormap
import pandas as pd
import pickle5 as pickle
import  matplotlib as mpl
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
    experiments = getInfoFromListOfDict(inputs, 'experiments')
    plot_width = 5880
    data = {}
    for experiment in experiments:
        
        print('\n *** Loading datashader, please wait... ***')
        data[experiment] , folderpath = load_datashader_maps(experiment, plot_width)

    print(' *** Datashader data successfully loaded ***\n')
    return data

def computeSimilarity(data, inputs):
    '''
    

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    inputs : TYPE
        DESCRIPTION.

    Returns
    -------
    result : TYPE
        DESCRIPTION.
    imA : TYPE
        DESCRIPTION.
    imB : TYPE
        DESCRIPTION.

    '''
    experiments = getInfoFromListOfDict(inputs, 'experiments')
    
    #%% ------ SIMIALRITY 
    result = {}
    imA = {}
    imB = {}
    for experiment in experiments:
        print('Processing experiment: ', experiment)
        classes = list(data[experiment][0].keys())         # List of population classes
        class_sel = [i for i in range(len(classes)) if classes[i] in ['scouts','swarms']]
        arr = data[experiment][0] # Dictionnary with tracks arrays per class
        
        from skimage.metrics import structural_similarity
        
        # Get tracks images and normalize them
        im1 = arr[classes[class_sel[0]]]
        im1 = (im1 - im1.min()) / (im1.max() - im1.min())
        im2 = arr[classes[class_sel[1]]]
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
        result[experiment] = abs(diff*(im1>im1.min()))
        result[experiment] = (result[experiment] * 2**16).astype("uint16")
       
        imA[experiment] = im1
        imB[experiment] = im2 
    return result, imA, imB

def dispSimilarity(result, imA, imB, experiment, inputs):
    '''
    

    Parameters
    ----------
    result : TYPE
        DESCRIPTION.
    imA : TYPE
        DESCRIPTION.
    imB : TYPE
        DESCRIPTION.
    experiment : TYPE
        DESCRIPTION.
    inputs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    xstart = getInfoFromListOfDict(inputs, 'xstart')
    ystart = getInfoFromListOfDict(inputs, 'ystart')
    Zoom_Width = getInfoFromListOfDict(inputs, 'Zoom_Width')
    
    # Make custom colormap
    cmp = cm.get_cmap('Greens', 512)
    newcolors = cmp(np.linspace(0,1, 512))
    newcolors[0, :] = np.array([1, 1, 1, 1])
    newcmp = ListedColormap(newcolors)
    
    # ------ PLOT SPECIFIC EXPERIMENTS CORRELATION RESULTS
    classes = list(data[experiment][0].keys())
    class_sel = [i for i in range(len(classes)) if classes[i] in ['scouts','swarms']]
    
    im1 = imA[experiment]
    im2 = imB[experiment]
    corr = result[experiment]
    corr = corr-corr.min()
    corr = corr/corr.max()
    
    if experiment == 'exp8':
        from scipy.ndimage import rotate
        im1 = rotate(im1, -90)
        im2 = rotate(im2, -90)
        corr = rotate(corr, -90)
        
    
    interp = 'antialiased'
    fig, axs = plt.subplots(1,3)
    axs[0] = plt.subplot(1,3,1)
    axs[0].matshow(im2, cmap=newcmp, vmin=0, vmax=im2.max(), interpolation=interp)
    axs[0].invert_yaxis()
    axs[0].set_title(classes[class_sel[1]], fontsize=20)
    axs[1] = plt.subplot(1,3,2, sharex=axs[0], sharey=axs[0])
    axs[1].matshow(im1, cmap=newcmp, vmin=0, vmax=im1.max(), interpolation=interp)
    axs[1].invert_yaxis()
    axs[1].set_title(classes[class_sel[0]], fontsize=20)
    axs[2] = plt.subplot(1,3,3, sharex=axs[0], sharey=axs[0])
    axs[2].invert_yaxis()
    axs[2].matshow(corr, cmap='seismic', vmin=-0.1, vmax=0.1, interpolation=interp)
    axs[2].set_title('Correlation map', fontsize=20)
    
    for i in range(3):
        for tick in axs[i].xaxis.get_major_ticks():
            tick.tick1line.set_visible(False)
            tick.tick2line.set_visible(False)
            tick.label1.set_visible(False)
            tick.label2.set_visible(False)
        for tick in axs[i].yaxis.get_major_ticks():
            tick.tick1line.set_visible(False)
            tick.tick2line.set_visible(False)
            tick.label1.set_visible(False)
            tick.label2.set_visible(False)
        axs[i].set_xticklabels([])
        axs[i].set_yticklabels([])
        
    # Run those two lines separatly after figure is shown!
    plt.tight_layout()
    plt.subplots_adjust(wspace = .01)
        
    
    axs[0].set_xlim(xstart, xstart+Zoom_Width)
    axs[0].set_ylim(ystart, ystart-Zoom_Width)
        
        
    plt.figure()
    ax1 = plt.subplot(1,3,1)
    ax1.imshow(im2, cmap=newcmp, vmin=0, vmax=im2.max())
    ax1.set_title(classes[class_sel[1]], fontsize=20)
    ax2 = plt.subplot(1,3,2, sharex=ax1, sharey=ax1)
    ax2.imshow(im1, cmap=newcmp, vmin=0, vmax=im1.max())
    ax2.set_title(classes[class_sel[0]], fontsize=20)
    ax3 = plt.subplot(1,3,3, sharex=ax1, sharey=ax1)
    ax3.imshow(corr, cmap='seismic', vmin=-0.1, vmax=0.1)
    ax3.set_title('Correlation map', fontsize=20)
    
    
def quantifySimilarity(result, inputs):
    '''
    

    Parameters
    ----------
    result : TYPE
        DESCRIPTION.
    inputs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    experiments = getInfoFromListOfDict(inputs, 'experiments')
    
    # - give colors to experiments plots
    colorList = ["blue","orange","green","magenta","yellow","red"]
    color_dic = {}
    color_dic = dict(zip(experiments, colorList[:len(experiments)]))
    
    # ------ SIMIALRITY QUANTIFICATION
    field = 'area' 
    style = 'bar0'
    scaling = 'log'
    norm = False
    nb_bins = 50
    
    for experiment in experiments:
  
        # =============================================================================
        #         ## COMPUTE STATS WITH CORRELATION MAP
        # =============================================================================
        label_result = label(result[experiment]>0, connectivity=result[experiment].ndim, background=0)
        props_result = regionprops_table(label_result,result[experiment], properties=[field], cache = True)
        df_props_result = pd.DataFrame(props_result)
        ax2 = plt.subplot(1,1,1)

        # # log-scaled bins
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
        
        # plot it!
        xdata = bins[:-1]
        ydata = hist_norm
        
        xdata = xdata[np.nonzero(ydata)]
        ydata = ydata[np.nonzero(ydata)]
        
        if style == 'bar':
            ax2.bar(bins[:-1], hist_norm, widths2, linewidth=0, facecolor=color_dic[experiment],alpha=0.5)
        else:
            ax2.plot(xdata, ydata, linewidth=2, color=color_dic[experiment])
        ax2.set_yscale('linear')
        ax2.set_xscale(scaling)
        ax2.set_title(field, fontsize = 20)
        ax2.legend(experiments)
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
        print('mean correlation : ', df_props_result[field].mean())
        
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
    
    # Computes similarity index maps
    result, imA, imB = computeSimilarity(data, inputs)
    
    # Display one experiment
    experiment = getInfoFromListOfDict(inputs, 'experiments')[0]
    dispSimilarity(result, imA, imB, experiment, inputs)

    # Quantify correlated tracks for all experiments
    quantifySimilarity(result, inputs)
    
    plt.gcf().tight_layout()
    plt.gca().set_xlim(0,1000)
    plt.gca().set_ylim(0,2000)
    plt.gca().set_aspect(1.0/plt.gca().get_data_ratio(), adjustable='box')
        
    
