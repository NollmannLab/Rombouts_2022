#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 16:10:44 2022

@author: legall
"""

''' 
Loads dataframe file with all features for one experiment and computes the 
2D density histograms of trajectories using datashader. Export those maps.

Usage:
- modify DATA PATH to match the dataframe file of the experiment

'''

#%% IMPORTS ###################################################################
import os
import yaml
import sys
import numpy as np
import datashader as ds
import pandas as pd
import pickle5 as pickle
import time

#%%=====================================================================================================================
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
    print(' *** Tracks data successfully loaded ***\n')
    return data, folderpath


def formatDataframe(data):
    # -- Modify fields
    data['Nb_cells_mean'] = data.groupby(['TrackID'])['Nb_cells'].transform('mean') 
    
    
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
    # -- Sort dataframe
    # Sort dataframe rows by TrackID then Time
    data = data.sort_values(["TrackID", "TimeStamp"], ascending = (True, True))  
    
    # Find ends of tracks (change of track ID) and add NAN to be recognized as individual tracks with Datashader
    # Datashader recognizes ends of tracks when X and Y coord = NAN
    data["diff"] = data["TrackID"].diff(periods=-1)
    RowsToCopy= data.loc[(data["diff"] != 0) , :].copy()
    
    # Set values to np.nan in fake row except in ['TimeStamp','TrackID']
    cols = [col for col in RowsToCopy.columns if col not in ['TimeStamp','TrackID']]
    RowsToCopy[cols] = np.nan
    RowsToCopy['TimeStamp'] = RowsToCopy['TimeStamp']+1
    RowsToCopy['Total_Density'] = -1
    RowsToCopy['Nb_cells'] = -1
    data = pd.concat([data, RowsToCopy])

    # Sort again dataframe rows by TrackID then Time
    data = data.sort_values(["TrackID", "TimeStamp"], ascending = (True, True))  
    
    return data

def generateClassesTrackMaps(data, inputs):
    print('\n *** Computing Datashader matrices, please wait... ***')
    t = time.time()
    classes = ['scouts','loners', 'herd','rafts','swarms','all']
    class_sel = [0, 1, 2, 3, 4, 5]
    arr = {}
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
            density_max = 4.5
            Nb_cells_min = 3
            Nb_cells_max = 600
        elif classes[class_sel[i]] == 'swarms' :
            density_min = 0
            density_max = 4.5
            Nb_cells_min = 601
            Nb_cells_max = 10**6
        elif classes[class_sel[i]] == 'all' :
            density_min = 0
            density_max = 10**6
            Nb_cells_min = 0
            Nb_cells_max = 10**6
        speed_min = getInfoFromListOfDict(inputs, 'speed_min')
        speed_max = getInfoFromListOfDict(inputs, 'speed_max')
        distance_norm_min = getInfoFromListOfDict(inputs, 'distance_norm_min')
        distance_norm_max = getInfoFromListOfDict(inputs, 'distance_norm_max')
        Mask_width = getInfoFromListOfDict(inputs, 'Mask_width')
        t1 = getInfoFromListOfDict(inputs, 't1')
        t2 = getInfoFromListOfDict(inputs, 't2')
        
        if experiment == 'exp71':
            # This experiment has misconnected tracks at frame 389. These two lines break trajectories at this frame
            ind_to_break = data.index[ (data['TimeStamp'] >= 389) & (data['TimeStamp'] <= 389)].tolist()
            data.loc[ind_to_break, ['Total_Col','Total_Row']] = np.nan
            print('!!! Exp71 exception : Trajectories broke at frame 389. !!!')
        
        ind = data.index[ (data['Total_Density_mean'] >= 10**density_min) & (data['Total_Density_mean'] <= 10**density_max)
                         & (data['Nb_cells_mean'] >= Nb_cells_min) & (data['Nb_cells_mean'] <= Nb_cells_max)
                         & (data['track_length'] >= 0) & (data['track_length'] <= 70000)
                         & (data['instant_speed'] >= speed_min) & (data['instant_speed'] <= speed_max)
                         & (data['distance_norm'] >= distance_norm_min) & (data['distance_norm'] <= distance_norm_max)
                         & (data['Total_Col'] >= Mask_width) & (data['Total_Col'] <= 5880- Mask_width)
                         & (data['Total_Row'] >= Mask_width) & (data['Total_Row'] <= 5880- Mask_width)
                         & (data['TimeStamp'] >= t1) & (data['TimeStamp'] <= t2)].tolist()
        ind_endTracks = data.index[ (data['Total_Density'] == -1) & (data['Nb_cells'] == -1)].tolist()
        ind = np.sort(ind+ind_endTracks)

        
        #% ------ DATASHADER
        # Get fields
        if experiment == 'exp8' or experiment =='exp9':
            X_field = 'Total_Col'
            Y_field = 'Total_Row'
        else:
            X_field = 'Total_Row'
            Y_field = 'Total_Col'
        C_field = 'TimeStamp'
        
        df = data.loc[ind, [X_field, Y_field, C_field]]
        
        # ------ BUILD DATASHADER AGGREGATOR
        plot_width = getInfoFromListOfDict(inputs, 'plot_width')
        plot_height = getInfoFromListOfDict(inputs, 'plot_height')
        
        # # Default plot ranges:
        x_range = (0 , plot_width-1)
        y_range = (0 , plot_height-1)
        cvs = ds.Canvas(x_range=x_range, y_range=y_range, plot_width=round(plot_width/1), plot_height=round(plot_height/1))
        
        agg = cvs.line(df,X_field, Y_field,  ds.count()) #<-- USE THIS ONE FOR TRACK SIMILARITY
        arr[classes[class_sel[i]]] = agg.data   #img.values     #np.dot(ds_array[...,:3], [0.299, 0.587, 0.114])
    
    print('--> Done \n')
    print('Elapsed: %s' % (time.time() - t))
    
    
    return arr
# ======================================================================================================================
#                                                   MAIN
# ======================================================================================================================

if __name__ == '__main__':


    scriptname = os.path.basename(sys.argv[0])
    
    # Reads the yaml document with scripts parameters
    yamlDoc = open("ScriptsInfos.yaml", 'r')
    
    
    inputs, outputs = readYamlDoc(yamlDoc, scriptname)

    
    # load dataframe
    experiment = getInfoFromListOfDict(inputs, 'experiment')
    data, folderpath = load_dataframe_exp(experiment)
    
    print(data.columns.tolist())
    
    # format the data (add new fields and clean it a bit)
    data = formatDataframe(data)  
   
    # generate maps of trajectories for each population class
    arr = generateClassesTrackMaps(data, inputs)
    
    
    # ------ Save output datashader matrices
    SaveResults = getInfoFromListOfDict(outputs, 'SaveResults')
    plot_width = getInfoFromListOfDict(inputs, 'plot_width')
    if SaveResults:
        print('--> Saving file, please wait ...\n')
        
        folder_out = folderpath+'datashader/'
        if not os.path.exists(folder_out):
            # Make root destination folder
            os.makedirs(folder_out)
        file_to_write = open(folder_out+ experiment+'_tracks_classes' +'_'+str(plot_width)+ '_v2.pkl', "wb")
        pickle.dump(arr, file_to_write)
        file_to_write.close()
        
        print('--> Saving done \n')
        
    
    # Use this script to  continue  -->       Figure3_load_and_analyse_tracks_datashader_maps_v20210901.py
  
    
    
