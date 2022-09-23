#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 09:49:19 2021

@author: legall
"""


"""
Converts Matlab mat files to python dataframe format.

1. Loads tracking data (mat files, MSD & multiscale segmentation)
2. Computes additional features (cell growth, speed, some stats, ...)
3. Change data format
4. Exports as a single dataframe python file.

Usage:
-  modify PreprocessingInfos.yaml file to specify the desired experiment to process
-  modify DATA PATH in function "experimentDataPath" to match the location of you matlab mat files
"""



#%% IMPORTS ###################################################################
import os
import yaml
import sys
import numpy as np
import numpy.matlib
from scipy.io import loadmat
import pandas as pd
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


def loadmatfiles(folderpath, remove_coli):
    '''
    Loads all matlab .mat files (matrices) within a folder into a single dictionnary. Each matrix is a feature 

    :param folderpath: folder path that contains mat files
    :return: my_dict: dictionnary with each field containing one matrix
    '''
    print('\n *** Data loading ***')
    if remove_coli :
        #- Loads list of rows corresponding to Ecoli cell that need to be deleted
        rows_to_del = np.load(folderpath+'Ecoli_rows_to_del.npy')

    my_dict = {}
    for file in os.listdir(folderpath):
        if file.endswith(".mat"):
            
            data = loadmat(folderpath+file)
            print(folderpath+file)
            filename, ext = os.path.splitext(file)
            print(filename, ' loading...')
            data = data.get(filename).astype('float32')
            if filename=='Total_Col' or filename=='Total_Row':
                col2add = np.zeros((data.shape[0],1))
                col2add[col2add == 0] = 'nan'
                data = np.concatenate((data, col2add), axis=1)
            else:
                data = np.concatenate((data, data[:,-1][:,None]), axis=1)
            if remove_coli :
                data = np.delete(data, np.where(rows_to_del), axis=0)
                
            my_dict[filename] = data
            print(filename, ' loaded')
    print('--> All files successfully loaded.')

    # Add a field with tracks ID to control results
    ID = np.arange(my_dict.get('Total_tracks').shape[0])  # Create an array with as many tracks ID as there are tracks
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
    track_length = np.matlib.repmat(track_length, my_dict.get('Total_tracks').shape[1], 1)  # reshape track_length to match tracks length
    track_length = np.swapaxes(track_length, 0, 1)  # swap row and col
    my_dict['track_length'] = track_length
    print('--> New fields TrackID and TimeStamp successfully added \n')

    return my_dict


def compute_field_stats(data_dict, field = 'Total_Density') :
    '''Takes a dictionnary and a field to compute statistics
    :param data_dict: a dictionnary with different matrices
    :param field is a
    :return: a dictionnary with an additional field for the stats
    '''
    # - compute stats
    field_data = data_dict[field]
    area = data_dict['Total_Area']  # Field used to discriminate elements that are not detected cells (area = 0)
    mean = np.nanmean(np.where(area != 0, field_data, np.nan), axis = 1)
    min = np.nanmin(np.where(area != 0, field_data, np.nan), axis = 1)
    max = np.nanmax(np.where(area != 0, field_data, np.nan), axis = 1)
    std = np.nanstd(np.where(area != 0, field_data, np.nan), axis = 1)

    # - reshape stats values to match input dictionnary shape
    mean_2 = np.matlib.repmat(mean, field_data.shape[1], 1)  # reshape stats to match tracks length
    mean_2 = np.swapaxes(mean_2, 0, 1)  # swap row and col
    data_dict[field+'_mean'] = mean_2
    min_2 = np.matlib.repmat(min, field_data.shape[1], 1)  # reshape stats to match tracks length
    min_2 = np.swapaxes(min_2, 0, 1)  # swap row and col
    data_dict[field + '_min'] = min_2
    max_2 = np.matlib.repmat(max, field_data.shape[1], 1)  # reshape stats to match tracks length
    max_2 = np.swapaxes(max_2, 0, 1)  # swap row and col
    data_dict[field + '_max'] = max_2
    std_2 = np.matlib.repmat(std, field_data.shape[1], 1)  # reshape stats to match tracks length
    std_2 = np.swapaxes(std_2, 0, 1)  # swap row and col
    data_dict[field + '_std'] = std_2

    print('Statistics were added on field : ', field)

    return data_dict

def compute_instant_speed_and_direction(data_dict, win_size = 1,**kwargs):
    '''Takes a dictionnary and a list of fields and outputs those into a matrix
    :param data_dict: a dictionnary with different matrices
    :param win_size is a value for the Nth point in the track to take to compute speed
    :return: a dictionnary with an additional field "instant_speed"
    '''

    #- compute speed from previous and next points locations
    frame_num = data_dict['TimeStamp']
    rows = data_dict['Total_Row']
    cols = data_dict['Total_Col']

    frame_num = frame_num.astype(np.double)
    rows = rows.astype(np.double)
    cols = cols.astype(np.double)

    drows = np.roll(rows, -win_size, axis=1) - rows
    dcols = np.roll(cols, -win_size, axis=1) - cols
    dt = np.roll(frame_num, -win_size, axis=1) - frame_num
    instant_speed_to = np.sqrt(drows**2+dcols**2) / dt
    disp_to = np.sqrt(drows**2+dcols**2)
    direction_to = np.degrees(np.arctan2(dcols, drows))

    drows = np.roll(rows, win_size, axis=1) - rows
    dcols = np.roll(cols, win_size, axis=1) - cols
    dt = frame_num - np.roll(frame_num, win_size, axis=1)
    instant_speed_from = np.sqrt(drows ** 2 + dcols ** 2) / dt
    disp_from = np.sqrt(drows ** 2 + dcols ** 2)
    direction_from = np.degrees(np.arctan2(dcols, drows))

    #- set area to 0 where instant speed is erronous (first/last columns and where there weren't any prev/next points
    #- Note : area is the criteria I use later on to select valid datapoints (cells with null area have no meaning)
    area = data_dict['Total_Area']
    valid = np.zeros(area.shape)
    valid[np.where(area > 0)] = 1

    valid = valid * np.roll(valid, -win_size, axis=1) * np.roll(valid, win_size, axis=1)
    valid[:, :win_size-1] = 0
    valid[:, -win_size:] = 0

    area = area * valid
    instant_speed_from = instant_speed_from * valid
    instant_speed_to = instant_speed_to * valid
    disp_to = disp_to  * valid
    disp_from = disp_from * valid
    direction_to = direction_to * valid
    direction_from = direction_from * valid

    direction_to = np.nan_to_num(direction_to)
    direction_from = np.nan_to_num(direction_from)


    data_dict["direction_from"] = direction_from
    data_dict["direction_to"] = direction_to

    #- Computes relative direction angle (turns right or left)
    # converts atan2 results into
    direction_to = direction_to + np.logical_and(direction_to < 0, direction_to > -np.inf) * 360
    direction_from = direction_from + np.logical_and(direction_from < 0, direction_from > -np.inf) * 360
    # Compute relative angle from -180° to +180°
    rel = 180 - direction_to + direction_from
    rel = rel + np.logical_and(rel < -180, rel > -np.inf) * 360
    rel = rel - np.logical_and(rel > 180, rel < np.inf) * 360

    data_dict["direction_rel"] = rel
    data_dict["direction_abs"] = abs(rel)
    data_dict["disp_to"] = disp_to
    data_dict["disp_from"] = disp_from
    data_dict["instant_speed_from"] = instant_speed_from
    data_dict["instant_speed_to"] = instant_speed_to
    data_dict["instant_speed"] = (instant_speed_from+instant_speed_to)/2
    data_dict["Total_Area"] = area

    print('Fields instant_speed, displacement and directions were added to data_dict')

    return data_dict

def dict_2_mat(data_dict, fields={}, **kwargs):
    '''Takes a dictionnary and a list of fields and outputs those into a matrix
    :param data_dict: a dictionnary with different matrices
    :return: a matrix for dimensionnality reduction
    '''
    print('\n *** Breaking down track into single data points ***')
    if not fields:  # if provided dictionnary fields is empty
        nb_k = data_dict.__len__()  # number of fields
        nb_points = np.count_nonzero(data_dict.get('Total_tracks'))  # number of cells detected

        orderedNames = data_dict.keys()
        dataMatrix_All = np.array([data_dict[i] for i in orderedNames]) # 3D matrix with rows=fields, col=track, slice=frame_num
        dataMatrix_All_flat = np.reshape(dataMatrix_All, (nb_k, dataMatrix_All.shape[1]*dataMatrix_All.shape[2])) # 2D matrix with rows=fields, col=cell
        dataMatrix_All_flat_swap = np.swapaxes(dataMatrix_All_flat, 0, 1)  # swap row and col
        i_area = list(data_dict.keys()).index("Total_Area") #   finds area key index
        data = dataMatrix_All_flat_swap[~(dataMatrix_All_flat_swap[:,i_area] == 0)]  # removes points with empty area
    else:
        print('Not yet implemented...')
    print('--> Done \n')

    return (data, list(orderedNames))

def fit_linear(data_dict, field_y='Total_Area', field_x='TimeStamp') :
    '''
    Takes a dictionnary and fits selected fields with linear regression.
    '''
    my = data_dict[field_y]
    mx = data_dict[field_x]
    slope = []
    for row_y, row_x in zip(my, mx):
        out_y = row_y[np.nonzero(row_y)]
        out_x = row_x[np.nonzero(row_y)]
        if out_y.size == 0 or out_x.size ==0:
            a = np.nan
        else:
            a, b = np.polyfit(out_x, out_y, 1) # y = ax + b.
        slope.append(a)
    return slope

def compute_field_log10(data_dict, fields_to_log):
    '''
    Takes the log10 of selected fields.
    Replaces -Inf from log10 with zero.
    '''
    for f in fields_to_log:
        data = data_dict[f]
        data_dict[f] = np.log10(data)
        data_dict[f][np.isneginf(data_dict[f])] = 0 # replace -Inf from log10 with zero
    
    return data_dict

def mergeMatlabfiles2Dataframe(folderpath, inputs):
    # LOADS TRACKS DATA INTO DICTIONNARY - ############################################################################
    remove_coli = False # in case user wants to remove some EC cells wrongly identified as MX cells.
    data_dict = loadmatfiles(folderpath, remove_coli)
     
    # LOADS ADDTIONNAL FEATURES AND PUT THEM INTO DICTIONNARY - ######################################################
    data = np.load(folderpath+'MSD_n_5frames.npy').astype('float32')
    data = np.concatenate((data, data[:,-1][:,None]), axis=1)
    data_dict['MSD_n_5frames'] = data #np.load(folderpath+'MSD_n_5frames.npy').astype('float32')
    del data
    data = np.load(folderpath+'MSD_D_5frames.npy').astype('float32')
    data = np.concatenate((data, data[:,-1][:,None]), axis=1)
    data_dict['MSD_D_5frames'] = data #np.load(folderpath+'MSD_D_5frames.npy').astype('float32')
    del data
    
    # COMPUTES CELL GROWTH - #########################################################################################
    print('\n *** Computing growth rate...')
    growth = np.array(fit_linear(data_dict, field_y='Total_Backbone', field_x='TimeStamp'))
    growth = np.matlib.repmat(growth, data_dict[next(iter(data_dict))].shape[1],
                              1)  # reshape stats to match tracks length
    growth = np.swapaxes(growth, 0, 1)  # swap row and col
    data_dict['growth_rate'] = growth
    print('--> Done \n')

    # COMPUTES INSTANTANEOUS SPEED (PIX/FRAME) FROM PREVIOUS AND NEXT DATAPOINTS-  ###################################
    # - COMPUTES MOVEMENT DIRECTION (FROM AND TO) FROM PREVIOUS AND NEXT DATAPOINTS
    win_size = getInfoFromListOfDict(inputs, 'win_size')
    data_dict = compute_instant_speed_and_direction(data_dict, win_size)

    # - COMPUTES FEATURES MEAN, MIN, MAX AND VARIATION ALONG TRACKS - ###################################################
    data_dict = compute_field_stats(data_dict, 'Total_Density')
    data_dict = compute_field_stats(data_dict, 'instant_speed')
    data_dict = compute_field_stats(data_dict, 'Total_Area')
    data_dict = compute_field_stats(data_dict, 'Total_Backbone')
    data_dict = compute_field_stats(data_dict, 'direction_rel')

    # TAKES THE LOG10 OF SOME SELECTED FIELDS - #############################
    fields_to_log = ['Total_Density','Total_Area','instant_speed','instant_speed_from','instant_speed_to','Total_DensityEc','Total_DensityEcMyxo']
    print('\n *** Computing log10 ...')
    data_dict = compute_field_log10(data_dict, fields_to_log)
    print('--> Done \n')
    
    
    # BREAKS TRACKS INTO SINGLE DATA POINTS. OUTPUTS A MATRIX "DATA" AND A LIST "FIELDS" WHICH TELLS WHAT DATA COLUMNS
    # CORRESPOND TO. - #################################################################################################
    (data, fields) = dict_2_mat(data_dict)

    # Get clusters info
    print('\n *** Appending clusters data...')

    import pickle5 as pickle
    
    fh = folderpath+"dataframe/"+"clusters_df.pkl"
    with open(fh, "rb") as fh:
      clusters_df = pickle.load(fh)
    fh = folderpath+"dataframe/"+"cells_df.pkl"
    with open(fh, "rb") as fh:
      cells_df = pickle.load(fh)
    
    # Load cluster dataframes
    clusters_df = clusters_df.rename(columns = {'area':'clust_area'})
    cells_df = cells_df.merge(clusters_df[['clust_area', 'cluster_id','Frame_num']],
                              left_on=['cluster_id','Frame_num'], 
                              right_on=['cluster_id','Frame_num'], 
                              how='left')
    
    # Get clusters key info to find matching indices between datasets
    cell_ids = cells_df['cell_id']+1
    Frame_num = cells_df['Frame_num']-1
    
    # Format both relevant data into 2D array for comparison
    A = np.append(data[:,fields.index("Total_tracks")][:,None], data[:,fields.index("TimeStamp")][:,None], axis=1).astype('int64')
    B = np.append(cell_ids[:, None], Frame_num[:, None], axis=1)
    
    # Find common rows indices
    dims = B.max(0)+1
    X1D = np.ravel_multi_index(B.T,dims)
    searched_valuesID = np.ravel_multi_index(A.T,dims)
    sidx = X1D.argsort()
    indices = sidx[np.searchsorted(X1D,searched_valuesID,sorter=sidx)]
    
    # Merge cluster and cells dataframes info
    clusters_df = clusters_df.rename(columns = {'area':'clust_area',
                                            'orientation':'clust_orientation',
                                            'perimeter':'clust_perimeter',
                                            'eccentricity':'clust_eccentricity',
                                            'solidity':'clust_solidity',
                                            })


    cells_df = cells_df.merge(clusters_df[['clust_area', 'cluster_id','Frame_num',
                                            'clust_perimeter','clust_eccentricity',
                                            'clust_solidity','clust_orientation']],
                              left_on=['cluster_id','Frame_num'], 
                              right_on=['cluster_id','Frame_num'], 
                              how='left')
    # Store new clusters features into "data"
    fields_to_exp = list(cells_df.columns.values[3:-1])
    data_to_exp = cells_df[cells_df.columns[3:-1]].values[indices,:]

    fields = fields + fields_to_exp
    data = np.append(data, data_to_exp, axis=1)
    df = pd.DataFrame(data, columns=fields)
    
    return df
    
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

    # Loads all individual files (matlab, MSD, clusters) and merge everything into one single dataframe
    df = mergeMatlabfiles2Dataframe(folderpath, inputs)
    print('--> Done \n')
    
    #%% SAVES DATA AS DATAFRAME FORMAT - #############################
    saveResults = getInfoFromListOfDict(outputs, 'saveResults')
    if saveResults:
        print('Saving')
        win_size = getInfoFromListOfDict(inputs, 'win_size')
        df.to_pickle(folderpath+'dataframe/'+"all_features_WithNAN_win_size_"+str(win_size)+'.pkl')
        print('Data saved to dataframe format')
    
    