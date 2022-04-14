#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 11:56:44 2020

@author: marcnol
"""

import os
import numpy as np

import glob
import csv

import matplotlib 


# replace the following with the folders with your data
rootFolder = '/home/marcnol/grey/rawData_2021/Experiment_50_Sara_WT/003_FastTimeLapse_RAMM_Test/Analyzed/'

#%% loads all tracks in CSV format

tracksFolder = rootFolder+'Track_IDs/'
files = [x for x in glob.glob(tracksFolder +'*.csv')]
print(f'{len(files)} files to process')

tracks = list()
# col 1: frames, col2: cell IDs in those frames.
for x in files:
    file = open(x)
    csvreader = csv.reader(file)
    rows = []
    for row in csvreader:
        rows.append(row)
    number_points = len(rows)
    track = np.zeros((number_points ,2))
    for i, row in enumerate(rows):
        for idx in range(len(row)):
            track[i,idx] = float(row[idx])

    tracks.append(track)
    print(f'loaded and processed file: {os.path.basename(x)} with {number_points} data points')

print(f'loaded {len(tracks)} tracks')

#%% converts mask files to dictionary

imageFolder= rootFolder + 'masks/'

files = [x for x in glob.glob(imageFolder + '*.csv')]
print(f"{len(files)} files found ")

keys = [int(os.path.basename(file).split("_")[1].split('.')[0]) for file in files]  

mask_dict = {}
# mask file format: maskID, x-coord, y-coord

for key,file in zip(keys,files):
    
    newfile = open(file)
    csvreader = csv.reader(newfile)
    rows = []
    for row in csvreader:
        rows.append(row)
    number_points = len(rows)
    track = np.zeros((number_points,3))
    for i, row in enumerate(rows):
        for idx in range(len(row)):
            track[i,idx] = float(row[idx])

    mask_dict[key]={}
    mask_dict[key]['track_xy'] = track
    
    print(f'loaded and processed file: {os.path.basename(file)} for key {key} with {number_points} data points')
    
#%% builds tracks

track_TXY = list()
# track_TXY: each list item is an array with rows containing (frame number, x-centroid, y-centroid)

# iterate over tracks
for track in tracks:
    track_TXY_array = np.zeros((len(track),3))
    
    # iterate over frames
    for idx in range(len(track)):
        frame_number = int(track[idx,0])
        mask_number = int(track[idx,1])
        item = mask_dict[frame_number]['track_xy'][mask_number-1,:]

        track_TXY_array[idx,:] = (frame_number,item[1],item[2])
        
    track_TXY.append(track_TXY_array)

print('done')
    
#%% saves in csv format

# track number, frame number, x centroid (px), y centroid (px)

# converts track_TXY to list
output_data = list()
for track_number, track in enumerate(track_TXY):
    for row_number in range(len(track)):
        row = [track_number, int(track[row_number,0]), track[row_number,1], track[row_number,2]] 
        output_data.append(row) 

myFile = open(rootFolder+'track_TXY.csv', 'w')

with myFile:
   writer = csv.writer(myFile)
   writer.writerows(output_data)
   
