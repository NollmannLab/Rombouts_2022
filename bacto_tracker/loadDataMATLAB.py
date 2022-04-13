#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 11:56:44 2020

@author: marcnol
"""

from scipy.io import loadmat
import os
import h5py
import numpy as np
from matplotlib.pylab import plt

import glob
import csv
import os

from matplotlib.pyplot import figure
import matplotlib 
matplotlib.rc('xtick', labelsize=15) 
matplotlib.rc('ytick', labelsize=15) 

rootFolder = '/home/marcnol/grey/rawData_2021/Experiment_50_Sara_WT/003_FastTimeLapse_RAMM_Test/Analyzed/'
tracksFolder = rootFolder+'Track_IDs/'

outputFigures = '/home/marcnol/gdrive/papers/2021.Methods_Sara/version_1/Figures/Figure_4/'

#%% loads all tracks in CSV format

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

    
#%% calculates lengths

track_lengths = list()
dt = 1 # in frames
for track in track_TXY:
    track_lengths.append(dt*len(track))

figure(figsize=(10, 8), dpi=300)
plt.hist(track_lengths, bins=range(dt,100*dt,10*dt),alpha=.7, label='N='+str(len(track_TXY)))
plt.xlabel('Track length, frames',fontsize = 20)
plt.ylabel('Counts', fontsize = 20)
plt.xlim([10,90])
plt.legend(fontsize = 20)
plt.savefig(outputFigures + "length_histogram.svg")

#%% calculates speeds

track_speeds = list()
dt = 0.5 # in secs
pixel_size = 0.1 # in microns

for track in track_TXY[0:100]:
    track_length = len(track)
    for idx in range(track_length-1):
        dx = track[idx,0] - track[idx+1,0]
        dy = track[idx,1] - track[idx+1,1]        
        v = pixel_size*np.sqrt( dx**2 + dy**2 )/dt
        track_speeds.append(v)

figure(figsize=(10, 8), dpi=300)
plt.hist(track_speeds, bins=range(0,10),alpha=0.7, label='N='+str(len(track_TXY)))
plt.xlabel('Velocity, microns / min',fontsize = 20)
plt.ylabel('Counts',fontsize = 20)
plt.xlim([0,10])
plt.legend(fontsize = 20)
plt.savefig(outputFigures + "speed_histogram.svg")


#%% examples of some tracks
figure(figsize=(10, 10), dpi=300)

for track in track_TXY[0:1000]:
    plt.plot(track[:,2],track[:,1],'-')

plt.ylim([0,5600])
plt.xlim([0,5600])
plt.xlabel('x, px',fontsize = 20)
plt.ylabel('y, px',fontsize = 20)
plt.savefig(outputFigures + "trajectories.svg")

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
