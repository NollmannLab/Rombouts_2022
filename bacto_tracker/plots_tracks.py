#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 11:16:55 2022

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


# replace the following with the output folder
outputFigures = '/home/marcnol/gdrive/papers/2021.Methods_Sara/version_1/Figures/Figure_4/'


#%% loads CSV track file

# replace the following with the folders with your data
rootFolder = '/home/marcnol/grey/rawData_2021/Experiment_50_Sara_WT/003_FastTimeLapse_RAMM_Test/Analyzed/'

myFile = rootFolder+'track_TXY.csv'

with open(myFile, newline='') as f:
    reader = csv.reader(f)
    data = list(reader)

# track number, frame number, x centroid (px), y centroid (px)

print(data)

#%%
# track_TXY: each list item is an array with rows containing (frame number, x-centroid, y-centroid)


# creates list with time, x, y
track_ID='0'
new_track = list()
track_list = list()
for row in data:
    new_track_ID = row[0]

    if new_track_ID != track_ID:
        track_list.append(new_track) 
        new_track = list()
        track_ID= new_track_ID

    new_track.append(row[1:4])
        
#%%    
# converts this list into track_TXY format
track_TXY = list()
for track in track_list:
    track_TXY_array = np.zeros((len(track),3))
    
    # iterate over frames
    for idx, row in enumerate(track):
        frame_number = int(row[0])
        x = float(row[1])
        y = float(row[2])        
        track_TXY_array[idx,:] = (frame_number,x,y)
    
    track_TXY.append(track_TXY_array)

print('done)')


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

