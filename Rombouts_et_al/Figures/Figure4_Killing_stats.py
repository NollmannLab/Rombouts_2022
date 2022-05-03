#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 09:52:45 2022

@author: legall
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import  matplotlib as mpl

#%%
Mutant = np.array(['A+S-', 'A-S+', 'WT', 'A+S-', 'A-S+', 'WT', 'A+S-', 'A-S+', 
                   'WT', 'WT', 'WT', 'WT', 'A+S-', 'A-S+'], dtype=object)

Kill_time = np.array([2262.0, 1015.0, 1104, 3856.0, 1981.0, 840, 2583.0, 2530.0, 185,
       1124, 1528, 673, 3381, 3110], dtype=object)

Control = np.array([4405.0, 6598.0, 2801, 3599.0, 3271.0, 4243, 2462.0, 3055.0,
                    6313, 1723, 3972, 3935, 4052, 4915], dtype=object)

Experiment = np.array(['exp12', 'exp11', 'exp6', 'exp51', 'exp54', 'exp8', 'exp53',
       'exp55', 'exp9', 'exp50', 'exp66', 'exp70', 'exp71', 'exp72'], dtype=object)

TimeStamp = np.array([33.9, 33.489, 40.898, 35.507, 37.7878, 37.25, 37.764,
       35.356, 35.33, 33.54, 35, 35.794, 31.711, 31.912], dtype=object)


#%% Remvove specific experiment from lists
exp_to_remove = np.array(['exp50'])
i_remove = np.where(exp_to_remove == Experiment)

Kill_time = np.delete(Kill_time, i_remove)
Experiment = np.delete(Experiment, i_remove)
Control = np.delete(Control, i_remove)
TimeStamp = np.delete(TimeStamp, i_remove)
Mutant = np.delete(Mutant, i_remove)

#%%
# Merge all data into one single array
data = np.vstack((Mutant,   Kill_time, Control,    Experiment, TimeStamp))
data = np.swapaxes(data,0,1)

# Convert it to dataframe and change data type to be readable by sns
df = pd.DataFrame(data, columns= ['Mutant',   'Kill_time', 'Control',    'Experiment', 'TimeStamp'])
df["Kill_time"] = df.Kill_time.astype(float)
df["Control"] = df.Control.astype(float)
df["TimeStamp"] = df.TimeStamp.astype(float)
#%%
# Format Dataframe to plot control and experiment side by side
df1 = df[['Mutant',   'Kill_time',    'Experiment', 'TimeStamp']]
df2 = df[['Mutant',   'Control',    'Experiment', 'TimeStamp']]

df_cat = pd.concat([df1, df2], ignore_index=True)
df_cat = df_cat.replace(np.nan,0)
df_cat["Time"] = df_cat["Kill_time"] + df_cat["Control"]    
df_cat["Time"] = df_cat["Time"] * df_cat["TimeStamp"]   # Time in seconds
df_cat["Time"] = df_cat["Time"] /3600                   # Time in hours
df_cat["Rate"] = 1/df_cat["Time"]
df_cat["Rate_log"] = np.log(df_cat["Rate"])
df_cat['Control'] = df_cat['Control'] > 0
df_cat.drop('Kill_time', axis=1, inplace=True)  # No longer needed

#%% Removes some experiments
df_cat = df_cat[ ~( (df_cat['Control'] == True) & (df_cat['Experiment'] == 'exp50') )]
df_cat = df_cat[ ~( (df_cat['Control'] == False) & (df_cat['Experiment'] == 'exp66') )]

#%% Sort dataframe index according to mutant custom list
ordered_classes = ['WT', 'A-S+', 'A+S-']

df_list = []

for i in ordered_classes:
   df_list.append(df_cat[df_cat['Mutant']==i])

ordered_df = pd.concat(df_list)

#%%     Display results
fig5, ax5 = plt.subplots(figsize = (9, 9))
exp_df = ordered_df[ordered_df['Control'] == False]
cont_df = ordered_df[ordered_df['Control'] == True]

# Make sure to remove the 'facecolor': 'w' property here, otherwise
# the palette gets overrided
boxprops = {'edgecolor': 'k', 'linewidth': 0, 'alpha': 0.15, 'facecolor': 'r'}
lineprops = {'color': 'k', 'linewidth': 0}

boxplot_kwargs = {'boxprops': boxprops, 'medianprops': lineprops,
                  'whiskerprops': lineprops, 'capprops': lineprops,
                  'width': 1}

ax5 = sns.boxplot(x="Mutant", y="Time", data=cont_df,  **boxplot_kwargs)


# Make sure to remove the 'facecolor': 'w' property here, otherwise
# the palette gets overrided
boxprops = {'alpha': 0.75, 'facecolor': 'b'}
lineprops = {'color': 'k', 'linewidth': 1}

boxplot_kwargs = {'boxprops': boxprops, 'medianprops': lineprops,
                  'whiskerprops': lineprops, 'capprops': lineprops,
                  'width': 0.3}

ax5 = sns.boxplot(x="Mutant", y="Time", data=exp_df,  **boxplot_kwargs)


ax5.set(ylim=(0, 50))
ax5.set_aspect(1./ax5.get_data_ratio())


[ax5.axvline(x+.5,color='k') for x in ax5.get_xticks()]
ax5.grid(True, axis = 'y', linestyle='--')
ax5.set_axisbelow(True)

for axis in ['top','bottom','left','right']:
    ax5.spines[axis].set_linewidth(2)

plt.tight_layout()
plt.box(on=True)
plt.ylabel('Fluorescence decay time, hours', fontsize = 18)
plt.xlabel('Mutants', fontsize = 18)