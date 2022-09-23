#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 16:08:12 2022

@author: legall
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import  matplotlib as mpl
from statannotations.Annotator import Annotator
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
exp_to_remove = np.array(['exp50']) # --> 
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


# fig, ax = plt.subplots()
# ax = sns.violinplot(x='Mutant', y="Kill_time", data=df)
# fig2, ax2 = plt.subplots()
# ax2 = sns.violinplot(x='Mutant', y="Control", data=df)

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

#%% Removes some values

df_cat = df_cat[ ~( (df_cat['Control'] == True) & (df_cat['Experiment'] == 'exp50') )]
# df_cat = df_cat[ ~( (df_cat['Control'] == True) & (df_cat['Experiment'] == 'exp66') )]
df_cat = df_cat[ ~( (df_cat['Control'] == False) & (df_cat['Experiment'] == 'exp66') )]

#%%     Sort dataframe index according to mutant custom list
ordered_classes = ['WT', 'A-S+', 'A+S-']

df_list = []

for i in ordered_classes:
   df_list.append(df_cat[df_cat['Mutant']==i])

ordered_df = pd.concat(df_list)


#%% Check if synergy
WT_time = ordered_df[ (ordered_df['Mutant'] == 'WT') & (ordered_df['Control'] == False)]['Time']
aS_time = ordered_df[ (ordered_df['Mutant'] == 'A-S+') & (ordered_df['Control'] == False)]['Time']
As_time = ordered_df[ (ordered_df['Mutant'] == 'A+S-') & (ordered_df['Control'] == False)]['Time']

from scipy.spatial.distance import cdist
def ratio(a, b):
    return 1/(1/a+1/b)

d = cdist(np.array(aS_time)[:, np.newaxis], np.array(As_time)[:, np.newaxis], lambda a, b: ratio(a[0], b[0]))
Mutants = d.flatten()

combined =  pd.DataFrame({"Mutant":"Combined", "Time": Mutants})
df = pd.concat([ordered_df[ordered_df["Control"]==False][["Mutant", "Time"]], combined])


x = "Mutant"
y = "Time"
order = ['A-S+', 'A+S-', 'Combined', 'WT']
configuration = {'test':'t-test_ind',   #'t-test_ind' or 'Mann-Whitney'
                 'comparisons_correction':None,
                 'text_format':'star'}

ax = sns.boxplot(data=df, x=x, y=y, order=order)

annot = Annotator(ax, [("A-S+", "WT"), ("A+S-", "WT"), ("Combined", "WT"), ("A+S-", "A-S+")], data=df, x=x, y=y, order=order)
annot.configure(**configuration, loc='inside', verbose=2)
annot.apply_test()
ax, test_results = annot.annotate()

bplot=sns.stripplot(data=df, x=x, y=y, order=order,
                   jitter=True, 
                   size=10,
                   alpha=0.5,
                   linewidth=1)

# from scipy import stats
# stats.ttest_ind(df[df["Mutant"]=='WT']['Time'], df[df["Mutant"]=='Combined']['Time'], equal_var=False)
# stats.ttest_ind(df[df["Mutant"]=='Combined']['Time'], df[df["Mutant"]=='WT']['Time'], equal_var=False)
#%%
# # Build WT array to the dimension
# WT_filled = np.empty(Mutants.shape[0])
# WT_filled.fill(np.nan)
# WT_filled[:WT_time.shape[0]] = WT_time

# # Build aS array to the dimension
# aS_filled = np.empty(Mutants.shape[0])
# aS_filled.fill(np.nan)
# aS_filled[:aS_time.shape[0]] = aS_time

# # Build As array to the dimension
# As_filled = np.empty(Mutants.shape[0])
# As_filled.fill(np.nan)
# As_filled[:As_time.shape[0]] = As_time

# df_ratios = pd.DataFrame( {"WT":WT_filled, "Combined":Mutants, "A-S+":aS_filled, "A+S-":As_filled})

# # make boxplot with Seaborn
# bplot=sns.boxplot(data=df_ratios, 
#                  width=0.3)
 
# # add stripplot to boxplot with Seaborn
# bplot=sns.stripplot(data=df_ratios, 
#                    jitter=True, 
#                    size=10,
#                    alpha=0.5,
#                    linewidth=1)

# # from statannot import add_stat_annotation
# # # statistical annotation
# # add_stat_annotation(bplot, data=df_ratios,
# #                     box_pairs=[("WT", "Combined"), ("WT", "A-S+"), ("WT", "A+S-")],
# #                     test='Mann-Whitney', text_format='star', loc='inside', verbose=2)

# from statannotations.Annotator import Annotator
# box_pairs=[("WT", "Combined"), ("WT", "A-S+"), ("WT", "A+S-")]
# annotator = Annotator(bplot, box_pairs, data=df_ratios)
# annotator.configure(test='t-test_ind', text_format='star', loc='inside')
# annotator.apply_and_annotate()

# from scipy.stats import mannwhitneyu
# U1, p = mannwhitneyu(df_ratios['WT'], df_ratios['Combined'])
# print(U1)
# from scipy import stats
# stats.ttest_ind(df_ratios['WT'].dropna(), df_ratios['Combined'], equal_var=False)