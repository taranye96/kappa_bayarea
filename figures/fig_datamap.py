#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 11:10:29 2021

@author: tnye
"""

###############################################################################
# Script that makes the subplots for the datamap figure.
    # Figure with rhyp vs magnitude
    # Histogram of depths
###############################################################################

# Imports
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

model_name = 'model1'

df = pd.read_csv(f'/Users/tnye/kappa/data/flatfiles/final_dataset.csv')

# idx = np.where((np.array(df['Mag']) >= 3.5) & (np.array(df['Mag']) <= 5.5) &
#                (np.array(df['Rhyp']) <= 250))[0]
# mag = df['Mag'][idx]
# rhyp = df['Rhyp'][idx]
# depth = df['Depth'][idx]

mag = np.array(df['Mag'])
rhyp = np.array(df['Rhyp'])
depth = np.array(df['Depth'])

################################ Set up figure ################################

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'
fig = plt.figure(figsize=(8,6))

# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.02
hist_width = 0.2

# Set up axes
scatter = [left, bottom, width, height]
rhyp_hist = [left, bottom + height + spacing, width, hist_width]
mag_hist = [left + width + spacing, bottom, hist_width, width]
depth_hist = [left, bottom - height - spacing - hist_width, width, hist_width]

ax_plot = plt.axes(scatter)
ax_plot.tick_params(direction='in', left=False, bottom=False, labelsize=14)
ax_plot.scatter(rhyp, mag, s=5, alpha=0.5, c=depth, cmap='RdPu_r')
ax_plot.grid(linestyle='-',alpha=0.5)
ax_plot.set_axisbelow(True)
ax_plot.set_xlabel('Record R$_{hyp}$ (km)',fontsize=14)
ax_plot.set_ylabel('Magnitude',fontsize=14)
ax_plot.spines['bottom'].set_color('gray')
ax_plot.spines['top'].set_color('gray') 
ax_plot.spines['right'].set_color('gray')
ax_plot.spines['left'].set_color('gray')

ax_rhyp = plt.axes(rhyp_hist)
ax_rhyp.tick_params(left=False,right=False,bottom=False,top=False,labelleft=False,labelright=False,labelbottom=False,labeltop=False)
ax_rhyp.hist(rhyp, bins=np.arange(0,250,10), color='lightgray', edgecolor='darkgray', linewidth=.5)
ax_rhyp.grid(linestyle='-',alpha=0.5)
ax_rhyp.set_axisbelow(True)
ax_rhyp.spines['bottom'].set_color('gray')
ax_rhyp.spines['top'].set_color('gray') 
ax_rhyp.spines['right'].set_color('gray')
ax_rhyp.spines['left'].set_color('gray')

ax_mag = plt.axes(mag_hist)
ax_mag.tick_params(left=False,right=False,bottom=False,top=False,labelleft=False,labelright=False,labelbottom=False,labeltop=False)
ax_mag.hist(mag, bins=np.arange(3.5,5.75,0.25), orientation='horizontal', color='lightgray', edgecolor='darkgray', linewidth=.5)
ax_mag.grid(linestyle='-',alpha=0.5)
ax_mag.set_axisbelow(True)
ax_mag.spines['bottom'].set_color('gray')
ax_mag.spines['top'].set_color('gray') 
ax_mag.spines['right'].set_color('gray')
ax_mag.spines['left'].set_color('gray')
plt.show()

plt.savefig('/Users/tnye/kappa/plots/paper/data_stats.png', dpi=300)
plt.close()


################################ Set up figure ################################

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'
fig = plt.figure(figsize=(8,4))

# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.15, 0.65
spacing = 0.02
hist_width = 0.4

# Set up axes
depth_hist = [left, bottom, width, hist_width]
ax_depth = plt.axes(depth_hist)
ax_depth.tick_params(left=False,right=False,bottom=False,top=False,labelright=False,labeltop=False,labelsize=14)
bins_list = np.arange(0,21,1)
n, bins, patches = ax_depth.hist(depth, edgecolor='k', bins=bins_list, linewidth=.25)
for i in range(len(patches)):
    patches[i].set_facecolor(plt.cm.RdPu_r((bins_list[i]+bins_list[i+1])/2/20))
ax_depth.grid(linestyle='-',alpha=0.5)
ax_depth.set_axisbelow(True)
ax_depth.spines['bottom'].set_color('gray')
ax_depth.spines['top'].set_color('gray') 
ax_depth.spines['right'].set_color('gray')
ax_depth.spines['left'].set_color('gray')
ax_depth.set_xlabel('Depth (km)',fontsize=14)
ax_depth.set_ylabel('Counts',fontsize=14)
ax_depth.set_xlim(0,20)
ax_depth.tick_params(axis='y', labelrotation=45)
# plt.yticks(np.arange(0, 150, 50))
plt.show()

plt.savefig('/Users/tnye/kappa/plots/paper/depth.png', dpi=300)
plt.close()


