#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 12:31:21 2022

@author: tnye
"""

# Imports
import numpy as np
import pandas as pd
from glob import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

stn = 'BJOB'
# stn = 'N001'
model_name = 'final_model'
home_dir = f'/Users/tnye/kappa/traditional_method/models/{model_name}'
spectra_dir = '/Users/tnye/kappa/data/spectra/acc'
min_events = 10
min_mag = 3.5
max_mag = 5.5
max_dist = 250

# Read in dataframe 
event_file = '/Users/tnye/kappa/data/flatfiles/full_bay_flatfile.csv'
lst_str_cols = ['Name','Channel']
dict_dtypes = {x : 'str'  for x in lst_str_cols}
event_df = pd.read_csv(event_file, dtype=dict_dtypes)

# Initial station dataframe
stn_df0 = event_df[(event_df['Name'] == stn)].reset_index(drop=True)
stn_df0 = stn_df0.sort_values(by=['OrgT'],ascending=False)

# Get list of event acceleration files for given station
channel = stn_df0['Channel'].iloc[0]
event_list = glob(f'{spectra_dir}/*/*_{stn}_*')
event_list = sorted(event_list,reverse=True)

# Make a list of the event origin times
record_origins = []
for i in range(len(event_list)):
    yyyy,mth,dd,hr,mm,sec = event_list[i].split('/')[-1].split('_')[4:10]
    sec = sec.split('.')[0]
    org = f'{yyyy}-{mth}-{dd} {hr}:{mm}:{sec}' 
    record_origins.append(org)

# Remove duplicate origin times for separate events
    # I'm assuming the 2nd is the one we have data for since it likely 
    # overwrote the first one
stn_df1 = stn_df0.copy().reset_index(drop=True)
ex_ind = [idx for idx, val in enumerate(np.array(stn_df0['OrgT'])) if val in np.array(stn_df0['OrgT'])[:idx]]
stn_df1 = stn_df1.drop(ex_ind)

# Remove events from event files that are not in the dataframe (perhaps strong motion but I wanted the bb data for that location)
missing = []
missing_ind = []
for i in range(len(record_origins)):
    if record_origins[i] not in np.array(stn_df1['OrgT']):
        missing.append(record_origins[i])
        missing_ind.append(i)
for index in sorted(missing_ind, reverse=True):
    del event_list[index]

# Remove events from df that are not in the event files
missing = []
missing_ind = []
stn_df2 = stn_df1.copy().reset_index()
for i in range(len(stn_df1)):
    if stn_df1['OrgT'].iloc[i] not in record_origins:
        a=stn_df1['OrgT'].iloc[i]
        missing.append(stn_df1['OrgT'].iloc[i])
        missing_ind.append(i)
        stn_df2 = stn_df2.drop([i])

   
# Remove events from df that do not fit the magnitude and dist criteria
ex = []
ex_ind = []
stn_df3 = stn_df2.copy()
for i in range(len(stn_df2)):
    rhyp = stn_df2['rhyp'].iloc[i]
    mag = stn_df2['Mag'].iloc[i]
    org = stn_df2['OrgT'].iloc[i]
    depth = stn_df2['Qdep'].iloc[i]/1000
    rrup = np.tan(np.arccos(depth/rhyp))*depth
    if rrup >= max_dist or mag <= min_mag or mag >= max_mag:
        ex_ind.append(i)
        ex.append(org)
        stn_df3 = stn_df3.drop(index=stn_df3.index[np.where(stn_df3['OrgT']==org)[0][0]])

# Remove events from event list that do not fit the magnitude and dist criteria
for index in sorted(ex_ind, reverse=True):
    del event_list[index]



################################# Make Figure #################################

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'

fig = plt.figure(figsize=(8,6))  

# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.075, 0.65

# Set up axes
rect_plot = [left, bottom, width, height]

ax_plot = plt.axes(rect_plot)
ax_plot.tick_params(direction='in', top=True, right=True)

# Set up colormap
viridis = plt.get_cmap('viridis_r') 
cNorm  = colors.Normalize(vmin=3.5, vmax=5.5)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=viridis)
    
# Loop through events
for i in range(len(stn_df3)):
    
    orgt = stn_df3['OrgT'].iloc[i]
    yyyy,mth,dd = orgt.split(' ')[0].split('-')
    hh,mm,sec = orgt.split(' ')[1].split(':')
    event = f'Event_{yyyy}_{mth}_{dd}_{hh}_{mm}_{sec}'
    mag = stn_df3['Mag'].iloc[i]
    
    files = glob(f'{spectra_dir}/{event}/*_{stn}_*')
    if len(files) > 0:
        event_file = glob(f'{spectra_dir}/{event}/*_{stn}_*')[0]
        # Read in data
        data = np.genfromtxt(event_file, comments = '#', dtype = float)
        freq = data.T[0]
        amp = data.T[1]
        
        colorVal = scalarMap.to_rgba(mag)
        ax_plot.semilogy(freq,amp,c=colorVal,alpha=0.5,lw=.8,label=org)

ax_plot.set_xlabel('Freq (Hz)', fontsize=12)
ax_plot.set_ylabel('Amp (m)', fontsize=12)
ax_plot.set_xlim(xmax=35)
plt.title(f'{stn}', fontsize=12)
cax = fig.add_axes([0.77, 0.075, 0.025, 0.65])
fig.colorbar(scalarMap, label='Magnitude', cax=cax)

plt.show()
plt.savefig(f'/Users/tnye/kappa/plots/paper/{stn}_spectra.png', dpi=300)
# plt.close()
