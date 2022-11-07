#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 13:53:35 2021

@author: tnye
"""

###############################################################################
# Script that gets Vs30 values for the final set of stations using the global 
# Vs30 model. This script is called by get_Vs30.sh
###############################################################################

# Imports
import numpy as np
import pandas as pd
import sys

model_name = sys.argv[1]

kappa_file = f'/Users/tnye/kappa/traditional_method/models/{model_name}/{model_name}_kappa_culled.out'

lst_str_cols = ['Name','Channel']
dict_dtypes = {x : 'str'  for x in lst_str_cols}
event_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/data_download/bay_2000-2021_SNR_irfilt.csv', dtype=dict_dtypes)

stns = []
lon = []
lat = []
with open(kappa_file) as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
        if i > 0:
            stn = line.split(',')[0]
            stns.append(stn)
            lon.append(event_df['Slon'].iloc[np.where(event_df['Name']==stn)[0][0]])
            lat.append(event_df['Slat'].iloc[np.where(event_df['Name']==stn)[0][0]])
    
# Save Vs30 to text file
outfile = open(f'/Users/tnye/kappa/data/Vs30/{model_name}_stn_coords_culled.txt', 'w')
out = (np.array([lon, lat], dtype=object)).T
np.savetxt(outfile, out, fmt='%1.3f, %1.3f', delimiter=',')
outfile.close()

# Save Stations to text file
outfile = open(f'/Users/tnye/kappa/data/Vs30/{model_name}_stations_culled.txt', 'w')
out = (np.array([stns], dtype=object)).T
np.savetxt(outfile, out, fmt='%s', delimiter=',')
outfile.close()