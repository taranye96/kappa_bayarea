#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 13:53:35 2021

@author: tnye
"""

###############################################################################
# Script that gets Vs30 values for the final set of stations using the global 
# Vs30 model. 
###############################################################################

# Imports
import numpy as np
import pandas as pd
from glob import glob
from scipy.optimize import curve_fit
import sys

model_name = sys.argv[1]

kappa_file = f'/Users/tnye/kappa/traditional_method/models/{model_name}/{model_name}_kappa.out'

lst_str_cols = ['Name','Channel']
dict_dtypes = {x : 'str'  for x in lst_str_cols}
event_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/full_bay_flatfile.csv', dtype=dict_dtypes)

stns = []
lon = []
lat = []
with open(kappa_file) as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
        if i > 0:
            stn = line.split('\t')[0]
            stns.append(stn)
            lon.append(event_df['Slon'].iloc[np.where(event_df['Name']==stn)[0][0]])
            lat.append(event_df['Slat'].iloc[np.where(event_df['Name']==stn)[0][0]])
    
# Save Vs30 to text file
outfile = open(f'/Users/tnye/kappa/data/Vs30/{model_name}_stn_coords.txt', 'w')
out = (np.array([lon, lat], dtype=object)).T
np.savetxt(outfile, out, fmt='%1.3f, %1.3f', delimiter=',')
outfile.close()

# Save Stations to text file
outfile = open(f'/Users/tnye/kappa/data/Vs30/{model_name}_stations.txt', 'w')
out = (np.array([stns], dtype=object)).T
np.savetxt(outfile, out, fmt='%s', delimiter=',')
outfile.close()