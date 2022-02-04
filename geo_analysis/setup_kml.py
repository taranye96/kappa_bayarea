#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 13:42:48 2022

@author: tnye
"""

###############################################################################
#Script that makes a text file for Google Earth.
###############################################################################

import numpy as np
import pandas as pd

model_name = 'model1'
kappa_file = f'/Users/tnye/kappa/traditional_method/models/model1/{model_name}_kappa.out'

lst_str_cols = ['Name','Channel']
dict_dtypes = {x : 'str'  for x in lst_str_cols}
event_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/full_bay_flatfile.csv', dtype=dict_dtypes)

stns = []
kappa = []
lon = []
lat = []
with open(kappa_file) as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
        if i > 0:
            stn = line.split('\t')[0]
            stns.append(stn)
            kappa.append(float(line.split('\t')[1]))
            lon.append(event_df['Slon'].iloc[np.where(event_df['Name']==stn)[0][0]])
            lat.append(event_df['Slat'].iloc[np.where(event_df['Name']==stn)[0][0]])
    
# Save kappa to text file
outfilename = f'/Users/tnye/kappa/data/google_earth/{model_name}_kappa.txt'
outfile = open(outfilename, 'w')
out = (np.array([stns, lon, lat, kappa], dtype=object)).T

# Write value to file and save
outfile.write('#Station,Longitude,Latitude,Kappa(s)\n')
np.savetxt(outfile, out, fmt='%s, %1.3f, %1.3f, %s', delimiter=',')
outfile.close()