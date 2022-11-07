#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 22 18:01:16 2022

@author: tnye
"""

###############################################################################
# Makes .csv file with station kappa estiamtes to be used to make a kml.
###############################################################################

# Imports 
import numpy as np
import pandas as pd

model_name = 'final_model3'

lst_str_cols = ['Name']
dict_dtypes = {x : 'str'  for x in lst_str_cols}

kappa_file = f'/Users/tnye/kappa/traditional_method/models/final_model3/{model_name}_kappa.out'
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
outfilename = f'/Users/tnye/kappa/data/google_earth/Bay_kappa.txt'
outfile = open(outfilename, 'w')
out = (np.array([stns, lon, lat, kappa], dtype=object)).T

# Write value to file and save
outfile.write('#Station,Longitude,Latitude,Kappa(s)\n')
np.savetxt(outfile, out, fmt='%s, %1.3f, %1.3f, %s', delimiter=',')
outfile.close()
