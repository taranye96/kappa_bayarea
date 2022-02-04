#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 22 11:47:28 2022

@author: tnye
"""

###############################################################################
# Script that makes text files for GMT with coordinates for final set of evetns
# and stations. 
###############################################################################

# Imports
import pandas as pd
import numpy as np

model_name = 'model1'

lst_str_cols = ['Station','Name']
dict_dtypes = {x : 'str'  for x in lst_str_cols}
final_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/final_dataset.csv', dtype=dict_dtypes)
catalog = pd.read_csv('/Users/tnye/kappa/data/flatfiles/full_bay_flatfile.csv', dtype=dict_dtypes)

events = np.unique(np.array(final_df['Event']))
stations = np.unique(np.array(final_df['Station']))

slon = []
slat = []
qlon = []
qlat = []
mag = []

for stn in stations:
    ind = np.where(np.array(catalog['Name'])==stn)[0][0]
    slon.append(catalog['Slon'].iloc[ind])
    slat.append(catalog['Slat'].iloc[ind])

for event in events:
    yyyy,mth,dd,hh,mm,sec = event.split('_')[1:]
    org = f'{yyyy}-{mth}-{dd} {hh}:{mm}:{sec}'
    ind = np.where(np.array(catalog['OrgT'])==org)[0][0]
    qlon.append(catalog['Qlon'].iloc[ind])
    qlat.append(catalog['Qlat'].iloc[ind])
    mag.append(catalog['Mag'].iloc[ind]/10)

# Make station file
outfile = open(f'/Users/tnye/kappa/GMT/data/sites/{model_name}_stns.txt', 'w')
out = np.column_stack((slon, slat))
outfile.write('#Longitude,Latitude\n')
np.savetxt(outfile, out, fmt='%s', delimiter=',')
outfile.close()

# Make station file
outfile = open(f'/Users/tnye/kappa/GMT/data/events/{model_name}_events.txt', 'w')
out = np.column_stack((qlon, qlat, mag))
outfile.write('#Longitude,Latitude,Scaled Magnitude\n')
np.savetxt(outfile, out, fmt='%s', delimiter=',')
outfile.close()