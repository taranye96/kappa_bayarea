#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 11:35:52 2021

@author: tnye
"""

###############################################################################
# Creates text file for All Bay area event coordinates to be used in GMT figure.
###############################################################################

# Imports 
import numpy as np
import pandas as pd

catalog = pd.read_csv('/Users/tnye/kappa/data/flatfiles/bay_events.csv')

lon_list = np.array(catalog['longitude'])
lat_list = np.array(catalog['latitude'])
mag_list = np.array(catalog['mag'])

for i, mag in enumerate(mag_list):
    mag_list[i] = mag/10
    

######################## Station Coords Text File #############################

# Combine coordinates into 1 array
coords = np.array([lon_list, lat_list, mag_list])
coords = coords.T

# Save coordinates
path = '/Users/tnye/kappa/GMT/all_bay_events.txt'
with open(path, 'w+') as datafile_id:
    np.savetxt(datafile_id, np.array(['Longitude','Latitude','Magnitude'])[np.newaxis],
               fmt='%s',delimiter=',')
    np.savetxt(datafile_id, coords, fmt='%1.3f',delimiter=',')
