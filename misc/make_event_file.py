#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 15:59:34 2021

@author: tnye
"""

###############################################################################
# Creates text file for Bay Area event coordinates to be used in GMT figure.
###############################################################################

# Imports 
import numpy as np
import pandas as pd

catalog = pd.read_csv('/Users/tnye/kappa/GMT/data/events/events_m3.5.csv')
lon = catalog['Longitude']
lat = catalog['Latitude']
mag = catalog['Magnitude']

for i in range(len(mag)):
    
    mag[i] = mag[i]/10
    

######################## Station Coords Text File #############################

# Combine coordinates into 1 array
coords = np.array([lon, lat, mag])
coords = coords.T

# Save coordinates
path = '/Users/tnye/kappa/GMT/data/events/events_m3.5.txt'
with open(path, 'w+') as datafile_id:
    np.savetxt(datafile_id, np.array(['Longitude','Latitude','Magnitude'])[np.newaxis],
               fmt='%s',delimiter=',')
    np.savetxt(datafile_id, coords, fmt='%1.3f',delimiter=',')
