#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 10:24:49 2021

@author: tnye
"""

###############################################################################
# Creates text files for Bay Area broadband station coordinates to be used in
# GMT figure.
###############################################################################

# Imports 
import numpy as np
import pandas as pd

catalog = pd.read_csv('/Users/tnye/kappa/data/flatfiles/SNR_5_file.csv')

stations = np.unique(np.array(catalog['Name']))

lon_list = np.array([])
lat_list = np.array([])

for stn in stations:
    
    # Fist occurence of station in catalog
    ind = np.where(catalog==stn)[0][0]
   
    # Get station coordiantes
    lon = round(catalog['Slon'][ind],3)
    lat = round(catalog['Slat'][ind],3)

    # Append lon and lat to lists
    lon_list = np.append(lon_list,lon)
    lat_list = np.append(lat_list,lat)


######################### Station Data Text File ##############################

# Combine coordinates into 1 array
stn_info = np.array([stations, lon_list, lat_list])
stn_info = stn_info.T

# Save coordinates
path = '/Users/tnye/kappa/data/stations.txt'
with open(path, 'w+') as datafile_id:
    np.savetxt(datafile_id, np.array(['Station','Longitude','Latitude'])[np.newaxis],
               fmt='%s',delimiter=',')
    np.savetxt(datafile_id, stn_info, fmt=['%s','%1.3f','%1.3f'],delimiter=',')
    

######################## Station Coords Text File #############################

# Combine coordinates into 1 array
coords = np.array([lon_list, lat_list])
coords = coords.T

# Save coordinates
path = '/Users/tnye/kappa/code/GMT/stn_coords.txt'
with open(path, 'w+') as datafile_id:
    np.savetxt(datafile_id, np.array(['Longitude','Latitude'])[np.newaxis],
               fmt='%s',delimiter=',')
    np.savetxt(datafile_id, coords, fmt='%1.3f',delimiter=',')


######################## Station Labels Text File #############################

# Combine coordinates into 1 array
labels = np.array([lon_list, lat_list, stations])
labels = labels.T

# Save coordinates
path = '/Users/tnye/kappa/code/GMT/stn_labels.txt'
with open(path, 'w+') as datafile_id:
    np.savetxt(datafile_id, np.array(['Longitude','Latitude','Station'])[np.newaxis],
               fmt='%s',delimiter=',')
    np.savetxt(datafile_id, labels, fmt=['%1.3f','%1.3f','%s'],delimiter=',')