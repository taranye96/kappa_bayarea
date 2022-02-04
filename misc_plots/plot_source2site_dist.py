#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 10:39:25 2021

@author: tnye
"""

###############################################################################
# Script that makes a histogram of source to station distances. 
###############################################################################

# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Flatfile with all of the event and station metadata
event_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/SNR_5_file.csv')

stations = np.unique(event_df['Name'])

med_dist = []
avg_dist = []

for stn in stations:
    rhyp_list = []
    
    # Get indices for where station occurs in dataframe
    ind = np.where(event_df['Name']==stn)
    
    # Loop through indexes for this station and get rhyp
    for i in ind:
        rhyp_list.append(event_df['rhyp'][i])
    
    med_dist.append(np.median(rhyp_list))
    avg_dist.append(np.mean(rhyp_list))

# Get station coordinates
lon_list = np.array([])
lat_list = np.array([])

for stn in stations:
    
    # Fist occurence of station in catalog
    ind = np.where(event_df==stn)[0][0]
   
    # Get station coordiantes
    lon = round(event_df['Slon'][ind],3)
    lat = round(event_df['Slat'][ind],3)

    # Append lon and lat to lists
    lon_list = np.append(lon_list,lon)
    lat_list = np.append(lat_list,lat)

# Combine coordinates and distances into 1 array
stn_med = np.array([lon_list, lat_list, med_dist])
stn_med = stn_med.T
stn_avg = np.array([lon_list, lat_list, avg_dist])
stn_avg = stn_avg.T

# Save median distanced information 
path = '/Users/tnye/kappa/data/bay_med_dist.txt'
with open(path, 'w+') as datafile_id:
    np.savetxt(datafile_id, np.array(['Longtidude', 'Latitude', 'MedianDist(km)'])[np.newaxis],
               fmt='%s',delimiter=',')
    np.savetxt(datafile_id, stn_med, fmt='%1.3f',delimiter=',')
    
# Save average distance information 
path = '/Users/tnye/kappa/data/bay_avg_dist.txt'
with open(path, 'w+') as datafile_id:
    np.savetxt(datafile_id, np.array(['MedianDist(km)','AverageDist(km)'])[np.newaxis],
               fmt='%s',delimiter=',')
    np.savetxt(datafile_id, stn_avg, fmt='%1.3f',delimiter=',')
    

# Plot histogram of median distances
plt.hist(med_dist,12)
plt.xlabel('Median Distance from Source (km)')
plt.savefig('/Users/tnye/kappa/plots/distance/median_dist.png', dpi=300)
plt.close()

# Plot histogram of average distances
plt.hist(avg_dist,12)
plt.xlabel('Average Distance from Source (km)')
plt.savefig('/Users/tnye/kappa/plots/distance/average_dist.png', dpi=300)
plt.close()