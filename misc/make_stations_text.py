#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 15:19:43 2021

@author: tnye
"""

###############################################################################
# Script that makes a text file with the Bay Area stations and station
# coordinates 
###############################################################################

# Imports
import numpy as np
import pandas as pd
from obspy import read
from glob import glob

model_name = 'model1'

# Read in datafile
main_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/SNR_5_file.csv')

# Get list of stations
stations = np.unique(np.array(main_df['Name']))

# Initialize text file
log = f'/Users/tnye/kappa/data/stations_{model_name}.txt'
f = open(log,'w+')
f.write('# Station, Samprate, Longitude, Latitude, Rhyp(km)\n')

# Loop through stations
for stn in stations:
    
    # Get station info
    stn_ind = np.where(main_df['Name']==stn)[0][0]
    network = np.array(main_df['Network'])[stn_ind]
    st_lon = round(np.array(main_df['Slon'])[stn_ind],3)
    st_lat = round(np.array(main_df['Slat'])[stn_ind],3)
    rhyp = round(np.array(main_df['rhyp'])[stn_ind],3)
    
    # Get samprate
    st = read(glob(f'/Users/tnye/kappa/data/waveforms/acc/filtered/*/{network}_{stn}*.sac')[0])
    samprate = st[0].stats.sampling_rate
    # Save to file
    f.write(f'{network} {stn}\t{samprate}\t{st_lon}\t{st_lat}\t{rhyp}\n')

f.close()