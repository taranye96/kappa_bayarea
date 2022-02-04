#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 12:04:17 2021

@author: tnye
"""

###############################################################################
# Script that makes text file with list of stations and how many events they 
# recorded to be plotted in GMT. 
###############################################################################

# Imports
import numpy as np
import pandas as pd
from glob import glob

model_name = 'model4'
home_dir = f'/Users/tnye/kappa/traditional_method/models/{model_name}'

outfilename = f'/Users/tnye/kappa/GMT/data/kappa/m4_density.txt'

files = sorted(glob(f'{home_dir}/stn_flatfiles/*.csv'))

log_stns = []
log_density = []
with open('/Users/tnye/kappa/traditional_method/logs/model4_event_density.log') as density_log:
    for line in density_log:
        log_stns.append(line.split(':')[0])
        log_density.append(int(line.split(' ')[1]))
    
        
lon = []
lat = []
num_events = []

for file in files:
    df = pd.read_csv(file,converters = {'Name':str})
    
    if len(df)>0:
        stn = df['Name'][0]
        slon = df['Slon'][0]
        slat = df['Slat'][0]
        lon.append(slon)
        lat.append(slat)
        
        stn_ind = np.where(np.array(log_stns)==stn)[0][0]
        num_events.append(log_density[stn_ind])

    
# Save kappa to text file
outfile = open(outfilename, 'w')
out = (np.array([lon, lat, num_events], dtype=object)).T

# Write value to file and save
outfile.write('#lon,lat,events\n')
np.savetxt(outfile, out, fmt='%1.3f, %1.3f, %d', delimiter=',')
outfile.close()