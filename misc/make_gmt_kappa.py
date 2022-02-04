#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 12:00:49 2021

@author: tnye
"""

###############################################################################
# Module with function to make station GMT file using results from the 
# traditional (Anderson and Hough) approach. 
###############################################################################

# Imports
import numpy as np
import pandas as pd
from glob import glob


def make_gmt(home_dir, model_name, kappa_file, event_df):
    stns = []
    kappa = []
    lon = []
    lat = []
    logstd = []
    std = []
    with open(kappa_file) as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if i > 0:
                stn = line.split('\t')[0]
                stns.append(stn)
                kappa.append(float(line.split('\t')[1]))
                logstd.append(np.log10(float(line.split('\t')[2])))
                std.append(float(line.split('\t')[2]))
                lon.append(event_df['Slon'].iloc[np.where(event_df['Name']==stn)[0][0]])
                lat.append(event_df['Slat'].iloc[np.where(event_df['Name']==stn)[0][0]])
        
    # Save kappa to text file
    outfilename = f'/Users/tnye/kappa/GMT/data/kappa/{model_name}.txt'
    outfile = open(outfilename, 'w')
    out = (np.array([lon, lat, kappa], dtype=object)).T
    
    # Write value to file and save
    outfile.write('#Longitude,Latitude,Kappa(s)\n')
    np.savetxt(outfile, out, fmt='%1.3f, %1.3f, %s', delimiter=',')
    outfile.close()
    
    # Standard deviation to text file
    outfilename = f'/Users/tnye/kappa/GMT/data/kappa/{model_name}_std.txt'
    outfile = open(outfilename, 'w')
    out = (np.array([lon, lat, std], dtype=object)).T
    
    # Write value to file and save
    outfile.write('#Longitude,Latitude,Kappa Std Dev\n')
    np.savetxt(outfile, out, fmt='%1.3f, %1.3f, %s', delimiter=',')
    outfile.close()
    
    # Log Standard deviation to text file
    outfilename = f'/Users/tnye/kappa/GMT/data/kappa/{model_name}_logstd.txt'
    outfile = open(outfilename, 'w')
    out = (np.array([lon, lat, logstd], dtype=object)).T
    
    # Write value to file and save
    outfile.write('#Longitude,Latitude,Kappa log10 Std Dev\n')
    np.savetxt(outfile, out, fmt='%1.3f, %1.3f, %s', delimiter=',')
    outfile.close()
    
    return()

