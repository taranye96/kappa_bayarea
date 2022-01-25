#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 21:40:55 2021

@author: tnye
"""

###############################################################################
# Script that makes input files for the USGS CENCALVM (central California
# velociy model) using all of the broadband and strong motion stations in my 
# kappa dataset. 
###############################################################################

# Imports
import numpy as np
import pandas as pd

stn_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/all_stations.csv')

for i in range(len(stn_df)):
    stn = stn_df['name'].iloc[i]
    lon = stn_df['longitude'].iloc[i]
    lat = stn_df['latitude'].iloc[i]
    
    log = f'/Users/tnye/kappa/gmm_parameters/etree/input/{stn}.txt'
    f = open(log,'w+')

    f.write(f'{lon}\t{lat}\t0\n')
    f.write(f'{lon}\t{lat}\t-100\n')
    f.write(f'{lon}\t{lat}\t-150\n')
    f.write(f'{lon}\t{lat}\t-200\n')
    f.write(f'{lon}\t{lat}\t-250\n')
    f.write(f'{lon}\t{lat}\t-300\n')
    f.write(f'{lon}\t{lat}\t-350\n')
    f.write(f'{lon}\t{lat}\t-400\n')
    f.write(f'{lon}\t{lat}\t-450\n')
    f.write(f'{lon}\t{lat}\t-500\n')
    f.write(f'{lon}\t{lat}\t-550\n')
    f.write(f'{lon}\t{lat}\t-600\n')
    f.write(f'{lon}\t{lat}\t-650\n')
    f.write(f'{lon}\t{lat}\t-700\n')
    f.write(f'{lon}\t{lat}\t-750\n')
    f.write(f'{lon}\t{lat}\t-800\n')
    f.write(f'{lon}\t{lat}\t-850\n')
    f.write(f'{lon}\t{lat}\t-900\n')
    f.write(f'{lon}\t{lat}\t-950\n')
    f.write(f'{lon}\t{lat}\t-1000\n')
    f.write(f'{lon}\t{lat}\t-1050\n')
    f.write(f'{lon}\t{lat}\t-1100\n')
    f.write(f'{lon}\t{lat}\t-1150\n')
    f.write(f'{lon}\t{lat}\t-1200\n')
    f.write(f'{lon}\t{lat}\t-1250\n')
    f.write(f'{lon}\t{lat}\t-1300\n')
    f.write(f'{lon}\t{lat}\t-1350\n')
    f.write(f'{lon}\t{lat}\t-1400\n')
    f.write(f'{lon}\t{lat}\t-1450\n')
    f.write(f'{lon}\t{lat}\t-1500')
    
    f.close()