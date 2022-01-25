#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 15:41:43 2021

@author: tnye
"""

###############################################################################
# Script that makes text filke with number of events each station recorded. 
###############################################################################

import numpy as np
import pandas as pd

cat = pd.read_csv('/Users/tnye/kappa/data/flatfiles/SNR_5_file.csv')

stns = np.unique(np.array(cat['Name']))

stn_num = []
stn_lon = []
stn_lat = []

for stn in stns:
    num = len(np.where(cat['Name']==stn)[0])
    stn_num.append(num)
    
    lon = np.array(cat['Slon'])[np.where(cat['Name']==stn)[0][0]]
    stn_lon.append(lon)
    
    lat = np.array(cat['Slat'])[np.where(cat['Name']==stn)[0][0]]
    stn_lat.append(lat)
    
f = open('/Users/tnye/kappa/data/stn_density.txt',"w")
f.write('Longitude,Latitude,Number of Events Recorded\n')
for i in range(len(stns)):
    f.write(f'{stn_lon[i]},{stn_lat[i]},{stn_num[i]}\n')
f.close()