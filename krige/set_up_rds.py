#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 12 06:06:13 2021

@author: tnye
"""

###############################################################################
# Script that makes a coords csv that is used to make the .rds file in R for 
# the residual decomposition. 
###############################################################################

# Imports 
import pandas as pd
import numpy as np
from pyproj import Proj

model_name = 'model2'
home_dir = f'/Users/tnye/kappa/traditional_method/models/{model_name}'

kappa_df = pd.read_csv(f'{home_dir}/{model_name}_kappa.out',delimiter='\t')
coords_df = pd.read_csv(f'/Users/tnye/kappa/GMT/data/kappa/{model_name}.txt')

stns = kappa_df['#site ']
kappa = kappa_df[' kappa(s) ']
logkappa = []
for k in kappa:
    logkappa.append(np.log10(k))

lon = np.array(coords_df['#Longitude'])
lat = np.array(coords_df['Latitude'])

myProj = Proj("+proj=utm +zone=10, +north +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
UTMx,UTMy = myProj(lon,lat, inverse=False)
  
# Creates DataFrame
data = {'Station':stns, 'UTMx':UTMx, 'UTMy':UTMy, 'Longitude':lon, 'Latitude':lat, 'log10_kappa':logkappa, 'kappa':kappa}
df = pd.DataFrame(data)

df.to_csv(f'/Users/tnye/kappa/krige/{model_name}_coords.csv')
