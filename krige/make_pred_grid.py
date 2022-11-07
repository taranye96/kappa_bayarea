#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 11:45:00 2021

@author: tnye
"""

###############################################################################
# Script that makes a text file of grid points for kriging of k0.
###############################################################################

# Imports 
import numpy as np

lonmin = -123.5
lonmax = -120.55
latmin = 36.0
latmax = 39.4

x = np.arange(lonmin, lonmax, 0.005)
y = np.arange(latmin, latmax, 0.005)

# Make grid
lon = []
lat = []
for i in range(len(x)):
    for j in range(len(y)):
        lon.append(x[i])
        lat.append(y[j])

f = open('/Users/tnye/kappa/krige/grid_lonlat2.txt','w')
f.write('Longitude,Latitude\n')
for i in range(len(lon)):
    f.write(f"{lon[i]},{lat[i]}\n")
f.close()