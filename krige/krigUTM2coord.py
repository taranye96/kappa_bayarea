#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 11:59:05 2021

@author: tnye
"""

# Imports
import numpy as np
import pandas as pd
from pandas import DataFrame
from pyproj import Proj

###################################### UTM ###################################

# df = pd.read_csv('/Users/tnye/kappa/krige/krige_k0.csv')
# UTMx = df['UTMx']
# UTMy = df['UTMy']
# kappa = df['pred_k0']

# myProj = Proj("+proj=utm +zone=10, +north +ellps=WGS84 +datum=WGS84 +units=km +no_defs")

# lon, lat = myProj(UTMx,UTMy, inverse=True)

# # Save kappa to text file
# # outfile = open('/Users/tnye/kappa/GMT/data/kappa/pred_k0_lonlat.txt', 'w')
# outfile = open('/Users/tnye/kappa/krige/pred_k0_lonlat.txt', 'w')
# out = (np.array([lon, lat, kappa], dtype=object)).T

# # Write value to file and save
# outfile.write('#Longitude,Latitude,Kappa(s)\n')
# np.savetxt(outfile, out, fmt='%1.3f, %1.3f, %s', delimiter=',')
# outfile.close()

# # Save kappa to text file
# # outfile = open('/Users/tnye/kappa/GMT/data/kappa/pred_k0_utm.txt', 'w')
# outfile = open('/Users/tnye/kappa/krige/pred_k0_utm.txt', 'w')
# out = (np.array([UTMx, UTMy, kappa], dtype=object)).T

# # Write value to file and save
# outfile.write('#UTMx,UTMy,Kappa(s)\n')
# np.savetxt(outfile, out, fmt='%1.3f, %1.3f, %s', delimiter=',')
# outfile.close()


############################## Longitude Latitude ############################

df = pd.read_csv('/Users/tnye/kappa/krige/krige_k0_lonlat.csv')

lon = df['Longitude']
lat = df['Latitude']
kappa = df['pred_k0']
k0_stddev = df['k0_stddev']

# Save kappa to text file
outfile = open('/Users/tnye/kappa/krige/pred_k0_lonlat.txt', 'w')
out = (np.array([lon, lat, kappa], dtype=object)).T

# Write value to file and save
outfile.write('#Longitude,Latitude,Kappa(s)\n')
np.savetxt(outfile, out, fmt='%1.3f, %1.3f, %1.4f', delimiter=',')
outfile.close()

# Save standard deviation to text file
outfile = open('/Users/tnye/kappa/krige/pred_stddev_lonlat.txt', 'w')
out = (np.array([lon, lat, k0_stddev], dtype=object)).T

# Write value to file and save
outfile.write('#Longitude,Latitude,Kappa(s)\n')
np.savetxt(outfile, out, fmt='%1.3f, %1.3f, %1.4f', delimiter=',')
outfile.close()


