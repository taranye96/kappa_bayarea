#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 09:38:37 2021

@author: tnye
"""

# Imports
import numpy as np
from pandas import DataFrame
from pyproj import Proj

file = '/Users/tnye/kappa/GMT/data/kappa/rrup.txt'     
lon = np.genfromtxt(file,skip_header=1,dtype=None,delimiter=', ')[:,0]
lat = np.genfromtxt(file,skip_header=1,dtype=None,delimiter=', ')[:,1]
kappa = np.genfromtxt(file,skip_header=1,dtype=None,delimiter=', ')[:,2]
stns = []

stn_file = '/Users/tnye/kappa/traditional_method/models/test/rrup/rrup_kappa.out'

with open(stn_file) as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
        if i > 0:
            stns.append(line.split('\t')[0])

# stn_array = np.array([])
# for i in range(len(stns)):
#     stn_array
        
logkappa = []
for i in range(len(kappa)):
    logkappa.append(np.log10(kappa[i]))

myProj = Proj("+proj=utm +zone=10, +north +ellps=WGS84 +datum=WGS84 +units=km +no_defs")

UTMx, UTMy = myProj(lon,lat, inverse=False)

# df = DataFrame(np.c_[UTMx, UTMy, lon, lat], columns=['UTMx', 'UTMy', 'Lon', 'Lat'])
# df = DataFrame(np.c_[stns, UTMx, UTMy, kappa], columns=['Station', 'UTMx', 'UTMy', 'log10_k0'])
df = DataFrame(np.c_[stns, UTMx, UTMy, lon, lat, logkappa, kappa], columns=['Station', 'UTMx', 'UTMy', 'Longitude', 'Latitude', 'log10_k0', 'kappa'])
df.to_csv('/Users/tnye/kappa/krige/coords.csv')

f = open('/Users/tnye/kappa/krige/ne.txt','w')
f.write('Northing\tEasting\n')
for i in range(len(UTMx)):
    f.write(f"{UTMy[i]*1000}\t{UTMx[i]*1000}\n")
f.close()

f = open('/Users/tnye/kappa/krige/lonlat.txt','w')
f.write('Latitude\tLongitude\n')
for i in range(len(lat)):
    f.write(f"{lat[i]}\t{lon[i]}\n")
f.close()