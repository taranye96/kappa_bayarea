#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 11:45:00 2021

@author: tnye
"""

# Imports 
import numpy as np
import pandas as pd
import math
from pandas import DataFrame
from pyproj import Proj

lonmin = -123.5
lonmax = -120.55
latmin = 36.05
latmax = 39.34

x = np.arange(lonmin, lonmax, 0.01)
y = np.arange(latmin, latmax, 0.01)

# Make grid
lon = []
lat = []
for i in range(len(x)):
    for j in range(len(y)):
        lon.append(x[i])
        lat.append(y[j])

# Rotate 30 degrees counter clockwise
def rotate(p, origin=(0, 0), degrees=0):
    angle = np.deg2rad(degrees)
    R = np.array([[np.cos(angle), -np.sin(angle)],
                  [np.sin(angle),  np.cos(angle)]])
    o = np.atleast_2d(origin)
    p = np.atleast_2d(p)
    return np.squeeze((R @ (p.T-o.T) + o.T).T)

lonlat_center=[(lonmin+lonmax)/2, (latmin+latmax)/2]
# lonlat_center=[-122, 37.5]

rot_lon = []
rot_lat = []
for i in range(len(x)):
    for j in range(len(y)):
    # lonr, latr = rotate([lon[i],lat[i]], lonlat_center, np.deg2rad(30))
        lonr = (x[i]-lonlat_center[0])*np.cos(np.deg2rad(30)) - (y[j]-lonlat_center[1])*np.sin(np.deg2rad(30)) + lonlat_center[0]
        latr = (x[i]-lonlat_center[0])*np.sin(np.deg2rad(30)) + (y[j]-lonlat_center[1])*np.cos(np.deg2rad(30)) + lonlat_center[1]
        rot_lon.append(lonr)
        rot_lat.append(latr)

f = open('/Users/tnye/kappa/krige/grid_lonlat.txt','w')
f.write('Longitude,Latitude\n')
for i in range(len(lon)):
    f.write(f"{lon[i]},{lat[i]}\n")
f.close()

f = open('/Users/tnye/kappa/krige/rotgrid_lonlat.txt','w')
f.write('Longitude,Latitude\n')
for i in range(len(rot_lon)):
    f.write(f"{rot_lon[i]},{rot_lat[i]}\n")
f.close()

