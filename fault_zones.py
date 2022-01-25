#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 16:24:46 2022

@author: tnye
"""

###############################################################################
# Script used to determine whether a station is located within or outside of 
# a major fault zone in the Bay area. 
###############################################################################

import numpy as np
import pandas as pd
from glob import glob
import geopandas as gpd
from shapely.geometry import Point, Polygon

kml = '/Users/tnye/kappa/data/google_earth/damage_zones.kml'

# Get list of station coordinates
locations = pd.read_csv('/Users/tnye/kappa/GMT/data/kappa/updated_stns.txt')

kappa_file = pd.read_csv('/Users/tnye/kappa/traditional_method/models/final_model/final_model_kappa.out',delimiter='\t')

# Initialize dataframe
stn_df = pd.DataFrame(np.nan, np.arange(len(locations)), columns=['Station', 'Longitude', 'Latitude', 'Kappa', 'Fault_Zone'])
stn_df['Longitude'] = locations['#Longitude']
stn_df['Latitude'] = locations['Latitude']
stn_df['Station'] = kappa_file['#site ']
stn_df['Kappa'] = kappa_file[' kappa(s) ']

# Turn coordinates into tuple array and GeoDataFrame
geometric_points = []
for xy in zip(locations['#Longitude'], locations['Latitude']):
   geometric_points.append(Point(xy))
geo_locations = gpd.GeoDataFrame(locations,
                                 crs = {'init': 'epsg:4326'},
                                 geometry = geometric_points)
# Loop through stations
for i in range(0, len(geo_locations)):
    point = geo_locations.geometry.loc[i]
    
    print(stn_df['Station'].iloc[i])
    
    gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
    my_map = gpd.read_file(kml, driver='KML')
    my_map['geometry'] = my_map['geometry'].buffer(0)
   
    # Check if point is in any of the polygons
    try:
        check = my_map.contains(point).values
        if True in check:
            stn_df['Fault_Zone'].iloc[i] = 'within'
    
        else:
            stn_df['Fault_Zone'].iloc[i] = 'outside'
    except:
        continue

stn_df.to_csv('/Users/tnye/kappa/data/fault_zone_stns.csv')

