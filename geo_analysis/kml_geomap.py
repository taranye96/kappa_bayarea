#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 15:34:11 2021

@author: tnye
"""

###############################################################################
# Script used to determine what unit a station is in using a geologic map of 
# the Bay area. 
###############################################################################

import numpy as np
import pandas as pd
from glob import glob
import geopandas as gpd
from shapely.geometry import Point

kml_list = glob('/Users/tnye/kappa/data/google_earth/bay_kmls/*.kml')

# Get list of station coordinates
locations = pd.read_csv('/Users/tnye/kappa/GMT/data/kappa/updated_stns.txt')

kappa_file = pd.read_csv('/Users/tnye/kappa/traditional_method/models/final_model/final_model_kappa.out',delimiter='\t')

# Initialize dataframe
stn_df = pd.DataFrame(np.nan, np.arange(len(locations)), columns=['Station', 'Longitude', 'Latitude', 'Kappa', 'Geology', 'Unit', 'Description'])
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
    
    # Loop through kml files
    for kml in kml_list:
        gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
        my_map = gpd.read_file(kml, driver='KML')
        my_map['geometry'] = my_map['geometry'].buffer(0)
   
        # Check if point is in any of the polygons
        try:
            check = my_map.contains(point).values
            if True in check:
                idx = np.where(check==True)[0][0]
                unit = my_map.iloc[idx].Name
                description = my_map.iloc[idx].Description
               
                if unit[0] == 'a' or unit[0] == 'Q':
                    # if stn_df['Latitude'].iloc[i] >= 37.45 and stn_df['Longitude'] >= -122.38:
                    #     stn_df['Geology'].iloc[i] = 'East Quaternary'
                    #     stn_df['Geology'].iloc[i] = 'West Quaternary'
                    # elif stn_df['Latitude'].iloc[i] <= 37.45 and stn_df['Longitude'] <= -122.38:
                        if 'alluvium' in description.lower():
                            stn_df['Geology'].iloc[i] = 'Alluvium'
                        elif 'mud' in description.lower():
                            stn_df['Geology'].iloc[i] = 'Mud'
                        elif 'artificial' in description.lower():
                            stn_df['Geology'].iloc[i] = 'Artificial Fill'
                        elif 'sand' in description.lower():
                            stn_df['Geology'].iloc[i] = 'Sand'
                        else:
                            stn_df['Geology'].iloc[i] = 'Misc. Quaternary'
                elif 'sediment' in description.lower():
                    stn_df['Geology'].iloc[i] = 'Sedimentary'
                elif 'volcanic' in description.lower():
                    stn_df['Geology'].iloc[i] = 'Volcanic'
                elif 'meta' in description.lower() or 'ite' in description.lower():
                    stn_df['Geology'].iloc[i] = 'Metamorphic'
                else:
                    stn_df['Geology'].iloc[i] = 'unknown'
                
                stn_df['Unit'].iloc[i] = unit
                stn_df['Description'].iloc[i] = description
                
                break
        
            else:
                continue
        except:
            continue

stn_df.to_csv('/Users/tnye/kappa/data/stn_geology2.csv')
                

