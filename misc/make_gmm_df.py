#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 22:29:30 2021

@author: tnye
"""

###############################################################################
# Script that sets up an initial dataframe with parameters needed to call the 
# ASK14 Ground-Motion Model for the San Francisco Bay Area dataset. 
###############################################################################

# Imports
import numpy as np
import pandas as pd
from os import path
from glob import glob

model_name = 'model1'

lst_str_cols = ['Name','Channel']
dict_dtypes = {x : 'str'  for x in lst_str_cols}
df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/full_bay_flatfile.csv', dtype=dict_dtypes)

with open(f'/Users/tnye/kappa/data/Vs30/{model_name}_stations.txt') as file:
    lines = file.readlines()
    stns = np.array([line.rstrip() for line in lines])


vs30_file = np.genfromtxt(f'/Users/tnye/kappa/data/Vs30/{model_name}_vs30.txt')
nga_events_df = pd.read_csv('/USers/tnye/kappa/data/flatfiles/NGA_West2_Finite_Fault_Info_050142013.csv')
nga_stns_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/NGA_West2_SiteDatabase_V032.csv')
focal_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/event_focals.csv')

names = []
events = []
origins = []
Qlon = []
Qlat = []
Qdep = []
mags = []
dists = []
Slon = []
Slat = []
elev = []
vs30 = []
ztor = []
rakes = []
dips = []
widths = []
z1pt0 = []

for i in range(len(df)):
    
    stn = df['Name'].iloc[i]
    if stn in stns:
        org = df['OrgT'].iloc[i]
        yyyy,mth,dd = org.split(' ')[0].split('-')
        hh,mm,ss = org.split(' ')[1].split(':')
        event = f'Event_{yyyy}_{mth}_{dd}_{hh}_{mm}_{ss}'
        
        # Add to dataframe if waveforms exist 
        if len(glob(f'/Users/tnye/kappa/data/waveforms/acc_all/Vs3.1/cut2_20/{event}/*_{stn}_H*')) > 0:
        
            # Only use events M3.5+
            mag = df['Mag'].iloc[i]
            if mag >= 3.5:
        
                # Get basic event and site info
                stn = df['Name'].iloc[i]
                quake = df['Quake#'].iloc[i]
                names.append(stn)
                events.append(quake)
                Qlon.append(df['Qlon'].iloc[i])
                Qlat.append(df['Qlat'].iloc[i])
                Qdep.append(df['Qdep'].iloc[i])
                origins.append(df['OrgT'].iloc[i])
                mags.append(df['Mag'].iloc[i])
                dists.append(df['rhyp'].iloc[i])
                Slon.append(df['Slon'].iloc[i])
                Slat.append(df['Slat'].iloc[i])
                elev.append(df['Selv'].iloc[i])
                ztor.append(df['Qdep'].iloc[i]/1000) #in km
                
                qlatr = round(df['Qlat'].iloc[i],5)
                qlonr = round(df['Qlon'].iloc[i],5)
                qlat = df['Qlat'].iloc[i]
                qlon = df['Qlon'].iloc[i]
                
                # Rake, dip, width
                if quake[2:] in np.array(nga_events_df['Earthquake Name']):
                    ind = np.where(nga_events_df['Earthquake Name']==quake[2:])[0][0]
                    rake = nga_events_df['Rake (deg)'].iloc[ind]
                    dip = nga_events_df['Dip  (deg)'].iloc[ind]
                    width = nga_events_df['Width (km)'].iloc[ind]
                elif quake in np.array(focal_df['id']):
                    ind = np.where(focal_df['id']==quake)[0][0]
                    rake = focal_df['nc_np1_rake'].iloc[ind]
                    dip = focal_df['nc_np1_dip'].iloc[ind]
                    w = 'calc'
                else:
                    rake = 0
                    dip = 90
                    w = 'calc'
                # Get rid of faulty data (some had rake, dip = -999)
                if np.abs(rake)>180:
                        rake = 0
                        dip = 90
                
                if w == 'calc':
                    if mag < 5:
                        stress = 5*10**6
                        M0 = 10**((3/2)*mag + 9.1)
                        width = np.cbrt((7/16)*M0/stress)/1000
                    else:
                        # Thingbaijam scaling law for reverse earthquakes
                        b = 0.435
                        a = -1.669
                        width = 10**(a + (b*mag))/1000
            
                widths.append(width)
                rakes.append(rake)
                dips.append(dip)
                # Vs30
                if stn in np.array(nga_stns_df['Station ID']):
                    ind = np.where(nga_stns_df['Station ID']==stn)[0][0]
                    vs30.append(nga_stns_df['Vs30 for Analysis (m/s)'].iloc[ind])
                else:
                    vs30.append(vs30_file[np.where(stns==stn)[0][0]])
                
                # Z1pt0
                z_path = f'/Users/tnye/kappa/gmm_parameters/etree/output/{stn}.txt'
                fields = ['latitude', 'longitude', 'elevation', 'Vp', 'Vs', 'density', 'Qp',
                          'Qs', 'depth', 'fault id', 'zone id', 'elev of cell']
                if path.exists(z_path):
                    z_file = np.genfromtxt(z_path)
                    ind = min(range(len(z_file[:,4])), key=lambda i: abs(z_file[:,4][i]-1000))
                    z1pt0.append(np.abs(z_file[ind,2]))
                else:
                    z1pt0.append(np.nan)


# Crete dictionary 
dataset_dict = {'Station Name':names,'Station Longitude':Slon,'Station Latitude':Slat,
                'Station Elevation':elev,'Event ID':events,'Longitude':Qlon,'Latitude':Qlat,
                'Depth':Qdep,'OrgT':origins,'Magnitude':mags,'Ztor(km)':ztor,
                'Rake':rake,'Dip':dip,'Width(km)':width,'Rrup(km)':dists,'Vs30(m/s)':vs30,
                'Z1pt0(m)':z1pt0}
gmm_df = pd.DataFrame(data=dataset_dict)
gmm_df.to_csv(f'/Users/tnye/kappa/data/flatfiles/gmm_init_{model_name}.csv')
