#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 22:29:30 2021

@author: tnye
"""

###############################################################################
# Script that sets up an initial dataframe with parameters needed to call the 
# ASK14 and BA14 Ground-Motion Models for the San Francisco Bay Area dataset. 
###############################################################################

# Imports
import numpy as np
import pandas as pd
from os import path
from glob import glob

model_name = 'model4.6.7'

lst_str_cols = ['Name','Channel']
dict_dtypes = {x : 'str'  for x in lst_str_cols}
df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/data_download/bay_2000-2021_SNR_irfilt.csv', dtype=dict_dtypes)

wf_dir = '/Users/tnye/kappa/data/waveforms/bay_acc_2000-2021/instrument_corrected/cut_SNR_3_80%'

with open(f'/Users/tnye/kappa/data/Vs30/{model_name}_stations_culled.txt') as file:
    lines = file.readlines()
    culled_stns = np.array([line.rstrip() for line in lines])
all_stns = pd.read_csv(f'/Users/tnye/kappa/data/Vs30/{model_name}_stations_culled.txt')

vs30_file = np.genfromtxt(f'/Users/tnye/kappa/data/Vs30/{model_name}_vs30.txt')
meas_vs30_df = pd.read_csv('/Users/tnye/kappa/data/Vs30/VS30_mcphillips_2020.csv')
meas_vs30_list = np.array(meas_vs30_df['VS30__M_S_'])
meas_vs30_stns = np.array(meas_vs30_df['NETWORK_ST'])
for i in range(len(meas_vs30_stns)):
    try:
        meas_vs30_stns[i] = meas_vs30_stns[i].split('.')[1]
    except:
        continue

nga_events_df = pd.read_csv('/USers/tnye/kappa/data/flatfiles/NGA_West2_Finite_Fault_Info_050142013.csv')
nga_stns_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/NGA_West2_SiteDatabase_V032.csv')
focal_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/event_focals.csv')

# Get part of dataframe with time domain SNR>5
SNR_ind = np.where((df['SNR_E_time']>=3) & (df['SNR_N_time']>=3))[0]
df_SNR = df.loc[SNR_ind]

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
vs30_meas = []
vs30 = []
ztor = []
rakes = []
dips = []
widths = []
z1pt0 = []

for i in range(len(df_SNR)):
    
    stn = df_SNR['Name'].iloc[i]
    if stn in culled_stns:
        org = df_SNR['OrgT'].iloc[i]
        yyyy,mth,dd = org.split(' ')[0].split('-')
        hh,mm,ss = org.split(' ')[1].split(':')
        event = f'Event_{yyyy}_{mth}_{dd}_{hh}_{mm}_{ss}'
        
        # Add to dataframe if waveforms exist 
        if len(glob(f'{wf_dir}/{event}/*_{stn}_H*')) > 0:
        
            # Get basic event and site info
            stn = df_SNR['Name'].iloc[i]
            quake = df_SNR['Quake#'].iloc[i]
            names.append(stn)
            events.append(quake)
            Qlon.append(df_SNR['Qlon'].iloc[i])
            Qlat.append(df_SNR['Qlat'].iloc[i])
            Qdep.append(df_SNR['Qdep'].iloc[i])
            origins.append(df_SNR['OrgT'].iloc[i])
            mags.append(df_SNR['Mag'].iloc[i])
            dists.append(df_SNR['rhyp'].iloc[i])
            Slon.append(df_SNR['Slon'].iloc[i])
            Slat.append(df_SNR['Slat'].iloc[i])
            elev.append(df_SNR['Selv'].iloc[i])
            ztor.append(df_SNR['Qdep'].iloc[i]/1000) #in km
            
            qlatr = round(df_SNR['Qlat'].iloc[i],5)
            qlonr = round(df_SNR['Qlon'].iloc[i],5)
            qlat = df_SNR['Qlat'].iloc[i]
            qlon = df_SNR['Qlon'].iloc[i]
            
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
                mag = df_SNR['Mag'].iloc[i]
                if mag < 5:
                    # Wells and Coppersmith 1984 scaling law
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
            # Check if measured values exist:
            if stn != 'GDXB':
                if stn in meas_vs30_stns:
                    stn_idx = np.where(meas_vs30_stns==stn)[0][0]
                    v = meas_vs30_list[stn_idx]
                    vs30_meas.append(v)
                    vs30.append(v)
                
                # Check if in NGA flatfile
                elif stn in np.array(nga_stns_df['Station ID']):
                    ind = np.where(nga_stns_df['Station ID']==stn)[0][0]
                    v = nga_stns_df['Vs30 for Analysis (m/s)'].iloc[ind]
                    vs30.append(v)
                    if np.isnan(float(nga_stns_df['Measured  Vs30 (m/s) when zp > 30 m; inferred from Vsz otherwise '].iloc[ind])):
                        vs30_meas.append(0.0)
                    else:
                        vs30_meas.append(v)
                    
                # Else get from global grid file
                else:
                    vs30.append(vs30_file[np.where(all_stns==stn)[0][0]])
                    vs30_meas.append(0.0)
            else:
                # Check if in NGA flatfile
                if stn in np.array(nga_stns_df['Station ID']):
                    ind = np.where(nga_stns_df['Station ID']==stn)[0][0]
                    v = nga_stns_df['Vs30 for Analysis (m/s)'].iloc[ind]
                    vs30.append(v)
                    if np.isnan(float(nga_stns_df['Measured  Vs30 (m/s) when zp > 30 m; inferred from Vsz otherwise '].iloc[ind])):
                        vs30_meas.append(0.0)
                    else:
                        vs30_meas.append(v)
                    
                # Else get from global grid file
                else:
                    vs30.append(vs30_file[np.where(all_stns==stn)[0][0]])
                    vs30_meas.append(0.0)
            
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
            # z1ref = (1/1000)*np.exp((-7.67/4)*np.log()) #or you can use reference Z1.0


# Crete dictionary 
dataset_dict = {'Station Name':names,'Station Longitude':Slon,'Station Latitude':Slat,
                'Station Elevation':elev,'Event ID':events,'Longitude':Qlon,'Latitude':Qlat,
                'Depth':Qdep,'OrgT':origins,'Magnitude':mags,'Ztor(km)':ztor,
                'Rake':rake,'Dip':dip,'Width(km)':width,'Rrup(km)':dists,'Vs30_meas(m/s)':vs30_meas,
                'Vs30(m/s)':vs30,'Z1pt0(m)':z1pt0}
gmm_df = pd.DataFrame(data=dataset_dict)
gmm_df.to_csv(f'/Users/tnye/kappa/data/flatfiles/gmm_init_{model_name}_culled.csv')
