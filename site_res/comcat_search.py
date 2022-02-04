#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 11:31:07 2021

@author: tnye
"""

###############################################################################
# Script used to find earthquake source info using comcat. 
###############################################################################

# Imports 
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
from libcomcat.search import count, get_event_by_id, search
from libcomcat.dataframes import get_detail_data_frame

events = search(starttime=datetime(2000, 1, 4),endtime=datetime(2018, 4, 30),minlatitude=36,maxlatitude=39,minlongitude=-124,
       maxlongitude=-121,maxmagnitude=6.1,minmagnitude=2.5)

df = get_detail_data_frame(events, get_all_magnitudes=True)

df.to_csv('/Users/tnye/kappa/data/flatfiles/event_focals.csv')

########
df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/full_bay_flatfile.csv')
focal_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/event_focals.csv')
nga_events = pd.read_csv('/Users/tnye/kappa/data/flatfiles/NGA_West2_Finite_Fault_Info_050142013.csv')
nga_sites = pd.read_csv('/Users/tnye/kappa/data/flatfiles/NGA_West2_SiteDatabase_V032.csv')

events = np.unique(np.array(df['Quake#']))

missed_events = []
missed_mag = []
for event in events:
    if event not in np.array(focal_df['id']):
        ind = np.where(df['Quake#']==event)[0][0]
        missed_events.append(event)
        missed_mag.append(df['Mag'].iloc[ind])
