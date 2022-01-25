#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 11:56:54 2019

@author: eking
Revised by: Tara Nye

Cuts waveforms to 15 
"""

###############################################################################
# Script trims the waveforms to some time before and after S-arrival.
###############################################################################

# Standard Library Imports 
import numpy as np
import pandas as pd
from glob import glob
from os import path, makedirs
import datetime
# Local Imports
import kappa_utils as ku


################################ Parameters ###################################

# Working Directory 
working_dir = '/Users/tnye/kappa/data/waveforms'

# Main flatfile for events and stations
main_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/SNR_5_file.csv')

# Path to corrected seismograms
event_dirs = glob(working_dir + '/filtered/Event_*')

# Path to save cut waveforms
outpath = working_dir + '/cutWF_15'

# S-wave velocity (km/s)
Vs = 3.1 

# Number of seconds before and after the S-arrival to cut the record to
cutstart = 5
cutend = 15


############################### Cut waveforms #################################

# Get list of events 
events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))

# Make directories for cut waveforms for each event  
for i in range(len(events)):
    if not path.exists(outpath + '/' + events[i]):
        makedirs(outpath + '/'  + events[i])

# Loop through events 
for event in events:
    
    # Get list of records for this event
        # Just need one component because this is to get the list of stations
        # that recorded this event 
    recordpaths = glob(working_dir + '/filtered/' + event + '/*' + '_HHE*.sac')
    
    # Get event time
    yyyy, mth, dd, hr, mm, sec = event.split('_')[1:]
    event_time = f'{yyyy}_{mth}_{dd}_{hr}_{mm}_{sec}'
    
    # Get stations that recorded this event 
    stns = [(x.split('/')[-1]).split('_')[1] for x in recordpaths]
    
    # Loop through stations 
    for j in range(len(stns)):
        
        # Get North and East components 
        recordpath_N = glob(working_dir + '/filtered/' + event + '/*_' + stns[j] + '_HHN_' + event_time + '.sac')
        recordpath_E = glob(working_dir + '/filtered/' + event + '/*_' + stns[j] + '_HHE_' + event_time + '.sac')
       
        # Check that both a North and East component exist 
        if(len(recordpath_N) == 1 and len(recordpath_E) == 1):
            
            # Get hypocentral distance (km)
            stn_ind = np.where(main_df['Name']==stns[j])[0][0]
            rhyp = main_df['rhyp'].iloc[stn_ind]
            
            # Get origin time
            orig = datetime.datetime(int(yyyy),int(mth),int(dd),int(hr),int(mm),int(sec))
            
            # Calc S-wave arrival time in seconds after origin time 
            stime = rhyp/Vs
            
            # Calc S-wave arrival time as a datetime object
            Sarriv = orig + datetime.timedelta(0,stime)
            
            # Cut North component
            outpath_N = recordpath_N[0].replace('filtered','cutWF_' + np.str(cutend))
            ku.cut_swave(recordpath_N[0], outpath_N, Sarriv, cutstart, cutend)
             
            # Cut East component
            outpath_E = recordpath_E[0].replace('filtered','cutWF_' + np.str(cutend))
            ku.cut_swave(recordpath_E[0], outpath_E, Sarriv, cutstart, cutend)
    
    