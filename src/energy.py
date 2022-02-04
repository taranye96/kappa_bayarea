#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 13:18:30 2021

@author: tnye
"""

###############################################################################
# Script used to calculate how long it takes to reach 80% of Arias intensity
# (representative of total energy). This was to help motivate the waveform
# trimming windows. 
###############################################################################

# Standard Library Imports 
import numpy as np
import pandas as pd
from glob import glob
from os import path
from obspy import read
import datetime
from obspy.core.utcdatetime import UTCDateTime

# Local Imports
import kappa_utils as ku
import signal_average_fns as avg
import arias_intensity as AI

################################ Parameters ###################################

# Working Directory 
working_dir = '/Users/tnye/kappa/data/waveforms/acc'

# Main flatfile for events and stations
main_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/SNR_5_file.csv')

# Path to corrected seismograms
event_dirs = glob(working_dir + '/filtered/Event_*')

# S-wave velocity (km/s)
Vs = 3.1 

# Get list of events 
events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))
    
log = '/Users/tnye/kappa/data/energy.log'
f = open(log,'w+')

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
        
        # Read in records
        N_record = read(recordpath_N[0])
        E_record = read(recordpath_E[0])
        
        N_data = N_record[0].data
        E_data = E_record[0].data
        
        # Make sure records have same length
        if len(N_data) != len(E_data):
                
                # Check for gaps:
                st = E_record.append(N_record[0])
                gaps = st.get_gaps()
                
                if gaps == []:
                    # Minimum record length between the components
                    npts = np.min([len(N_data),len(E_data)])
                    
                    # Shorten the components to be the same length 
                    N_data = N_data[:npts]
                    E_data = E_data[:npts]
        
        NE = avg.get_eucl_norm_2comp(N_data, E_data)
        
        # Calculate Arias intensity 
        dt = 1/(N_record[0].stats.sampling_rate)
        arias, Narias = AI.get_arias_intensity(NE, dt)
        
        # Calculate time from start of record to reach 80% of total Arias Intensity 
        time_to_80 = AI.get_time_from_percent(Narias, 0.8, dt)
        time80 = N_record[0].stats.starttime + datetime.timedelta(0,time_to_80)
        
        # Get hypocentral distance (km)
        stn_ind = np.where((main_df['Name']==stns[j]) & (main_df['OrgT']==f'{yyyy}-{mth}-{dd} {hr}:{mm}:{sec}'))
        rhyp = np.array(main_df['rhyp'])[stn_ind][0]
        
        # Get origin time
        orig = datetime.datetime(int(yyyy),int(mth),int(dd),int(hr),int(mm),int(sec))
        
        # Calc S-wave arrival time in seconds after origin time 
        stime = rhyp/Vs
        
        # Calc S-wave arrival time as a datetime object
        Sarriv = orig + datetime.timedelta(0,stime)
        
        # Calculate time from S-arrival to 80% energy
        time_after_S = time80 - UTCDateTime(Sarriv)
        
        # Save to log file
        f.write(f'{event},{stns[j]},{time_after_S}\n')

f.close()
        