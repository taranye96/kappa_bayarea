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
import arias_intensity_fns as AI

################################ Parameters ###################################

# Working Directory 
working_dir = '/Users/tnye/kappa/data/waveforms/bay_acc_2000-2021/instrument_corrected/'

# Main flatfile for events and stations
main_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/data_download/bay_2000-2021_SNR_irfilt.csv')

wf_dir = 'SNR_3'

# Path to corrected seismograms
event_dirs = sorted(glob(f'{working_dir}/{wf_dir}/Event_*'))

# Get list of events 
events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))
    
log = '/Users/tnye/kappa/data/bay_2000-2021_energy_SNR_3_irfilt.log'
f = open(log,'w+')
f.write('Event\tStation\tSec after Sarr (N)\tSec after Sarr (E)\tSec after Sarr (Z)\n')

# Loop through events 
for event in events:
    
    # Get list of records for this event
        # Just need one component because this is to get the list of stations
        # that recorded this event 
    recordpaths = glob(f'{working_dir}/{wf_dir}/{event}/*_H*E*.mseed')
    
    # Get event time
    yyyy, mth, dd, hr, mm, sec = event.split('_')[1:]
    event_time = f'{yyyy}_{mth}_{dd}_{hr}_{mm}_{sec}'
    
    print(event_time)
    
    # Get stations that recorded this event 
    stns = [(x.split('/')[-1]).split('_')[1] for x in recordpaths]
    
    # Loop through stations 
    for j in range(len(stns)):
        
        # Get North and East components 
        recordpath_N = glob(f'{working_dir}/{wf_dir}/{event}/*_{stns[j]}_H*N_{event_time}.mseed')
        recordpath_E = glob(f'{working_dir}/{wf_dir}/{event}/*_{stns[j]}_H*E_{event_time}.mseed')
        recordpath_Z = glob(f'{working_dir}/{wf_dir}/{event}/*_{stns[j]}_H*Z_{event_time}.mseed')
        
        # Read in records
        N_record = read(recordpath_N[0])
        E_record = read(recordpath_E[0])
        Z_record = read(recordpath_Z[0])
        
        N_data = N_record[0].data
        E_data = E_record[0].data
        Z_data = Z_record[0].data
        
        # Calculate Arias intensity 
        dt = 1/(N_record[0].stats.sampling_rate)
        arias_N, normarias_N = AI.get_arias_intensity(N_data, dt)
        arias_E, normarias_E = AI.get_arias_intensity(E_data, dt)
        arias_Z, normarias_Z = AI.get_arias_intensity(Z_data, dt)
        
        # Calculate time from start of record to reach 80% of total Arias Intensity 
        time_to_80_N = AI.get_time_from_percent(normarias_N, 0.8, dt)
        time80_N = N_record[0].stats.starttime + datetime.timedelta(0,time_to_80_N)
        
        time_to_80_E = AI.get_time_from_percent(normarias_E, 0.8, dt)
        time80_E = E_record[0].stats.starttime + datetime.timedelta(0,time_to_80_E)
        
        time_to_80_Z = AI.get_time_from_percent(normarias_Z, 0.8, dt)
        time80_Z = Z_record[0].stats.starttime + datetime.timedelta(0,time_to_80_Z)
        
        ind = np.where((main_df['Name']==stns[j]) & (main_df['OrgT']==f'{yyyy}-{mth}-{dd} {hr}:{mm}:{sec}'))[0]
        Sarr = main_df['Sarr'].values[ind][0].replace(' ', 'T')
        
        # Calculate time from S-arrival to 80% energy
        time_after_S_N = round(time80_N - UTCDateTime(Sarr),3)
        time_after_S_E = round(time80_E - UTCDateTime(Sarr),3)
        time_after_S_Z = round(time80_Z - UTCDateTime(Sarr),3)
        
        # Save to log file
        f.write(f'{event}\t{stns[j]}\t{time_after_S_N}\t{time_after_S_E}\t{time_after_S_Z}\n')

f.close()
        