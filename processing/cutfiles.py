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
from obspy.core.utcdatetime import UTCDateTime

# Local Imports
import kappa_utils as ku


################################ Parameters ###################################

# Working Directory 
working_dir = '/Users/tnye/kappa/data/waveforms/bay_acc_2000-2021/instrument_corrected'

# Main flatfile for events and stations
main_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/data_download/bay_2000-2021_SNR_irfilt.csv')

energy_data = pd.read_csv('/Users/tnye/kappa/data/bay_2000-2021_energy_SNR_3_irfilt.log',delimiter="\t")

# Path to corrected seismograms
wf_dir = 'SNR_3_culled'
event_dirs = glob(working_dir + f'/{wf_dir}/Event_*')

# Path to save cut waveforms
outpath = 'cut_SNR_3_80%'

# Number of seconds before and after the S-arrival to cut the record to
cutstart = 2


############################### Cut waveforms #################################

# Get list of events 
events = []
for i in range(len(event_dirs)):
    events.append(path.basename(sorted(event_dirs)[i]))

# Make directories for cut waveforms for each event  
for i in range(len(events)):
    if not path.exists(working_dir + '/' +outpath + '/' + events[i]):
        makedirs(working_dir + '/' + outpath + '/'  + events[i])

# Loop through events 
for event in events:
    
    e = event.split('/')[-1]
    print(e)
    
    # Get list of records for this event
        # Just need one component because this is to get the list of stations
        # that recorded this event 
    recordpaths = glob(working_dir + f'/{wf_dir}/' + event + '/*' + '_H*E*.mseed')
    
    # Get event time
    yyyy, mth, dd, hr, mm, sec = event.split('_')[1:]
    event_time = f'{yyyy}_{mth}_{dd}_{hr}_{mm}_{sec}'
    
    # Get stations that recorded this event 
    stns = [(x.split('/')[-1]).split('_')[1] for x in recordpaths]
    
    # Loop through stations 
    for j in range(len(stns)):
        
        # Get North and East components 
        recordpath_N = glob(working_dir + f'/{wf_dir}/' + event + '/*_' + stns[j] + '_H*N_' + event_time + '.mseed')
        recordpath_E = glob(working_dir + f'/{wf_dir}/' + event + '/*_' + stns[j] + '_H*E_' + event_time + '.mseed')
        recordpath_Z = glob(working_dir + f'/{wf_dir}/' + event + '/*_' + stns[j] + '_H*Z_' + event_time + '.mseed')
       
        # Check that both a North and East component exist 
        if(len(recordpath_N) == 1 and len(recordpath_E) == 1):
            
            ind = np.where((main_df['Name']==stns[j]) & (main_df['OrgT']==f'{yyyy}-{mth}-{dd} {hr}:{mm}:{sec}'))[0]
            Sarr = UTCDateTime(main_df['Sarr'].values[ind][0].replace(' ', 'T'))
            
            energy_ind = np.where(energy_data['Event']==f'Event_{event_time}')[0][0]
            
            time80_E = energy_data['Sec after Sarr (E)'].values[energy_ind]
            time80_N = energy_data['Sec after Sarr (N)'].values[energy_ind]
            # time80_Z = energy_data['Sec after Sarr (Z)'].values[energy_ind]
            
            if time80_E > 0 and time80_N > 0:
                
                time80 = np.max([time80_E, time80_N])

                # Cut North component
                outpath_N = recordpath_N[0].replace(wf_dir,outpath)
                ku.cut_swave(recordpath_N[0], outpath_N, Sarr, cutstart, time80)
                 
                # Cut East component
                outpath_E = recordpath_E[0].replace(wf_dir,outpath)
                ku.cut_swave(recordpath_E[0], outpath_E, Sarr, cutstart, time80)
                
                # Cut Vertical component
                outpath_Z = recordpath_Z[0].replace(wf_dir,outpath)
                ku.cut_swave(recordpath_Z[0], outpath_Z, Sarr, cutstart, time80)
        
        