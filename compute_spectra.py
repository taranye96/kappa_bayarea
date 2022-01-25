#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 23:49:15 2022

@author: tnye
"""

###############################################################################
# Script used to obtain spectra from the seismic records. 
###############################################################################

# Standard Library Imports 
from obspy import read
from mtspec import mtspec
import os
import os.path as path
from glob import glob
import numpy as np
import pandas as pd


# Local Imports
from spec_func import bin_spec
from spec_func import bin_max_err
import eqs_main_fns as emf

######################### Set up paths and parameters #########################

min_mag = 0
max_dist = 400
nbins = 150
dtype = 'acc'

# Path to cut waveforms
events_path='/Users/tnye/kappa/data/waveforms/acc_all/Vs3.1/cut2_20'

# Event directories
event_dirs = sorted(glob(events_path + '/Event*'))

# Path to save spectra
outpath = '/Users/tnye/kappa/data/spectra/acc_35'

# Read in event catalog
catalog = '/Users/tnye/kappa/data/flatfiles/full_bay_flatfile.csv'
event_df = pd.read_csv(catalog)

# Flatfile station and event data
df_stns = np.array(event_df['Name'])
df_stlons = np.array(event_df['Slon'])
df_stlats = np.array(event_df['Slat'])
df_stelvs = np.array(event_df['Selv'])

df_origins = np.array(event_df['OrgT'])
df_hyplons = np.array(event_df['Qlon'])
df_hyplats = np.array(event_df['Qlat'])
df_hypdepths = np.array(event_df['Qdep'])

df_rhyp = np.array(event_df['rhyp'])
df_mag = np.array(event_df['Mag'])

# Make a directory for each event in outpath
events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))
for i in range(len(events)):
    if not path.exists(outpath + '/' + events[i]):
        os.makedirs(outpath + '/'  + events[i])


############################### Compute Spectra ###############################

# Loop through events 
for i, event in enumerate(events):
    
#    print('binning and fft of event: ' + event)
    
    # Get origin time
    yyyy,mth,dd,hr,mm,sec = event.split('_')[1:]
    time = f'{yyyy}_{mth}_{dd}_{hr}_{mm}_{sec}'
    
    # Get list of records for this event
        # Just need one component because this is to get the list of stations
        # that recorded this event 
    recordpaths = glob(events_path + '/' + event + '/*_*_H*N_' + time + '.sac')
    stns = [(x.split('/')[-1]).split('_')[1] for x in recordpaths]
        
    # Loop through stations 
    for j in range(len(stns)):
        
        ind = (np.where(event_df['Name']==stns[j])) and (np.where(df_origins==f'{yyyy}-{mth}-{dd} {hr}:{mm}:{sec}')[0][0])
        rhyp = event_df['rhyp'][ind]
        mag = event_df['Mag'][ind]
        depth = event_df['Qdep'][ind]/1000
        
        rrup = np.tan(np.arccos(depth/rhyp))*depth
        
        if rrup <= max_dist and mag >= min_mag:
            # Get East and North records
            recordpath_N = glob(events_path + '/' + event + '/*' + stns[j] + '_H*N_*.sac')
            recordpath_E = glob(events_path + '/' + event + '/*' + stns[j] + '_H*E_*.sac')
            
            # Make sure both a North and East component exist
            if(len(recordpath_E) == 1 and len(recordpath_N) == 1):
    
                # Get base name
                base_N = path.basename(recordpath_N[0])
                base_E = path.basename(recordpath_E[0])
        
                # Get sac info
                network = base_N.split('_')[0]
                station = base_N.split('_')[1]
                full_channel_N = base_N.split('_')[2]
                full_channel_E = base_E.split('_')[2]
                
                # Read in North stream
                N_stream = read(recordpath_N[0])
                N_tr = N_stream[0]
                N_data = N_tr.data
                
                # Read in East stream
                E_stream = read(recordpath_E[0])
                E_tr = E_stream[0]
                E_data = E_tr.data
                
                # Get channel
                channel = N_tr.stats.channel[:2]
                
                
                ###################### Check record lengths #######################
                
                # Some records have uneven horizontal components. Check if this is
                    # due to gaps in the data or a record being too short. If there
                    # are gaps, don't use the records.  If one is just shorter than 
                    # the other, shorten the long one to match. 
                    
                print(f'{N_tr.stats.station} {event}')
                
                if len(N_data)>0 and len(E_data)>0:
                
                    if len(N_data) != len(E_data):
                        
                        # Check for gaps:
                        st = E_stream.append(N_stream[0])
                        gaps = st.get_gaps()
                        
                        if gaps == []:
                            print(f'...Shortening record for {N_tr.stats.station} {event}')
                            print(f'......N={len(N_data)} E={len(E_data)}')
                        
                            # Minimum record length between the components
                            npts = np.min([len(N_data),len(E_data)])
                            
                            # Shorten the components to be the same length 
                            N_data = N_data[:npts]
                            E_data = E_data[:npts]
                            
                        else:
                            print(f'...Gaps in record for {N_tr.stats.station} {event}')
        
                    
                    ############## Begin fourier transform to get spectra #############
                    
                    # mtspec returns power spectra (square of spectra)
                    
                    #### North component 
                    N_spec_amp, N_freq , N_jack, N_fstat, N_dof =  mtspec(N_data,delta=N_tr.stats.delta, time_bandwidth=4, number_of_tapers=7, quadratic=True, statistics=True)
         			
                    # Find standard deviationnm 
                    sigmaN = (N_jack[:,1] - N_jack[:,0])/3.29
                   
                    # Get arrays of power spectra and frequencies 
                    spec_array_N = np.array(N_spec_amp)
                    freq_array_N = np.array(N_freq)
                    
                    
                    #### East component 
                    E_spec_amp, E_freq , E_jack, E_fstat, E_dof =  mtspec(E_data, delta=E_tr.stats.delta, time_bandwidth=4, number_of_tapers=7, quadratic=True, statistics=True)
               
                    # Find standard deviation
                    sigmaE = (E_jack[:,1] - E_jack[:,0])/3.29          
                    
                    # Get arrays of power spectra and frequencies 
                    spec_array_E = np.array(E_spec_amp)
                    freq_array_E = np.array(E_freq)
                    
                    # If evenly sampled
                    if(len(spec_array_E)==len(spec_array_N)):
                        
                        # spectra is power spectra so add the two components
                        data_NE_2 = spec_array_E + spec_array_N
                        freq_NE_2 = freq_array_E + freq_array_N
                        
                        # Convert from power spectra to amplitude spectra 
                        data_NE = np.sqrt(data_NE_2)
                        freq_NE = np.sqrt(freq_NE_2)
                        sigma = np.sqrt((spec_array_N/data_NE**2.)*sigmaN + ((spec_array_E/data_NE**2.)*sigmaE))
                        
                        # Bin spectra 0.1-end
                        bins, binned_data = bin_spec(data_NE, E_freq, num_bins=nbins)
                        bins_sig, binned_sig = bin_max_err(sigma, E_freq, num_bins=nbins)
        
                        # Make sure that all spectra are numbers
                        if (np.isnan(binned_data).any() == False):
                            
                            if dtype == 'vel':
                                units = 'm'
                            elif dtype == 'acc':
                                units = 'm/s'
                                
                            out_channel = f'{channel}NE'
                            print(f'channel = {out_channel}')
                            # Write to file
                            outfile = open(outpath + '/' + event + '/'+ network + '_' + station + '_' + out_channel + '_' + event + '.out', 'w')
                            data = np.array([bins, binned_data, binned_sig])
                            data = data.T
                            outfile.write(f'#bins \t \t {dtype}_spec_NE_{units} \t binned_sig \n')
                            np.savetxt(outfile, data, fmt=['%E', '%E', '%E'], delimiter='\t')
                            outfile.close()
                        else:
                            print('nans exist)')      
                        print('record exists')
                else:
                    print()
        else:
            print(f'{event} and {stns[j]} do not meet criteria')
            
