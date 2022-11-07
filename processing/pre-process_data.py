#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 12:19:05 2022

@author: tnye
"""

###############################################################################
# This script pre-process the records by de-meaning, performing base-line
# corrections, and eliminates records with SNR < 3. 
###############################################################################

# Imports
from glob import glob
from os import path, makedirs
import numpy as np
import pandas as pd
from obspy import read
from obspy.core import UTCDateTime
from mtspec import mtspec
import tsueqs_main_fns as tmf
import filt_fns as filt


wf_dir = '/Users/tnye/kappa/data/waveforms/bay_acc_2000-2021/instrument_corrected/unprocessed'

SNR = 3

events = sorted(glob(f'{wf_dir}/*'))

############################### Process data ##################################

catalog = pd.read_csv('/Users/tnye/kappa/data/flatfiles/data_download/bay_2000-2021_flatfile.csv')
catalog['Noise_source'] = ""
catalog['SNR Window (s)'] = ""
catalog['SNR_E_freq'] = ""
catalog['SNR_N_freq'] = ""
catalog['SNR_Z_freq'] = ""
catalog['SNR_E_HFfreq'] = ""
catalog['SNR_N_HFfreq'] = ""
catalog['SNR_Z_HFfreq'] = ""
catalog['SNR_E_time'] = ""
catalog['SNR_N_time'] = ""
catalog['SNR_Z_time'] = ""

rm_ind = []

for event in events:

    e_short = event.split('/')[-1]
    print(f'Working on {e_short}')
    
    # Get list of stations recorded by this event
    stn_files = glob(f'{event}/*')
    stn_list = []
    for file in stn_files:
        stn_list.append(file.split('/')[-1].split('_')[1])
    stns = np.unique(stn_list)
    
    _,yyyy,mth,dd,hh,mm,ss = e_short.split('_')
    orgt = f'{yyyy}-{mth}-{dd} {hh}:{mm}:{ss}'
    
    # Loop through stations
    for stn in stns:
        E_file = glob(f'{event}/*{stn}_H*E_*')[0]
        N_file = glob(f'{event}/*{stn}_H*N_*')[0]
        Z_file = glob(f'{event}/*{stn}_H*Z_*')[0]
        
        stE = read(E_file)
        stN = read(N_file)
        stZ = read(Z_file)
        
        channel = stE[0].stats.channel
        samprate = stE[0].stats.sampling_rate
        npts = stE[0].stats.npts
        
        # Filter microseism
        st_filtE = stE.copy()
        st_filtN = stN.copy()
        st_filtZ = stZ.copy()
        
        st_filtE[0].data = filt.butter_highpass_filter(st_filtE[0].data, 0.28, samprate, order=5)
        st_filtN[0].data = filt.butter_highpass_filter(st_filtN[0].data, 0.28, samprate, order=5)
        st_filtZ[0].data = filt.butter_highpass_filter(st_filtZ[0].data, 0.28, samprate, order=5)
    
        # Get mean of first 10 seconds of record
        N_noise_mean = np.mean(st_filtN[0].data[:int(samprate*10)]) # first 10 seconds of noise
        E_noise_mean = np.mean(st_filtE[0].data[:int(samprate*10)]) 
        Z_noise_mean = np.mean(st_filtZ[0].data[:int(samprate*10)]) 
        
        # Make copy of stream
        stN_demean = st_filtN.copy()
        stN_demean[0].data = st_filtN[0].data - N_noise_mean
        stE_demean = st_filtE.copy()
        stE_demean[0].data = st_filtE[0].data - E_noise_mean
        stZ_demean = st_filtZ.copy()
        stZ_demean[0].data = st_filtZ[0].data - Z_noise_mean
        
        # Get the pre-event baseline
        N_baseline = tmf.compute_baseline(stN_demean)
        E_baseline = tmf.compute_baseline(stE_demean)
        Z_baseline = tmf.compute_baseline(stZ_demean)
        
        # Get the baseline corrected stream object
        N_basecorr = tmf.correct_for_baseline(stN_demean,N_baseline)
        E_basecorr = tmf.correct_for_baseline(stE_demean,E_baseline)
        Z_basecorr = tmf.correct_for_baseline(stZ_demean,Z_baseline)
        
        # Calcualte SNR using first 8 seconds of non-tapered record as nosie and 8 seconds
            # following P-wave as signal
        st_ind = np.where((catalog['Name']==stn) & (catalog['OrgT']==orgt))[0][0]
        catalog['OrgT'].iloc[st_ind]=f'{yyyy}-{mth}-{dd}T{hh}:{mm}:{ss}'
        
        # Get index of P-arrival
        parr_str = catalog['Parr'].values[st_ind]
        p_yyyy,p_mth,p_dd = parr_str.split(' ')[0].split('-')
        p_hh,p_mm,p_ss = parr_str.split(' ')[1].split(':')
        parr = UTCDateTime(int(p_yyyy),int(p_mth),int(p_dd),int(p_hh),int(p_mm),int(p_ss))
        pind = np.absolute(stN[0].times('UTCDateTime')-parr).argmin()
        
        # Get index of S-arrival
        sarr_str = catalog['Sarr'].values[st_ind]
        s_yyyy,s_mth,s_dd = sarr_str.split(' ')[0].split('-')
        s_hh,s_mm,s_ss = sarr_str.split(' ')[1].split(':')
        sarr = UTCDateTime(int(s_yyyy),int(s_mth),int(s_dd),int(s_hh),int(s_mm),int(s_ss))
        sind = np.absolute(stN[0].times('UTCDateTime')-sarr).argmin()
        
        # Get index of noise following edge of cosine taper (frac = int(npts * 0.05 / 2.0 + 0.5))
        nind = round((npts * 0.025) + 1)
        
        # Get length of SNR window based on useable noise content before P-arrival
        if (pind - nind)/samprate > 8:
            noise = 'pre_noise'
            duration = 8
        elif (pind - nind)/samprate > 7 and (pind - nind)/samprate < 8:
            noise = 'pre_noise'
            duration = 7
        elif (pind - nind)/samprate > 6 and (pind - nind)/samprate < 7:
            noise = 'pre_noise'
            duration = 6
        else:
            noise = 'coda_noise'
            duration = 8 
        
        # North component
        if noise == 'pre_noise':
            N_noise = N_basecorr[0].data[nind:nind+int(duration*samprate)]
        
        elif noise == 'coda_noise':
            N_noise = N_basecorr[0].data[(npts-nind-1)-int(duration*samprate):(npts-nind-1)]
            
        N_signal = N_basecorr[0].data[sind:sind+int(duration*samprate)]
        
        
        N_sig_amp2, N_freq =  mtspec(N_signal,delta=stN[0].stats.delta, time_bandwidth=4, number_of_tapers=7, quadratic=True)
        N_noise_amp2, N_freq =  mtspec(N_noise,delta=stN[0].stats.delta, time_bandwidth=4, number_of_tapers=7, quadratic=True)
        
        N_sig_amp = np.sqrt(N_sig_amp2)
        N_noise_amp = np.sqrt(N_noise_amp2)
        
        freq_ind = np.where((N_freq>=10) & (N_freq<=30))
        
        if len(N_signal) == len(N_noise):
            sig_avg_N = np.std(N_signal)
            noise_avg_N = np.std(N_noise)
            
            if noise_avg_N == 0:
                noise_avg_N = 1E-10
            if sig_avg_N == 0:
                sig_avg_N = 1E-10
            
            for i in range(len(N_sig_amp)):
                if N_sig_amp[i] == 0:
                    N_sig_amp[i] = 1E-10
            
            SNR_N_freq = round(np.mean(N_sig_amp/N_noise_amp),2)
            SNR_N_HFfreq = round(np.mean(N_sig_amp[freq_ind]/N_noise_amp[freq_ind]),2)
            SNR_N_time = round(sig_avg_N/noise_avg_N,2)
        else:
            SNR_N_freq = 0
            SNR_N_HFfreq = 0
            SNR_N_time = 0
    
        # East component
        if noise == 'pre_noise':
            E_noise = E_basecorr[0].data[nind:nind+int(duration*samprate)]
        
        elif noise == 'coda_noise':
            E_noise = E_basecorr[0].data[(npts-nind-1)-int(duration*samprate):(npts-nind-1)]
            
        E_signal = E_basecorr[0].data[sind:sind+int(duration*samprate)]
        
        E_sig_amp2, E_freq =  mtspec(E_signal,delta=stE[0].stats.delta, time_bandwidth=4, number_of_tapers=7, quadratic=True)
        E_noise_amp2, E_freq =  mtspec(E_noise,delta=stE[0].stats.delta, time_bandwidth=4, number_of_tapers=7, quadratic=True)
        
        E_sig_amp = np.sqrt(E_sig_amp2)
        E_noise_amp = np.sqrt(E_noise_amp2)
        
        if len(E_signal) == len(E_noise):
            sig_avg_E = np.std(E_signal)
            noise_avg_E = np.std(E_noise)
            
            if noise_avg_E == 0:
                noise_avg_E = 1E-10
            if sig_avg_E == 0:
                sig_avg_E = 1E-10
                
            for i in range(len(E_sig_amp)):
                if E_sig_amp[i] == 0:
                    E_sig_amp[i] = 1E-10
            
            SNR_E_freq = round(np.mean(E_sig_amp/E_noise_amp),2)
            SNR_E_HFfreq = round(np.mean(E_sig_amp[freq_ind]/E_noise_amp[freq_ind]),2)
            SNR_E_time = round(sig_avg_E/noise_avg_E,2)
            
        else:
            SNR_E_freq = 0
            SNR_E_HFfreq = 0
            SNR_E_time = 0
    
        # Vertical component
        if noise == 'pre_noise':
            Z_noise = Z_basecorr[0].data[nind:nind+int(duration*samprate)]
        
        elif noise == 'coda_noise':
            Z_noise = Z_basecorr[0].data[(npts-nind-1)-int(duration*samprate):(npts-nind-1)]
            
        Z_signal = Z_basecorr[0].data[sind:sind+int(duration*samprate)]
        
        Z_sig_amp2, Z_freq =  mtspec(Z_signal,delta=stZ[0].stats.delta, time_bandwidth=4, number_of_tapers=7, quadratic=True)
        Z_noise_amp2, Z_freq =  mtspec(Z_noise,delta=stZ[0].stats.delta, time_bandwidth=4, number_of_tapers=7, quadratic=True)
        
        Z_sig_amp = np.sqrt(Z_sig_amp2)
        Z_noise_amp = np.sqrt(Z_noise_amp2)
        
        if len(Z_signal) == len(Z_noise):
            sig_avg_Z = np.std(Z_signal)
            noise_avg_Z = np.std(Z_noise)
            
            if noise_avg_Z == 0:
                noise_avg_Z = 1E-10
            if sig_avg_Z == 0:
                sig_avg_Z = 1E-10
            
            for i in range(len(Z_sig_amp)):
                if Z_sig_amp[i] == 0:
                    Z_sig_amp[i] = 1E-10
            
            SNR_Z_freq = round(np.mean(Z_sig_amp/Z_noise_amp),2)
            SNR_Z_HFfreq = round(np.mean(Z_sig_amp[freq_ind]/Z_noise_amp[freq_ind]),2)
            SNR_Z_time = round(sig_avg_Z/noise_avg_Z,2)
        else:
            SNR_Z_freq = 0
            SNR_Z_HFfreq = 0
            SNR_Z_time = 0
        
        # Save SNR to flatfile
        catalog['Noise_source'].iloc[st_ind] = noise
        catalog['SNR Window (s)'].iloc[st_ind] = duration
        catalog['SNR_E_freq'].iloc[st_ind] = SNR_E_freq
        catalog['SNR_N_freq'].iloc[st_ind] = SNR_N_freq
        catalog['SNR_Z_freq'].iloc[st_ind] = SNR_Z_freq
        catalog['SNR_E_HFfreq'].iloc[st_ind] = SNR_E_HFfreq
        catalog['SNR_N_HFfreq'].iloc[st_ind] = SNR_N_HFfreq
        catalog['SNR_Z_HFfreq'].iloc[st_ind] = SNR_Z_HFfreq
        catalog['SNR_E_time'].iloc[st_ind] = SNR_E_time
        catalog['SNR_N_time'].iloc[st_ind] = SNR_N_time
        catalog['SNR_Z_time'].iloc[st_ind] = SNR_Z_time
    
        # Save processed records if SNR greater than limit
        if not path.exists(event.replace(f"{wf_dir.split('/')[-1]}",'SNR_3')):
            makedirs(event.replace(f"{wf_dir.split('/')[-1]}",'SNR_3'))
        
        # Save processed records if SNR greater than limit
        if not path.exists(event.replace(f"{wf_dir.split('/')[-1]}",'SNR_5')):
            makedirs(event.replace(f"{wf_dir.split('/')[-1]}",'SNR_5'))
        
        # Save records with SNR > 3
        if SNR_N_freq >= 3 and SNR_E_freq >= 3 and SNR_Z_freq >= 3 and SNR_N_HFfreq >= 3 and SNR_E_HFfreq >= 3 and SNR_Z_HFfreq >= 3 and SNR_N_time >= 3 and SNR_E_time >= 3 and SNR_Z_time >= 3:
            
            trE = E_basecorr[0]
            filename = E_file.replace(f"{wf_dir.split('/')[-1]}",'SNR_3')
            trE.write(filename, format='MSEED')
            
            trN = N_basecorr[0]
            filename = N_file.replace(f"{wf_dir.split('/')[-1]}",'SNR_3')
            trN.write(filename, format='MSEED')
            
            trZ = Z_basecorr[0]
            filename = Z_file.replace(f"{wf_dir.split('/')[-1]}",'SNR_3')
            trZ.write(filename, format='MSEED')
        
        # Save records with SNR > 5
        if SNR_N_freq >= 5 and SNR_E_freq >= 5 and SNR_Z_freq >= 5 and SNR_N_HFfreq >= 5 and SNR_E_HFfreq >= 5 and SNR_Z_HFfreq >= 5 and SNR_N_time >= 5 and SNR_E_time >= 5 and SNR_Z_time >= 5:
            
            trE = E_basecorr[0]
            filename = E_file.replace(f"{wf_dir.split('/')[-1]}",'SNR_5')
            trE.write(filename, format='MSEED')
            
            trN = N_basecorr[0]
            filename = N_file.replace(f"{wf_dir.split('/')[-1]}",'SNR_5')
            trN.write(filename, format='MSEED')
            
            trZ = Z_basecorr[0]
            filename = Z_file.replace(f"{wf_dir.split('/')[-1]}",'SNR_5')
            trZ.write(filename, format='MSEED')

new_catalog = catalog.drop(index=rm_ind)
new_catalog.to_csv('/Users/tnye/kappa/data/flatfiles/data_download/bay_2000-2021_SNR_irfilt_fixedtime.csv')

