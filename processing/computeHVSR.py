#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 10:19:35 2022

@author: tnye
"""

###############################################################################
# Script used to compute HVSR
###############################################################################

# Imports
from os import path, makedirs
import numpy as np
import pandas as pd
from glob import glob
from obspy import read, Stream
from obspy.core.utcdatetime import UTCDateTime
import hvsr
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator

# Read in flatfile with event and station inforamtion 
catalog = pd.read_csv('/Users/tnye/kappa/data/flatfiles/data_download/bay_2000-2021_SNR_irfilt_fixedtime.csv')

# Read in dataframe with manually-selected frequency ranges
    # Note: This is for the final plots
freq_df = pd.read_csv('/Users/tnye/kappa/HVSR/HVSR_GM.csv')

# Get all mseed files
wf_dir = '/Users/tnye/kappa/data/waveforms/bay_acc_2000-2021/instrument_corrected/cut_SNR_3_80%'
mseeds = sorted(glob(f'{wf_dir}/*/*'))

# Make directory to store hvsr
if not path.exists('/Users/tnye/kappa/data/waveforms/bay_acc_2000-2021/hvsr-GM'):
    makedirs('/Users/tnye/kappa/data/waveforms/bay_acc_2000-2021/hvsr-GM')

# Get list of stations
stns = []
for file in mseeds:
    stn = file.split('/')[-1].split('_')[1]
    stns.append(stn)
stns = np.unique(stns)

# Loop through stations and compute HVSR
for stn in stns:
    
    print(stn)
    
    # Gather all mseed files for the each station
    stn_files = sorted(glob(f'{wf_dir}/*/*_{stn}_*'))
    
    tr_length = []
    for file in stn_files:
        tr_length.append(read(file)[0].stats.npts)
    
    max_length = np.max(tr_length)
    
    # Get list of events recorded on each station
    events = []
    for file in stn_files:
        events.append(file.split('/')[-2])
    events = np.unique(events)
    
    HV_j_list = []
    sigmaHV_j_list = []
    
    # Loop through events to create a stream for each station/event pair
    for i, event in enumerate(events):
        
        yyyy,mth,dd,hh,mm,ss = event.split('_')[1:]
        
        ind = np.where((catalog['Name']==stn) & (catalog['OrgT']==f'{yyyy}-{mth}-{dd}T{hh}:{mm}:{ss}'))[0][0]
        Sarr = catalog['Sarr'].iloc[ind]     
        start = UTCDateTime(Sarr)
        end = start + 10
        
        # Get list of mseed files for this station/event
        ev_files = sorted(glob(f'{wf_dir}/{event}/*_{stn}_*'))
        
        network = ev_files[0].split('/')[-1].split('_')[0]
        
        tr_E = read(ev_files[0])[0]
        tr_N = read(ev_files[1])[0]
        tr_Z = read(ev_files[2])[0]
        
        st = Stream()
        for ev in ev_files:
            tr = read(ev)[0]
            pad_length = max_length-tr.stats.npts
            end_pad = round(pad_length/2)
            start_pad = pad_length - end_pad
            # tr.data = np.pad(tr.data, (0,max_length-tr.stats.npts), 'constant')
            tr.data = np.pad(tr.data, (start_pad,end_pad), 'constant')
            st.append(tr)
        
        # Smooth time series with Hanning taper window
        st_smooth = hvsr.hann_filter_time(st)
        
        # Compute Geometric mean of spectra
        [NE_amp, sigmaNE], [Z_amp, sigmaZ], freq = hvsr.compute_spectra_GM(st_smooth,nfft=max_length)
        
        # Calculate hvsr
        HV_ji = NE_amp / Z_amp
        HV_j_list.append(HV_ji)
        sigmaHV_ji = HV_ji * np.sqrt((sigmaNE/NE_amp)**2 + (sigmaZ/Z_amp)**2)
        
        sigmaHV_j_list.append(sigmaHV_ji)
        
        # Save smoothed waveform
        filename = f'/Users/tnye/kappa/data/waveforms/bay_acc_2000-2021/hvsr-GM/{network}_{stn}_{event}.mseed'
        st_smooth.write(filename, format='MSEED')
    
    # Find index of frequency closest to 35 Hz
    f35 = np.where(np.abs(freq-35)==np.min(np.abs(freq-35)))[0][0]
    
    # Compute earthquake hvsr
    eHVSR = 10**(sum([np.log10(i) for i in HV_j_list])/len(HV_j_list))
    
    # Compute f0
    f0_ind = np.where(eHVSR[:f35]==np.max(eHVSR[:f35]))[0][0]
    f0 = freq[f0_ind]
    
    # Error propagation of logarithmic average of H/V at a site
    sigmaHV_j = np.log(10) * eHVSR * (1/len(sigmaHV_j_list)) * np.sqrt(sum([(0.434*sigma/HV_j_list[i])**2 for i, sigma in enumerate(sigmaHV_j_list)]))
    
    # Error bars
    pos_error = eHVSR + sigmaHV_j
    neg_error = eHVSR - sigmaHV_j
    
    
    ############################### Save plots ################################
    
    fig_dir = '/Users/tnye/kappa/HVSR/GM_2'
    if path.exists(fig_dir)==False:
        makedirs(fig_dir)
        
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.family'] = 'DejaVu Sans'
    
    ind = np.where(freq_df['Station']==stn)[0][0]
    lower_lim = freq_df['lower_lim'].iloc[ind]
    upper_lim = freq_df['upper_lim'].iloc[ind]
    
    if np.max(eHVSR)<=15:
        ymax = 15
    else:
        ymax = np.max(eHVSR+1)
    
    fig, ax = plt.subplots(1,1,figsize=(6,4))
    ax.plot(freq,eHVSR,c='k',lw=1)
    ax.plot(freq,pos_error,ls='--',c='gray',lw=0.5,label='Error')
    ax.plot(freq,neg_error,ls='--',c='gray',lw=0.5)
    ax.scatter(f0,eHVSR[f0_ind],marker='o',facecolor='none',edgecolor='r',label=r'$f_0$'+f'={round(f0,2)} Hz')
    ax.vlines(lower_lim,0,15,ls='--',lw=0.6,color='r',label='frequency limit')
    ax.vlines(upper_lim,0,15,ls='--',lw=0.6,color='r')
    ax.set_xlim(0,35.5)
    ax.set_ylim(0,ymax)
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.tick_params(axis='x', which='minor', top=True)
    ax.tick_params(axis='y', which='minor', right=True)
    ax.tick_params(direction="out",labelright=False,top=True,right=True,labelsize=14)
    ax.grid(linestyle='--',alpha=0.25)
    ax.set_xlabel('Frequency (Hz)',fontsize=14)
    ax.set_ylabel('Average eHVSR',fontsize=14)
    ax.set_title(f'{network} {stn}',fontsize=14)
    ax.legend()
    plt.subplots_adjust(wspace=0.35,hspace=0.35,right=0.95,bottom=0.15,left=0.125,top=0.925)
   
    plt.savefig(f'{fig_dir}/{stn}.png',dpi=300)
    plt.close()