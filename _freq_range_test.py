#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 12:28:05 2021

@author: tnye
"""

# Imports
import numpy as np 
import pandas as pd
from numpy.polynomial import polynomial as P
from glob import glob
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
import inv_fns as freqrange

# stn = 'JFP'
stn = 'SIA'

lst_str_cols = ['Name','Channel']
dict_dtypes = {x : 'str'  for x in lst_str_cols}
event_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/full_bay_flatfile.csv', dtype=dict_dtypes)

records = glob(f'/Users/tnye/kappa/traditional_method/spectra/sm_spectra/*/*_{stn}_*')
records = glob(f'/Users/tnye/kappa/data/spectra/acc/*/*_{stn}_*')

idx = 7
record_path = records[idx]
record = np.genfromtxt(records[idx])
freq = record[:,0]
amp = record[:,1]

event = record_path.split('/')[-2]
yyyy,mth,dd,hh,mm,sec = event.split('_')[1:]
M = event_df['Mag'].iloc[np.where(event_df['OrgT']==f'{yyyy}-{mth}-{dd} {hh}:{mm}:{sec}')[0][0]]

fe,fx,amp_pred,diff1,diff2 = freqrange.get_freq_range(freq,amp,M)

if fx-fe >=10:

    
    plt.figure()
    plt.semilogy(freq,amp,label='Spectra')
    plt.semilogy(freq,amp_pred,label='Polyfit Function')
    # plt.vlines(fe, 0.9*np.min(amp), 1.1*np.max(amp), color='green', ls='--', label='fe' )
    # plt.vlines(fx, 0.9*np.min(amp), 1.1*np.max(amp), color='purple', ls='--', label='fx' )
    plt.vlines(fe, 0.9*np.min(amp_pred), 1.1*np.max(amp), color='k', lw=1, ls='--')
    plt.vlines(fx, 0.9*np.min(amp_pred), 1.1*np.max(amp), color='k', lw=1, ls='--')
    plt.ylim(0.9*np.min(np.concatenate((amp, amp_pred), axis=None)), 1.1*np.max(np.concatenate((amp, amp_pred), axis=None)))
    plt.legend()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude (m/s)')
    plt.show()
    # plt.savefig('/Users/tnye/proposals/12NCEE/freq_rang.png', dpi=300)
    # plt.close()
    
    plt.figure()
    plt.plot(freq[:149],diff1,label='1st Derivative')
    plt.plot(freq[:148],diff2,label='2nd Derivative')
    # plt.vlines(fe, 1.1*np.min(diff), 1.1*np.max(diff), color='green', ls='--', label='fe' )
    # plt.vlines(fx, 1.1*np.min(diff), 1.1*np.max(diff), color='purple', ls='--', label='fx' )
    plt.vlines(fe, 1.1*np.min(diff1), 1.1*np.max(diff1), color='k', lw=1, ls='--')
    plt.vlines(fx, 1.1*np.min(diff1), 1.1*np.max(diff1), color='k', lw=1, ls='--')
    plt.ylim(1.1*np.min(np.concatenate((diff1, diff2), axis=None)), 1.1*np.max(np.concatenate((diff1, diff2), axis=None)))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.legend()
    plt.show()
    # plt.savefig('/Users/tnye/proposals/12NCEE/freq_rang2.png', dpi=300)
    # plt.close()
    
    plt.show()