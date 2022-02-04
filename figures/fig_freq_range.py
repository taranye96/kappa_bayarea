#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 15:20:13 2021

@author: tnye
"""

###############################################################################
# Script used to make the figure demonstrating the frequency range algorithm. 
###############################################################################

# Imports
import numpy as np 
import pandas as pd
from numpy.polynomial import polynomial as P
from glob import glob
from scipy.signal import argrelextrema
import matplotlib as mpl
import matplotlib.pyplot as plt
import inv_fns as freqrange
from numpy.polynomial import polynomial as P
from matplotlib.ticker import MultipleLocator, ScalarFormatter


model_name = 'model1'

# Station
stn = 'JFP'

# Read in dataframe 
lst_str_cols = ['Name','Channel']
dict_dtypes = {x : 'str'  for x in lst_str_cols}
event_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/full_bay_flatfile.csv', dtype=dict_dtypes)

# Read in a record for stn JFP
records = glob(f'/Users/tnye/kappa/traditional_method/spectra/sm_spectra/*/*_{stn}_*')
records = glob(f'/Users/tnye/kappa/data/spectra/acc/*/*_{stn}_*')
idx = 0
record_path = records[idx]
record = np.genfromtxt(records[idx])
freq = record[:,0]
amp = record[:,1]*100

# Get Magnitude for this event record
event = record_path.split('/')[-2]
yyyy,mth,dd,hh,mm,sec = event.split('_')[1:]
M = event_df['Mag'].iloc[np.where(event_df['OrgT']==f'{yyyy}-{mth}-{dd} {hh}:{mm}:{sec}')[0][0]]

# Get fe and fx
fe,fx,fc,_,_,_ = freqrange.get_freq_range(freq,amp,M)

# Get polyfit function
deg=10
coeff = P.polyfit(freq,np.log10(amp),deg)
log_pred = 0
for i in range(len(coeff)):
    log_pred += coeff[i]*(freq**i)
amp_pred = 10**log_pred

# Calcualte 1st and 2nd derivatives
diff1 = np.diff(amp_pred)
diff2 = np.diff(np.diff(amp_pred))

fc = freqrange.calc_fc(M)


################################# Make figure #################################

# Collect station flatfiles
stn_flatfiles = glob(f'/Users/tnye/kappa/traditional_method/models/{model_name}/stn_flatfiles/*.csv')

# Get fe, fx, and mag
fc_list = np.array([])
fe_list = np.array([])
fx_list = np.array([])
mag = np.array([])

for file in stn_flatfiles:
    df = pd.read_csv(file)
    fc_list = np.append(fc_list,np.array(df['fc']))
    fe_list = np.append(fe_list,np.array(df['fe']))
    fx_list = np.append(fx_list,np.array(df['fx']))
    mag = np.append(mag,np.array(df['Mag']))

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'

fig, axs = plt.subplots(2,2,figsize=(10,12))

axs[0,0].semilogy(freq,amp,label='Spectrum')
axs[0,0].semilogy(freq,amp_pred,ls='dashdot',label='Polyfit Function')
axs[0,0].vlines(fc, 0.9*np.min(amp_pred), 1.1*np.max(amp), color='k', lw=1, ls='--')
axs[0,0].vlines(fe, 0.9*np.min(amp_pred), 1.1*np.max(amp), color='k', lw=1, ls='--')
axs[0,0].vlines(fx, 0.9*np.min(amp_pred), 1.1*np.max(amp), color='k', lw=1, ls='--')
axs[0,0].text(fc+0.5,1.25*10**-2, '$f_c$',transform=axs[0,0].transData,fontsize=14,va='bottom',ha='left')
axs[0,0].text(fe+0.5,1.25*10**-2, '$f_e$',transform=axs[0,0].transData,fontsize=14,va='bottom',ha='left')
axs[0,0].text(fx+0.5,1.25*10**-2, '$f_x$',transform=axs[0,0].transData,fontsize=14,va='bottom',ha='left')
axs[0,0].set_xlim(0,35)
axs[0,0].set_ylim(0.9*np.min(np.concatenate((amp, amp_pred), axis=None)), 1.1*np.max(np.concatenate((amp, amp_pred), axis=None)))
axs[0,0].set_xlabel('Frequency (Hz)', fontsize=14)
axs[0,0].xaxis.set_minor_locator(MultipleLocator(5))
axs[0,0].tick_params(axis='x', which='minor', top=True)
axs[0,0].tick_params(axis='y', which='minor', right=True)
axs[0,0].set_ylabel('Acceleration Spectrum (cm/s)', fontsize=14)
axs[0,0].tick_params(direction="out",labelright=False,top=True,right=True,labelsize=14)
axs[0,0].grid(linestyle='--',alpha=0.25)
axs[0,0].legend(loc='upper right')
axs[0,0].text(-0.25,1,'(a)',transform=axs[0,0].transAxes,fontsize=14,va='top',ha='right')

axs[0,1].plot(freq[:149],diff1,color='purple',label='$1^{st}$ Derivative')
axs[0,1].plot(freq[:148],diff2,color='goldenrod',ls='dashdot',label='$2^{nd}$ Derivative')
axs[0,1].vlines(fe, 1.1*np.min(diff1), 1.1*np.max(diff1), color='k', lw=1, ls='--')
axs[0,1].vlines(fx, 1.1*np.min(diff1), 1.1*np.max(diff1), color='k', lw=1, ls='--')
axs[0,1].hlines(0, 0, 35, color='k', lw=0.8, ls='--', alpha=0.5)
axs[0,1].text(fe+0.5, 0.5*10**-3, '$f_e$',transform=axs[0,1].transData,fontsize=14,va='bottom',ha='left')
axs[0,1].text(fx+0.5, 0.5*10**-3, '$f_x$',transform=axs[0,1].transData,fontsize=14,va='bottom',ha='left')
axs[0,1].set_xlim(0,35)
axs[0,1].set_ylim(1.1*np.min(np.concatenate((diff1, diff2), axis=None)), 1.1*np.max(np.concatenate((diff1, diff2), axis=None)))
axs[0,1].set_xlabel('Frequency (Hz)', fontsize=14)
axs[0,1].set_ylabel('Amplitude', fontsize=14)
axs[0,1].xaxis.set_minor_locator(MultipleLocator(5))
axs[0,1].yaxis.set_minor_locator(MultipleLocator(0.001))
axs[0,1].tick_params(axis='x', which='minor', top=True)
axs[0,1].tick_params(axis='y', which='minor', right=True)
axs[0,1].tick_params(direction="out",labelright=False,top=True,right=True,labelsize=14)
axs[0,1].grid(linestyle='--',alpha=0.25)
axs[0,1].legend(loc='upper right')
axs[0,1].text(-0.2,1,'(b)',transform=axs[0,1].transAxes,fontsize=14,va='top',ha='right')
axs[0,1].ticklabel_format(axis="y", style="sci", scilimits=(0,0))

axs[1,0].scatter(fx_list, mag, marker='x', c='black', s=8, alpha=0.6, label='$f_x$')
axs[1,0].scatter(fe_list, mag, marker='^', c='steelblue', s=8, alpha=0.6, label='$f_e$')
axs[1,0].scatter(fc_list, mag, marker='.', c='gray', s=8, alpha=0.6, label='$f_c$')
axs[1,0].set_xlim(0,35)
axs[1,0].set_xlabel('Frequency (Hz)', fontsize=14)
axs[1,0].set_ylabel('Magnitude', fontsize=14)
axs[1,0].tick_params(direction="out",labelright=False,top=True,right=True,labelsize=14)
axs[1,0].xaxis.set_minor_locator(MultipleLocator(5))
axs[1,0].yaxis.set_minor_locator(MultipleLocator(0.25))
axs[1,0].tick_params(axis='x', which='minor', top=True)
axs[1,0].tick_params(axis='y', which='minor', right=True)
axs[1,0].text(-0.25,1,'(c)',transform=axs[1,0].transAxes,fontsize=14,va='top',ha='right')

handles, labels = axs[1,0].get_legend_handles_labels()
order = [2,1,0]
axs[1,0].legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc=(1.05,0.7))

axs[1,1].remove()

plt.subplots_adjust(wspace=0.35,hspace=0.35,right=0.98,bottom=0.1,left=0.125,top=0.95)

plt.show()

plt.savefig(f'/Users/tnye/kappa/plots/paper/freq_rang_{model_name}.png', dpi=300)
plt.close()


