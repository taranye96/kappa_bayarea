#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 13:29:56 2021

@author: tnye
"""

# Imports 
import numpy as np
import pandas as pd
from glob import glob
from os import path, makedirs
from numpy.polynomial import polynomial as P
from scipy.signal import argrelextrema
from sklearn.metrics import r2_score
import inv_fns as freqrange
from plot_stn_spectra import plot_spectra
from plot_kappa import plot_kappa
from make_gmt_kappa import make_gmt
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.colors as colors
import matplotlib.cm as cmx
from numpy.polynomial import polynomial as P

################################ Parameters ###################################

spectra_dir = '/Users/tnye/kappa/data/spectra/acc'
home_dir = '/Users/tnye/kappa/traditional_method/models/final_model'

# Read in dataframe 
event_file = '/Users/tnye/kappa/data/flatfiles/full_bay_flatfile.csv'


######################### Run Traditional Method ##############################

# Read in event df and make sure all names are strings
lst_str_cols = ['Name','Channel']
dict_dtypes = {x : 'str'  for x in lst_str_cols}
event_df = pd.read_csv(event_file, dtype=dict_dtypes)

# Get list of stations
stations = pd.read_csv('/Users/tnye/kappa/traditional_method/models/final_model/final_model_kappa.out', delimiter='\t')['#site ']
# stations = stations[:20]

# Set up figure
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'

nrows, ncols = 20, 14
dx, dy = 10,10 
# figsize = plt.figaspect(float(dy * nrows) / float(dx * ncols))
figsize = (10*ncols, 8*nrows)

fig, axs = plt.subplots(nrows, ncols, figsize=figsize, sharex='col')
# fig, axs = plt.subplots(len(stations), 2, sharex='col')

# Set up colormap
viridis = plt.get_cmap('viridis_r') 
cNorm  = colors.Normalize(vmin=3.5, vmax=5.5)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=viridis)

s = 0
i = 0
j = 0
for stn in stations:

 
    ############################### Kappa Regression ##############################
    
    # Initial station dataframe
    stn_df = pd.read_csv(f'{home_dir}/stn_flatfiles/{stn}.csv')
    
    def linearFunc(x,intercept,slope):
            y = intercept + slope * x
            return y
    
    # Kappa regression
    rhyp = np.array(stn_df['rhyp'])
    depth = np.array(stn_df['Qdep'])/1000
    kappa_list = np.array(stn_df['Kappa(s)'])
    kappa_std = np.array(stn_df['Kappa std dev'])
    
    ex_ind = []
    for k, l in enumerate(kappa_list):
        if l <=0:
            ex_ind.append(k)
    
    # Calculate epicentral distance 
    rrup = np.tan(np.arccos(depth/rhyp))*depth
       
    # Error fit 
    a_fit,cov=curve_fit(linearFunc,rrup,kappa_list,sigma=kappa_std)
    
    ## pull out intercept, slope, covariance:
    inter = a_fit[0]
    slope = a_fit[1]
    d_inter = np.sqrt(cov[0][0])
    d_slope = np.sqrt(cov[1][1])
    perr = np.sqrt(np.diag(cov))
    
    # Get R^2
    residuals = kappa_list - linearFunc(rrup, *a_fit)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((kappa_list-np.mean(kappa_list))**2)
    r2 = 1 - (ss_res / ss_tot)
    
    # prepare confidence level curves
    nstd = 3. # to draw 3-sigma intervals
    a_fit_up = a_fit + nstd * perr
    a_fit_dw = a_fit - nstd * perr
    
    x = np.linspace(0,np.max(rrup))
    fit = linearFunc(x, *a_fit)
    fit_up = linearFunc(x, *a_fit_up)
    fit_dw = linearFunc(x, *a_fit_dw)
    
    
    ################################# Make Figure #################################
    
    
    ### Record spectra subplot
    for r in range(len(stn_df)):
        
        orgt = stn_df['OrgT'].iloc[r]
        yyyy,mth,dd = orgt.split(' ')[0].split('-')
        hh,mm,sec = orgt.split(' ')[1].split(':')
        event = f'Event_{yyyy}_{mth}_{dd}_{hh}_{mm}_{sec}'
        
        files = glob(f'{spectra_dir}/{event}/*_{stn}_*')
        if len(files) > 0:
            event_file = glob(f'{spectra_dir}/{event}/*_{stn}_*')[0]
            # Read in data
            data = np.genfromtxt(event_file, comments = '#', dtype = float)
            freq = data.T[0]
            amp = data.T[1]
        
            # Obtain frequencies and amplitudes in frequency range used to obtain Kappa
            fe = stn_df['fe'].iloc[r]
            fx = stn_df['fx'].iloc[r]
            freq_ind = np.where((freq >= fe) & (freq <= fx))[0]
            org = stn_df['OrgT'].iloc[r]
            mag = stn_df['Mag'].iloc[r]
            
            A0 = stn_df['A0'].iloc[r]
            kappa = stn_df['Kappa(s)'].iloc[r]
            if kappa > 0:
                line_freq = freq[freq_ind]
                line_amp = A0*np.exp(-np.pi*kappa*line_freq)
                
                colorVal = scalarMap.to_rgba(mag)
                # axs[i,0].semilogy(freq,amp,c=colorVal,alpha=0.5,lw=.8,label=org)
                axs[i,j].semilogy(freq,amp,c='k',alpha=0.5,lw=.8)
                axs[i,j].semilogy(line_freq,line_amp,c='k',lw=.8)
    
    axs[i,j].set_xlabel('Frequency (Hz)', fontsize=12)
    axs[i,j].set_xlim(xmax=35)
    axs[i,j].text(0.95, 0.95, f'{stn}',transform=axs[i,0].transAxes,fontsize=12,va='top',ha='right')
    axs[i,j].grid(linestyle='-',alpha=0.5)
    axs[i,j].tick_params(direction="in",labelright=False,top=True,right=True,labelsize=12)
    # axs[i,0].set(aspect=1)
    
    if (j % 13 == 0) & (j != 0):
        i+=1
        j=0
    else:
        j+=1
    
    s+=1
    
    
    # ### Kappa regression subplot
    # axs[i,1].set_xlabel('R$_{epi}$ (km)', fontsize=12)
    # axs[i,1].plot(x, fit, 'r', lw=1, label='best fit curve')
    # axs[i,1].scatter(rrup, kappa_list, c='k', s=5, label='True curve')
    # axs[i,1].fill_between(x, fit_up, fit_dw, alpha=.25, label='3-sigma interval')
    # axs[i,1].set_xlim(xmin=0)
    # axs[i,1].text(0.05, 0.95, f'$\kappa_0$ = {round(inter,3)} s',transform=axs[i,1].transAxes,fontsize=12,va='top',ha='left')
    # axs[i,1].grid(linestyle='-',alpha=0.5)
    # axs[i,1].tick_params(direction="in",labelleft=False,labelright=True,top=True,right=True,labelsize=12)
    # axs[i,1].errorbar(rrup, kappa_list, yerr=kappa_std, xerr=0, ecolor='gray', alpha=.25, lw=1, fmt='none', label='data')
    # axs[i,1].set(aspect=1)
        
    # if i == length(stations)-1:
    #     cax = fig.add_axes([0.455, 0.1, 0.025, 0.825])
    #     cbar = fig.colorbar(scalarMap, label='Magnitude', cax=cax, ticks=[3.5, 4, 4.5, 5, 5.5])
    #     cbar.ax.set_yticklabels(['3.5', '4.0', '4.5', '5.0', '5.5'])
axs[0,1].legend(loc='lower right',fontsize=10)    
plt.subplots_adjust(wspace=0,hspace=0)
# plt.tight_layout()
    
plt.show()
plt.savefig('/Users/tnye/kappa/plots/paper/all_stns_spectra2.png', dpi=300)
# plt.close()
     
                
