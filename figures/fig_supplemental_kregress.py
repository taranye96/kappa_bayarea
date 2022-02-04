#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 15:51:26 2022

@author: tnye
"""

###############################################################################
# Script used to make the supplemental regression figures for each station.
###############################################################################

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
from plot_kappa_regression import plot_kappa
from make_gmt_kappa import make_gmt
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.colors as colors
import matplotlib.cm as cmx
from numpy.polynomial import polynomial as P
import matplotlib as mpl

################################ Parameters ###################################

model_name = 'model1'

spectra_dir = '/Users/tnye/kappa/data/spectra/acc'
home_dir = f'/Users/tnye/kappa/traditional_method/models/{model_name}'

# Read in dataframe 
event_file = '/Users/tnye/kappa/data/flatfiles/full_bay_flatfile.csv'


######################### Run Traditional Method ##############################

# Read in event df and make sure all names are strings
lst_str_cols = ['Name','Channel']
dict_dtypes = {x : 'str'  for x in lst_str_cols}
event_df = pd.read_csv(event_file, dtype=dict_dtypes)

# Get list of stations
stn_dfs = glob(f'/Users/tnye/kappa/traditional_method/models/{model_name}/stn_flatfiles/*.csv')

for df in stn_dfs:
 
    # Initial station dataframe
    stn_df = pd.read_csv(df)
    stn = stn_df['Name'][0]
    
    ############################# Do Regression ####################33#########
    
    def linearFunc(x,intercept,slope):
            y = intercept + slope * x
            return y
    
    # Kappa regression
    rhyp = np.array(stn_df['rhyp'])
    depth = np.array(stn_df['Qdep'])/1000
    kappa_list = np.array(stn_df['Kappa(s)'])
    kappa_std = np.array(stn_df['Kappa std dev'])
    
    ex_ind = []
    for i, k in enumerate(kappa_list):
        if k <=0:
            ex_ind.append(i)
    
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
    
    # Set up figure
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.family'] = 'Helvetica'
    fig, axs = plt.subplots(1,2,figsize=(14,5))
    
    # Set up colormap
    viridis = plt.get_cmap('viridis_r') 
    cNorm  = colors.Normalize(vmin=3.5, vmax=5.5)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=viridis)
    
    ### Record spectra subplot
    for i in range(len(stn_df)):
        
        orgt = stn_df['OrgT'].iloc[i]
        yyyy,mth,dd = orgt.split(' ')[0].split('-')
        hh,mm,sec = orgt.split(' ')[1].split(':')
        event = f'Event_{yyyy}_{mth}_{dd}_{hh}_{mm}_{sec}'
        
        files = glob(f'{spectra_dir}/{event}/*_{stn}_*')
        if len(files) > 0:
            event_file = glob(f'{spectra_dir}/{event}/*_{stn}_*')[0]
            # Read in data
            data = np.genfromtxt(event_file, comments = '#', dtype = float)
            freq = data.T[0]
            amp = data.T[1]*100
        
            # Obtain frequencies and amplitudes in frequency range used to obtain Kappa
            fe = stn_df['fe'].iloc[i]
            fx = stn_df['fx'].iloc[i]
            freq_ind = np.where((freq >= fe) & (freq <= fx))[0]
            org = stn_df['OrgT'].iloc[i]
            mag = stn_df['Mag'].iloc[i]
            
            A0 = stn_df['A0'].iloc[i]
            kappa = stn_df['Kappa(s)'].iloc[i]
            if kappa > 0:
                line_freq = freq[freq_ind]
                line_amp = A0*np.exp(-np.pi*kappa*line_freq)*100
                
                colorVal = scalarMap.to_rgba(mag)
                axs[0].semilogy(freq,amp,c=colorVal,alpha=0.5,lw=.8,label=org)
                axs[0].semilogy(line_freq,line_amp,c='k',lw=.8,)
    
    axs[0].set_xlabel('Frequency (Hz)', fontsize=12)
    axs[0].set_ylabel('Acceleration Spectrum (cm/s)', fontsize=12)
    axs[0].set_xlim(0,35)
    axs[0].set_title(f'{stn}', fontsize=12)
    axs[0].xaxis.set_minor_locator(MultipleLocator(5))
    axs[0].tick_params(which='minor', right=True, top=True)
    axs[0].tick_params(which='major', right=True, top=True)
    cax = fig.add_axes([0.47, 0.1, 0.025, 0.825])
    cbar = fig.colorbar(scalarMap, label='Magnitude', cax=cax, ticks=[3.5, 4, 4.5, 5, 5.5])
    cbar.ax.set_yticklabels(['3.5', '4.0', '4.5', '5.0', '5.5'])
    axs[0].text(-0.15,1.0,'(a)',transform=axs[0].transAxes,fontsize=14,va='bottom',ha='right')
    
    ### Kappa regression subplot
    axs[1].set_xlabel('$R_{epi}$ (km)', fontsize=12)
    axs[1].set_ylabel('Kappa (s)', fontsize=12)
    axs[1].set_title(f'{stn}', fontsize=12)
    axs[1].plot(x, fit, 'r', lw=1, label='best fit curve')
    axs[1].scatter(rrup, kappa_list, c='k', s=5, label='True curve')
    axs[1].fill_between(x, fit_up, fit_dw, alpha=.25, label='3-sigma interval')
    axs[1].set_xlim(xmin=0)
    axs[1].text(0.05, 0.95, f'$\kappa_0$ = {round(inter,3)} s',transform=axs[1].transAxes,fontsize=11,va='top',ha='left')
    axs[1].legend(loc='lower right',fontsize=10)
    axs[1].errorbar(rrup, kappa_list, yerr=kappa_std, xerr=0, ecolor='gray', alpha=.25, lw=1, fmt='none', label='data')
    
    axs[1].xaxis.set_minor_locator(MultipleLocator(25))
    axs[1].xaxis.set_major_locator(MultipleLocator(50))
    axs[1].yaxis.set_minor_locator(MultipleLocator(0.025))
    axs[1].yaxis.set_major_locator(MultipleLocator(0.05))
    axs[1].tick_params(which='minor', right=True, top=True)
    axs[1].tick_params(which='major', right=True, top=True)
    axs[1].text(-0.15,1.0,'(b)',transform=axs[1].transAxes,fontsize=14,va='bottom',ha='right')
     
    plt.subplots_adjust(wspace=0.425,right=0.99,bottom=0.1,left=0.08,top=0.95)
    
    # Make directories 
    if not path.exists(f'/Users/tnye/kappa/plots/paper/supplemental_{model_name}'):
        makedirs(f'/Users/tnye/kappa/plots/paper/supplemental_{model_name}')
    
    plt.savefig(f'/Users/tnye/kappa/plots/paper/supplemental_{model_name}/{stn}.png', dpi=300)
    plt.close()
     
                
