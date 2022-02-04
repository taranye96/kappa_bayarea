#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 22:34:07 2022

@author: tnye
"""

###############################################################################
# Script that makes the supplemental figure comparing the two algorithms.
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
from matplotlib.ticker import MultipleLocator, ScalarFormatter

################################ Parameters ###################################

models = ['model2', 'model1'] 

spectra_dir = '/Users/tnye/kappa/data/spectra/acc'

# Read in dataframe 
event_file = '/Users/tnye/kappa/data/flatfiles/full_bay_flatfile.csv'

# Read in event df and make sure all names are strings
lst_str_cols = ['Name','Channel']
dict_dtypes = {x : 'str'  for x in lst_str_cols}
event_df = pd.read_csv(event_file, dtype=dict_dtypes)


# Get list of stations
stns = ['C015', 'JBNB', 'NHV']

# Set up figure
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'
# fig, axs = plt.subplots(1,2,figsize2(12,5))
fig, axs = plt.subplots(3,2,figsize=(8,20))

# Set up colormap
viridis = plt.get_cmap('viridis_r') 
cNorm  = colors.Normalize(vmin=3.5, vmax=5.5)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=viridis)

spec_labels = ['(a)', '(c)','(e)', ]
reg_labels = ['(b)','(d)','(f)']
color = ['indigo','gold']
# color = ['#1f77b4', '#ff7f0e']
scatter_m = ['o', '^']
regress_c = ['#1f77b4', '#ff7f0e']
# regress_c = ['indigo', 'gold']
m_names = ['M2', 'M1']
for i, stn in enumerate(stns):
    
    m1_df = pd.read_csv(f'/Users/tnye/kappa/traditional_method/models/model1/stn_flatfiles/{stn}.csv')
    m2_df = pd.read_csv(f'/Users/tnye/kappa/traditional_method/models/model2/stn_flatfiles/{stn}.csv')

    events = np.unique(np.append(np.array(m1_df['OrgT']), np.array(m2_df['OrgT'])))
    max_rrup = np.max(np.append(np.array(m1_df['rhyp']), np.array(m2_df['rhyp'])))

    # Plot all spectra
    for orgt in events:
        yyyy,mth,dd = orgt.split(' ')[0].split('-')
        hh,mm,sec = orgt.split(' ')[1].split(':')
        event = f'Event_{yyyy}_{mth}_{dd}_{hh}_{mm}_{sec}'
        event_file = glob(f'{spectra_dir}/{event}/*_{stn}_*')[0]
        
        # Read in data
        data = np.genfromtxt(event_file, comments = '#', dtype = float)
        freq = data.T[0]
        amp = data.T[1]*100
        
        axs[i,0].semilogy(freq,amp,c='darkgray',alpha=0.5,lw=.8)
    
    # Plot model frequency ranges
    for m, model in enumerate(models):
        model_df = pd.read_csv(f'/Users/tnye/kappa/traditional_method/models/{model}/stn_flatfiles/{stn}.csv')
    
        for j in range(len(model_df)):
            orgt = model_df['OrgT'].iloc[j]
            yyyy,mth,dd = orgt.split(' ')[0].split('-')
            hh,mm,sec = orgt.split(' ')[1].split(':')
            event = f'Event_{yyyy}_{mth}_{dd}_{hh}_{mm}_{sec}'
            event_file = glob(f'{spectra_dir}/{event}/*_{stn}_*')[0]
            
            # Obtain frequencies and amplitudes in frequency range used to obtain Kappa
            fe = model_df['fe'].iloc[np.where(np.array(model_df['OrgT'])==orgt)[0][0]]
            fx = model_df['fx'].iloc[np.where(np.array(model_df['OrgT'])==orgt)[0][0]]
            freq_ind = np.where((freq >= fe) & (freq <= fx))[0]
            org = model_df['OrgT'].iloc[np.where(np.array(model_df['OrgT'])==orgt)[0][0]]
            mag = model_df['Mag'].iloc[np.where(np.array(model_df['OrgT'])==orgt)[0][0]]
            
            A0 = model_df['A0'].iloc[np.where(np.array(model_df['OrgT'])==orgt)[0][0]]
            kappa = model_df['Kappa(s)'].iloc[np.where(np.array(model_df['OrgT'])==orgt)[0][0]]
            if kappa > 0:
                line_freq = freq[freq_ind]
                line_amp = A0*np.exp(-np.pi*kappa*line_freq)*100
                if j==0:
                    axs[i,0].semilogy(line_freq,line_amp,alpha=0.5,c=color[m],lw=.8,label=m_names[m])
                else:
                    axs[i,0].semilogy(line_freq,line_amp,alpha=0.5,c=color[m],lw=.8)
            
    axs[i,0].set_xlim(0,35)
    axs[i,0].tick_params(axis='x', which='minor')
    if i == 2:
        axs[i,0].set_xlabel('Frequency (Hz)', fontsize=12)
    axs[i,0].set_ylabel('Acc. Spectrum (cm/s)', fontsize=12)
    axs[i,0].text(-0.25,1.0,spec_labels[i],transform=axs[i,0].transAxes,fontsize=14,va='bottom',ha='right')
    axs[i,0].set_title(f'{stn}', fontsize=12)
    axs[i,0].xaxis.set_minor_locator(MultipleLocator(5))
    axs[i,0].tick_params(which='minor', right=True, top=True)
    axs[i,0].tick_params(which='major', right=True, top=True)
    if i==2:
        order = [1,0]
        handles, labels = axs[i,0].get_legend_handles_labels()
        axs[i,0].legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc="lower center", prop={"size":10}, bbox_to_anchor=(0.5, -0.6), ncol=2)
        
         
        
    ######################### Plot kappa regression #######################
    
    def linearFunc(x,intercept,slope):
        y = intercept + slope * x
        return y
    
    k0 = []
        
    # Plot model frequency ranges
    for m, model in enumerate(models):
        model_df = pd.read_csv(f'/Users/tnye/kappa/traditional_method/models/{model}/stn_flatfiles/{stn}.csv')
        
        # Kappa regression
        rhyp = np.array(model_df['rhyp'])
        depth = np.array(model_df['Qdep'])/1000
        kappa_list = np.array(model_df['Kappa(s)'])
        kappa_std = np.array(model_df['Kappa std dev'])
        
        ex_ind = []
        for l, k in enumerate(kappa_list):
            if k <=0:
                ex_ind.append(l)
        
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
        
        k0.append(inter)
        
        # Get R^2
        residuals = kappa_list - linearFunc(rrup, *a_fit)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((kappa_list-np.mean(kappa_list))**2)
        r2 = 1 - (ss_res / ss_tot)
        
        # prepare confidence level curves
        nstd = 3. # to draw 3-sigma intervals
        a_fit_up = a_fit + nstd * perr
        a_fit_dw = a_fit - nstd * perr
        
        x = np.linspace(0,200)
        fit = linearFunc(x, *a_fit)
        fit_up = linearFunc(x, *a_fit_up)
        fit_dw = linearFunc(x, *a_fit_dw)
        
        if i==2:
            axs[i,1].set_xlabel('$R_{epi}$ (km)', fontsize=12)
        axs[i,1].set_ylabel(r'$\kappa$ (s)', fontsize=12)
        axs[i,1].set_title(f'{stn}', fontsize=12)
        if i==2:
            scatlabel=f'{m_names[m]} $\kappa$ estimates'
            fillabel=f'{m_names[m]} 3-$\sigma$ interval'
            fitlabel=f'{m_names[m]} best fit curve'
        else:
            scatlabel=None
            fillabel=None
            fitlabel=None
        axs[i,1].fill_between(x, fit_up, fit_dw, alpha=.25, facecolor=regress_c[m], label=fillabel)
        axs[i,1].plot(x, fit, c=regress_c[m], lw=1, label=fitlabel)
        axs[i,1].scatter(rrup, kappa_list, marker=scatter_m[m], c='k', s=5, label=scatlabel)
        axs[i,1].set_xlim(xmin=0)
        axs[i,1].text(-0.25,1.0,reg_labels[i],transform=axs[i,1].transAxes,fontsize=14,va='bottom',ha='right')
    axs[i,1].xaxis.set_minor_locator(MultipleLocator(25))
    axs[i,1].xaxis.set_major_locator(MultipleLocator(50))
    axs[i,1].yaxis.set_minor_locator(MultipleLocator(0.025))
    axs[i,1].yaxis.set_major_locator(MultipleLocator(0.05))
    axs[i,1].tick_params(which='minor', right=True, top=True)
    axs[i,1].tick_params(which='major', right=True, top=True)
    axs[i,1].set_xlim(0,max_rrup)
    if i==2:
        order = [5,1,4,3,0,2]
        handles, labels = axs[i,1].get_legend_handles_labels()
        axs[i,1].legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc="lower center", prop={"size":10}, bbox_to_anchor=(0.5, -0.85), ncol=2)
    
    textstr = '\n'.join((
        f'M1 $\kappa_0$: {round(k0[1],3)} s',
        f'M2 $\kappa_0$: {round(k0[0],3)} s'))
    props = dict(boxstyle='round', facecolor='white', alpha=0.7)
    axs[i,1].text(0.05, 0.725, textstr, transform=axs[i,1].transAxes, fontsize=10,
        verticalalignment='bottom', bbox=props)
plt.subplots_adjust(wspace=0.4,hspace=0.35,right=0.96,bottom=0.18,left=0.125,top=0.95)

plt.savefig(f'/Users/tnye/kappa/plots/paper/S_model_comp.png', dpi=300)
plt.close()
 
            
