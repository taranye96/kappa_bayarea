#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 12:29:58 2021

@author: tnye
"""

###############################################################################
# Script to make a figure that plots all the event spectra for a station,
# colored by magnitude. Lines are also plotted representing the exponential
# decay using the estimated kappa for each record. The figure includes a
# histogram for the number of events per magnitude bin. 
###############################################################################

# Imports 
import numpy as np
import pandas as pd
from glob import glob
from os import path, makedirs
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from numpy.polynomial import polynomial as P

def plot_spectra(stn,home_dir,min_mag,min_events,event_df,spectra_dir,figpath):
    
    # Initial station dataframe
    if path.exists(f'{home_dir}/stn_flatfiles/{stn}.csv'):
        stn_df = pd.read_csv(f'{home_dir}/stn_flatfiles/{stn}.csv')
    
        if len(stn_df) >= min_events:
            
            # Set up figure
            mpl.rcParams['pdf.fonttype'] = 42
            mpl.rcParams['ps.fonttype'] = 42
            mpl.rcParams['font.family'] = 'Helvetica'
            fig = plt.figure(figsize=(10,8))
            
            # definitions for the axes
            left, width = 0.1, 0.65
            bottom, height = 0.075, 0.65
            spacing = 0.065
            
            # Set up axes
            rect_plot = [left, bottom, width, height]
            rect_histx = [left, bottom + height + spacing, width, 0.2]
            
            ax_plot = plt.axes(rect_plot)
            ax_plot.tick_params(direction='in', top=True, right=True)
            ax_histx = plt.axes(rect_histx)
            ax_histx.tick_params(direction='in', labelbottom=True)
            ax_histx.hist(stn_df['Mag'], bins=4)
            
            # Set up colormap
            viridis = plt.get_cmap('viridis_r') 
            cNorm  = colors.Normalize(vmin=min_mag, vmax=5.5)
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=viridis)
                
            # Loop through events
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
                    amp = data.T[1]
                
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
                        line_amp = A0*np.exp(-np.pi*kappa*line_freq)
                        
                        colorVal = scalarMap.to_rgba(mag)
                        ax_plot.semilogy(freq,amp,c=colorVal,alpha=0.5,lw=.8,label=org)
                        ax_plot.semilogy(line_freq,line_amp,c='k',lw=.8,)
                        
                        # plt.figure()
                        # plt.title(f'{event}')
                        # plt.semilogy(freq,amp,c=colorVal,alpha=0.5,lw=.8,label=org)
                        # plt.semilogy(line_freq,line_amp,c='k',lw=.8,)
                        # plt.show()
            
            ax_plot.set_xlabel('Freq (Hz)', fontsize=12)
            ax_plot.set_ylabel('Amp (m)', fontsize=12)
            ax_plot.set_xlim(xmax=35)
            ax_histx.set_ylabel('Counts', fontsize=12)
            ax_histx.set_xlabel('Magnitude', fontsize=12)
            plt.title(f'{stn}', fontsize=12)
            ax_histx.set_xlim(min_mag,5.5)
            cax = fig.add_axes([0.77, 0.075, 0.025, 0.65])
            fig.colorbar(scalarMap, label='Magnitude', cax=cax)
            
            plt.savefig(figpath, dpi=300)
            plt.close()
        
        return()
    
        