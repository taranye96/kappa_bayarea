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
# histogram for the number of events per magnitude bin. This function is called
# in compute_kappa.py.
###############################################################################

# Imports 
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.ticker import MultipleLocator

def plot_spectra(freq_list,amp_list,fe_list,fx_list,mag,kappa_list,A0_list,stn_df,samprate,figpath):
    
    stn = stn_df['Name'].iloc[0]
    network = stn_df['Network'].iloc[0]
    
    # Set up figure
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.family'] = 'Helvetica'
    fig = plt.figure(figsize=(8,8))
    
    # definitions for the axes
    left, width = 0.1, 0.75
    bottom, height = 0.075, 0.6
    spacing = 0.065
    
    # Set up axes
    rect_plot = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    
    ax_plot = plt.axes(rect_plot)
    ax_plot.tick_params(direction='in', top=True, right=True)
    ax_histx = plt.axes(rect_histx)
    ax_histx.tick_params(direction='in', labelbottom=True)
    ax_histx.hist(stn_df['Mag'], range=(3.5,5.5), bins=8, color='lightgray', edgecolor='darkgray')
    
    # Set up colormap
    viridis = plt.get_cmap('viridis_r') 
    cNorm  = colors.Normalize(vmin=3.5, vmax=5.5)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=viridis)
        
    # Loop through events
    for i in range(len(freq_list)):
        
        amp = amp_list[i]
        freq = freq_list[i]
        fe = fe_list[i]
        fx = fx_list[i]
        freq_ind = np.where((freq >= fe) & (freq <= fx))[0]
        mag = stn_df['Mag'].iloc[i]
        
        A0 = A0_list[i]
        kappa = kappa_list[i]
        if kappa > 0:
            line_freq = freq[freq_ind]
            line_amp = A0*np.exp(-np.pi*kappa*line_freq)
            
            colorVal = scalarMap.to_rgba(mag)
            ax_plot.semilogy(freq,amp,c=colorVal,alpha=0.5,lw=.8)
            ax_plot.semilogy(line_freq,line_amp,c='k',lw=.8,)
    
    ax_plot.set_xlabel('Frequency (Hz)', fontsize=12)
    ax_plot.set_ylabel(r'Amplitude (cm/s$^2$)', fontsize=12)
    ax_plot.set_xlim(0,32)
    ax_plot.xaxis.set_major_locator(MultipleLocator(10))
    ax_plot.xaxis.set_minor_locator(MultipleLocator(5))
    ax_plot.tick_params(axis='y', which='minor', left=True, right=True)
    ax_plot.tick_params(which='minor', right=True, top=True, length=3)
    ax_plot.tick_params(which='major', right=True, top=True)
    ax_plot.tick_params(direction="out",labelright=False,top=True,right=True,labelsize=12)
    ax_plot.grid(linestyle='--',alpha=0.25)
    ax_plot.grid(linestyle='--',alpha=0.25)
    ax_histx.set_ylabel('No. Events', fontsize=12)
    ax_histx.set_xlabel('Magnitude', fontsize=12)
    ax_histx.set_xlim(3.5,5.5)
    ax_histx.tick_params(direction="out",labelright=False,top=True,right=True,labelsize=12)
    cax = fig.add_axes([width + .12, 0.075, 0.025, height])
    fig.colorbar(scalarMap, label='Magnitude', cax=cax)
    ax_histx.set_title(f'{network} {stn}', fontsize=14)
    
    plt.savefig(figpath, dpi=300)
    plt.close()
    
    return()
        