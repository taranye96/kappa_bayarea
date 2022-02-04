#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 22:41:11 2021

@author: tnye
"""

###############################################################################
# Script that makes a figure for each station plotting computed Kappa for each
# of the station's records against the record's hypocentral distance. The
# figue includes error bars for the data and a line is fit to the data with a 
# confidence band of 3 standard deviations. Where the line crosses 0 km distance
# is theoretically the kappa at the station. 
###############################################################################

# Imports
import numpy as np
import pandas as pd
from os import path, makedirs
from glob import glob
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from pylab import *


def plot_kappa(model_name,home_dir,df,min_events):
    
    # outfilename = f'{home_dir}/{model_name}_kappa.out'
    
    files = sorted(glob(f'{home_dir}/stn_flatfiles/*.csv'))
    
    # Make sure figure path exists
    figpath = f'{home_dir}/plots'
    if not path.exists(figpath):
        makedirs(figpath)
    
    def linearFunc(x,intercept,slope):
            y = intercept + slope * x
            return y
    
    if len(df)>=min_events:
        stn = df['Name'].iloc[0]
        mag = np.array(df['Mag'])
        rhyp = np.array(df['rhyp'])
        depth = np.array(df['Qdep'])/1000
        kappa = np.array(df['Kappa(s)'])
        kappa_std = np.array(df['Kappa std dev'])
        
        ex_ind = []
        for i, k in enumerate(kappa):
            if k <=0:
                ex_ind.append(i)
                
        # Remove events from event list that do not fit the magnitude and dist criteria
        mag = np.delete(mag, ex_ind)
        rhyp = np.delete(rhyp, ex_ind)
        kappa = np.delete(kappa, ex_ind)
        kappa_std = np.delete(kappa_std, ex_ind)
        depth = np.delete(depth, ex_ind)
        
        rrup = np.tan(np.arccos(depth/rhyp))*depth
       
        # Error fit 
        a_fit,cov=curve_fit(linearFunc,rrup,kappa,sigma=kappa_std)
        ## pull out intercept, slope, covariance:
        inter = a_fit[0]
        slope = a_fit[1]
        d_inter = np.sqrt(cov[0][0])
        d_slope = np.sqrt(cov[1][1])
        perr = np.sqrt(np.diag(cov))
        
        # Get R^2
        residuals = kappa - linearFunc(rrup, *a_fit)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((kappa-np.mean(kappa))**2)
        r2 = 1 - (ss_res / ss_tot)
        
        # stations.append(stn)
        
        # prepare confidence level curves
        nstd = 2. # to draw 2-sigma intervals (95% confidence)
        a_fit_up = a_fit + nstd * perr
        a_fit_dw = a_fit - nstd * perr
        
        x = np.linspace(0,np.max(rrup))
        fit = linearFunc(x, *a_fit)
        fit_up = linearFunc(x, *a_fit_up)
        fit_dw = linearFunc(x, *a_fit_dw)
        
        # Plot Kappa
        fig, ax = plt.subplots(1)
        rcParams['xtick.labelsize'] = 10
        rcParams['ytick.labelsize'] = 10
        rcParams['font.size']= 12
        errorbar(rrup, kappa, yerr=kappa_std, xerr=0, ecolor='gray', alpha=.25, lw=1, fmt='none', label='data')
        
        xlabel('Rrup (km)', fontsize=12)
        ylabel('Kappa (s)', fontsize=12)
        title(f'{stn}: Kappa={round(inter,4)}', fontsize=12)
        plot(x, fit, 'r', lw=1, label='best fit curve')
        scatter(rrup, kappa, c='k', s=5, label='True curve')
        ax.fill_between(x, fit_up, fit_dw, alpha=.25, label='95% Confidence Interval')
        ax.set_xlim(xmin=0)
        legend(loc='lower right',fontsize=8)
        plt.savefig(f'{figpath}/{stn}_kappa.png', dpi=300)
        plt.close()
        
        # # Plot Kappa vs magnitude
        # plt.scatter(mag, kappa)
        # plt.xlabel('Magnitude')
        # plt.ylabel('Kappa (s)')
        # plt.title(f'{stn}')
        # plt.savefig(f'{figpath}/{stn}_KvM.png', dpi=300)
        # plt.close()
    
    return(inter, d_inter, r2)
        
    # # Save kappa to text file
    # outfile = open(outfilename, 'w')
    # out = (np.array([stations, tstar, tstar_std, r_squared], dtype=object)).T
    
    # # Write value to file and save
    # outfile.write('#site \t tstar(s) \t tstar_std \t R-squared\n')
    # np.savetxt(outfile, out, fmt='%s', delimiter='\t')
    # outfile.close()
    
    
