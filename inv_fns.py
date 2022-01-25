#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 14:54:37 2021

@author: tnye
"""

###############################################################################
# Module with functions used for the frequency range selection algorithm. 
###############################################################################


def poly(coeff,freq):

    log_pred = 0
    for i in range(len(coeff)):
        log_pred += coeff[i]*(freq**i)
    amp_pred = 10**log_pred
    
    return(amp_pred)


def calc_fc(M,stressdrop=5e6,Vs=3100):
    
    M0 = 10.**((3./2.)*M + 9.1)
    fc = Vs*(stressdrop/(8.47*M0))**(1./3.)
    
    return(fc)


def get_freq_range(freq,amp,M,model='m1',deg=10):
    
    import numpy as np
    from numpy.polynomial import polynomial as P
    
    # Get polyfit function
    coeff = P.polyfit(freq,np.log10(amp),deg)
    
    # Get predicted amplitude from polyfit
    amp_pred = poly(coeff,freq)

    # Calcualte 1st and 2nd derivatives
    diff1 = np.diff(amp_pred)
    diff2 = np.diff(np.diff(amp_pred))
    
    # Get theoretical fc
    fc = calc_fc(M)
    
    # Get fe
    fe_ind = np.where(diff1==np.min(diff1[:146]))[0][0] # most negative slope [~0-35Hz] if using bins up to 40 Hz
    # fe_ind = np.where(diff1==np.min(diff1))[0][0] # most negative slope [~0-35Hz] if using bins up to 35 HZ
    fe = freq[fe_ind] 
    if fc > fe:
        fe = fc
    
    
    ################################# Model 1 #################################
    
    if model=='m1':
        # Get fx
            # If using spectra binned up to 40 Hz
        if len(np.where(np.abs(diff2[fe_ind:146]-diff1[fe_ind:146]) <= 2*10**-6)[0])>0:
            fx_ind = fe_ind+1 + np.where(np.abs(diff2[fe_ind:146]-diff1[fe_ind:146]) <= 2*10**-6)[0][0]
            fx = freq[fx_ind]
        else:
            fx_ind = np.where(np.abs(diff2-diff1[:148]) == np.min(np.abs(diff2[fe_ind:148]-diff1[fe_ind:148])))[0][0]
        fx = freq[fx_ind]
        
        # # Get fx
        #     # If using spectra binned up to 35 Hz
        # if len(np.where(np.abs(diff2[100:]-diff1[100:148]) <= 2*10**-6)[0])>0:
        #     fx_ind = 101 + np.where(np.abs(diff2[100:148]-diff1[100:148]) <= 2*10**-6)[0][0]
        #     fx = freq[fx_ind]
        # else:
        #     fx_ind = np.where(np.abs(diff2-diff1[:148]) == np.min(np.abs(diff2[100:148]-diff1[100:148])))[0][0]
        # fx = freq[fx_ind]
    
    ################################# Model 2 #################################
    
    elif model=='m2':
        # Get fx
            # # If using spectra binned up to 40 Hz
        # fx_ind = fe_ind+1 + np.where(np.abs(diff2[fe_ind:146] - 0))[0].argmin()
        if len(np.where(diff2[fe_ind:146] < 0)[0])>0:
            fx_ind = (fe_ind+1) + np.where(diff2[fe_ind:146] < 0)[0][0]
        else:
            fx_ind = np.where(np.abs(diff2[:146]) == np.min(np.abs(diff2[fe_ind:146])))[0][0]
        
        fx = freq[fx_ind]
        
        # # Get fx
        #     # If using spectra binned up to 35 Hz
        # fx_ind = 101 + np.where(np.abs(diff2[100:148] - 0))[0].argmin()
        # fx = freq[fx_ind]
    
    if fx > 35:
        fx = 35

    
    return(fe, fx, fc, amp_pred, diff1, diff2)
