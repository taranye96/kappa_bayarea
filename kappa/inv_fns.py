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
    '''
    Computes the predicted amplitude based on a polyfit function
    '''

    log_pred = 0
    for i in range(len(coeff)):
        log_pred += coeff[i]*(freq**i)
    amp_pred = 10**log_pred
    
    return(amp_pred)


def calc_fc(M,stressdrop=5e6,Vs=3500):
    '''
    Calculates theoretical corner frequenchy (fc)
    '''
    
    M0 = 10.**((3./2.)*M + 9.1)
    fc = Vs*(stressdrop/(8.47*M0))**(1./3.)
    
    return(fc)


def get_freq_range(stn,freq,amp,samprate,M,model,deg=10):
    '''
    Computes the frequency limits, fe and fx, which are used to solve for kappa
    '''
    
    import numpy as np
    import pandas as pd
    from numpy.polynomial import polynomial as P
    from sklearn.preprocessing import minmax_scale
    
    # Read in event df and make sure all names are strings
    hvsr_df = pd.read_csv('/Users/tnye/kappa/HVSR/HVSR_GM.csv')
    
    stn_ind = np.where(hvsr_df['Station']==stn)[0][0]
    lower_lim = hvsr_df['lower_lim'].iloc[stn_ind]
    upper_lim = hvsr_df['upper_lim'].iloc[stn_ind]
    
    # Get polyfit function
    coeff = P.polyfit(freq,np.log10(amp),deg)
    
    # Get predicted amplitude from polyfit
    amp_pred = poly(coeff,freq)

    # Calcualte 1st and 2nd derivatives
    diff1 = np.diff(amp_pred)
    diff2 = np.diff(np.diff(amp_pred))
    
    # Get theoretical fc
    fc = calc_fc(M)
    
    # Get frequency closest to f0
    lower_ind = np.where((np.abs(freq-lower_lim))==np.min(np.abs(freq-lower_lim)))[0][0]
    
    # Get fe
    fn80 = np.where((np.abs(freq-(0.8*samprate/2)))==np.min(np.abs(freq-(0.8*samprate/2))))[0][0]
    
    fe_ind = np.where(diff1==np.min(diff1[lower_ind:fn80]))[0][0] 
    min_ind = np.where(diff1==np.min(diff1))[0][0] # most negative slope 
    fe = freq[fe_ind] 
    if fe < 1.5*fc:  #ktenidou et al., 2016, PEER
        fe = 1.5*fc
    if fe < lower_lim:
        fe = lower_lim
    
    
    ################## Models for selecting frequency ranges ##################
    
    # Note: Various models were tested for selecting the upper limit, fx. The
        # preferred model is Model 4 (m4) where fx the frequency at which the
        # 1st derivative reaches a threshold value of -8*10**-5, assuming the 
        # acceleration spectra are in cm/s**2/
    
    if model=='m1':
        '''
        Fx = where 2nd derivative reaches 0
        '''
        
        # Get fx
        if len(np.where(diff2[fe_ind:fn80] < 0)[0])>0:
            fx_ind = (fe_ind+1) + np.where(diff2[fe_ind:fn80] < 0)[0][0]
        else:
            fx_ind = np.where(np.abs(diff2[:fn80]) == np.min(np.abs(diff2[fe_ind:fn80])))[0][0]
        
        fx = freq[fx_ind]
    
    elif model=='m2':
        '''
        Fx = some threshold based on difference between 1st and 2nd derivatives
        '''
        
        # Get fx
            # If using spectra binned up to 40 Hz
        if len(np.where(np.abs(diff2[fe_ind:fn80]-diff1[fe_ind:fn80]) <= 1*10**-10)[0])>0:
            fx_ind = fe_ind+1 + np.where(np.abs(diff2[fe_ind:fn80]-diff1[fe_ind:fn80]) <= 1*10**-10)[0][0]
            fx = freq[fx_ind]
        else:
            fx_ind = np.where(np.abs(diff2-diff1[:fn80]) == np.min(np.abs(diff2[fe_ind:fn80]-diff1[fe_ind:fn80])))[0][0]
        fx = freq[fx_ind]
    
    elif model=='m3':
        '''
        Fx = some threshold based on 2nd derivative only
        '''
        
        # Get fx
        thresh = 1*10**-12
        diff_peak_ind = fe_ind+1 +np.where(diff2[fe_ind:fn80]==np.max(diff2[fe_ind:fn80]))[0][0] # gets values past peak in 2nd derivative
        if len(np.where(diff2[diff_peak_ind:fn80] <= thresh)[0])>0:
            fx_ind = diff_peak_ind+1 + np.where(diff2[diff_peak_ind:fn80] <= thresh)[0][0]
            fx = freq[fx_ind]
        else:
            fx = np.min([freq[fn80],35])
    
    elif model=='m4':
        '''
        Fx = some threshold based on 1st derivative only
        ** This is the preferred model **
        '''
        
        # Get fx
        thresh = -8*10**-5
        if len(np.where(diff1[fe_ind:fn80] >= thresh)[0])>0:
            fx_ind = fe_ind+1 + np.where(diff1[fe_ind:fn80] >= thresh)[0][0]
            fx = freq[fx_ind]
        else:
            fx = np.min([freq[fn80],35])
    
    elif model=='m5':
        '''
        Fx = some threshold based on normalized 1st derivative only
        '''
        
        # Get fx
        thresh = 0.45
        diff1_norm = minmax_scale(diff1)
        if len(np.where(diff1_norm[fe_ind:fn80] >= thresh)[0])>0:
            fx_ind = fe_ind+1 + np.where(diff1_norm[fe_ind:fn80] >= thresh)[0][0]
            fx = freq[fx_ind]
        else:
            fx = np.min([freq[fn80],35])
    
    elif model=='m6':
        '''
        Fx = some % of 1st deriv min
        '''
        
        thresh = 0.004
        if len(np.where(diff1[fe_ind:fn80]/diff1[fe_ind] <= thresh)[0])>0:
            fx_ind = fe_ind+1 + np.where(diff1[fe_ind:fn80]/diff1[fe_ind] <= thresh)[0][0]
            fx = freq[fx_ind]
        else:
            fx = np.min([freq[fn80],35])
            
    elif model=='m7':
        '''
        Fx = 0.5 % of 1st deriv min
        Impose lower and upper limits based on HVSR log
        '''
        
        thresh = 0.0075
        if len(np.where(diff1[fe_ind:fn80]/diff1[fe_ind] <= thresh)[0])>0:
            fx_ind = fe_ind+1 + np.where(diff1[fe_ind:fn80]/diff1[min_ind] <= thresh)[0][0]
            fx = freq[fx_ind]
        else:
            fx = np.min([freq[fn80],35])
    
    if fx > upper_lim:
        fx = upper_lim

    
    return(fe, fx, fc, amp_pred, diff1, diff2)
