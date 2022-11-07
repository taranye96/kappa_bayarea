#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 14:10:57 2022

@author: tnye
"""

###############################################################################
# This module has a function to perform the inversion for kappa on a record. 
# This is called in kappa_methods.py.
###############################################################################

def run_kappa_inv(freq, amp, std):
     
    import numpy as np
    
    # Initialize d matrix for inversion 
        # Filled with ln amplitudes 
    d = np.zeros((len(freq), 1))
    
    # Initialize G matrix for inversion 
        # First column will be for pi*k*f
        # Second column will be 1's
    G = np.zeros((len(freq), 2))
    
    
    # Loop through frequency bins
    for j in range(len(freq)):
        d[j][0] = np.log(amp[j])#lin to ln
        G[j][0] = -1*np.pi*(freq[j])
        G[j][1] = 1
    
    
    ############################## Do inversion ###############################                          
    
    # Invert G matrix
    G_inv = np.linalg.pinv(G, rcond=1e-10)
    
    m = np.dot(G_inv,d)
    
    kappa = float(m[0])
    A0 = np.exp(float(m[1]))
    
    # Find covariance? Propagate error?
    covd = np.diag((std/amp)**2.)
    covm = np.dot((np.dot(G_inv, covd)), G_inv.T)
    var = (covm.diagonal())
    
    return(kappa, A0, var)
