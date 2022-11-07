#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 22:45:49 2022

@author: tnye
"""

###############################################################################
# Module with functions for computing kappa using differint horizontal component
# averaging methods.
###############################################################################

def L2_norm_kappa(stn, M, N_amp2, E_amp2, sigmaN2, sigmaE2, freq, freq_model, samprate, deg, delta_freq):
    
    import numpy as np
    import inv_fns as freqrange
    from kappa_inversion import run_kappa_inv
    
    # Get L2-norm average
    data_NE_2 = N_amp2 + E_amp2
    data_NE = np.sqrt(data_NE_2)
    
    # Get standard deviations
    sigma_NE = np.sqrt((N_amp2/data_NE**2.)*sigmaN2 + ((E_amp2/data_NE**2.)*sigmaE2))
    
    fe,fx,fc,_,_,_ = freqrange.get_freq_range(stn,freq,data_NE,samprate,M,model=freq_model,deg=deg)
    
    if fx-fe >=delta_freq:

        # Obtain frequencies and amplitudes in frequency range used to obtain Kappa
        freq_ind = np.where((freq >= fe) & (freq <= fx))[0]
        inv_freq = freq[freq_ind]
        inv_amp = (data_NE)[freq_ind] #data from log to lin
        std = sigma_NE[freq_ind]
        
        kappa, A0, var = run_kappa_inv(inv_freq, inv_amp, std)
        
        # Append standard deviations 
        kappa_std = np.sqrt(np.abs(var[0]))#error in kappa
        A_std = np.sqrt(np.abs(A0*var[1]))
    
        return(fe,fx,fc,data_NE,kappa,kappa_std,A0,A_std)

    else:
        return('')


def rotd50_kappa(stn, M, NE_amp2, sigmaNE2, freq, freq_model, samprate, deg, delta_freq):
    
    import numpy as np
    import inv_fns as freqrange
    from kappa_inversion import run_kappa_inv
    
    
    data_NE = np.sqrt(NE_amp2)
    
    sigma_NE = 0.5 * sigmaNE2 / data_NE
    
    fe,fx,fc,_,_,_ = freqrange.get_freq_range(stn,freq,data_NE,samprate,M,model=freq_model,deg=deg)
    
    if fx-fe >=delta_freq:

        # Obtain frequencies and amplitudes in frequency range used to obtain Kappa
        freq_ind = np.where((freq >= fe) & (freq <= fx))[0]
        inv_freq = freq[freq_ind]
        inv_amp = (data_NE)[freq_ind] #data from log to lin
        std = sigma_NE[freq_ind]
        
        kappa, A0, var = run_kappa_inv(inv_freq, inv_amp, std)
        
        # Append standard deviations 
        kappa_std = np.sqrt(np.abs(var[0]))#error in kappa
        A_std = np.sqrt(np.abs(A0*var[1]))
    
        return(fe,fx,fc,kappa,kappa_std,A0,A_std)
    
    else:
        return('','','','NaN','','','')


def GM_kappa_ind(stn, M, N_amp2, E_amp2, sigmaN2, sigmaE2, freq, freq_model, samprate, deg, delta_freq):
    
    import numpy as np
    import inv_fns as freqrange
    from kappa_inversion import run_kappa_inv
    
    data_N = np.sqrt(N_amp2)
    data_E = np.sqrt(E_amp2)
    
    # Get standard deviations
    sigmaN = 0.5 * sigmaN2 / data_N
    sigmaE = 0.5 * sigmaE2 / data_E
    
    fe_N,fx_N,fc,_,_,_ = freqrange.get_freq_range(stn,freq,data_N,samprate,M,model=freq_model,deg=deg)
    fe_E,fx_E,fc,_,_,_ = freqrange.get_freq_range(stn,freq,data_N,samprate,M,model=freq_model,deg=deg)
    
    if fx_N-fe_N >=delta_freq and fx_E-fe_E >=delta_freq:

        # North component 
        # Obtain frequencies and amplitudes in frequency range used to obtain Kappa
        freq_ind_N = np.where((freq >= fe_N) & (freq <= fx_N))[0]
        inv_freq_N = freq[freq_ind_N]
        inv_amp_N = (data_N)[freq_ind_N] #data from log to lin
        std_N = sigmaN[freq_ind_N]
        
        kappa_N, A0_N, var_N = run_kappa_inv(inv_freq_N, inv_amp_N, std_N)
        
        # Append standard deviations 
        kappa_std_N = np.sqrt(np.abs(var_N[0]))#error in kappa
        A_std_N = np.sqrt(np.abs(A0_N*var_N[1]))
        
        # East component 
        # Obtain frequencies and amplitudes in frequency range used to obtain Kappa
        freq_ind_E = np.where((freq >= fe_E) & (freq <= fx_E))[0]
        inv_freq_E = freq[freq_ind_E]
        inv_amp_E = (data_E)[freq_ind_E] #data from log to lin
        std_E = sigmaE[freq_ind_E]
        
        kappa_E, A0_E, var_E = run_kappa_inv(inv_freq_E, inv_amp_E, std_E)
        
        # Append standard deviations 
        kappa_std_E = np.sqrt(np.abs(var_E[0]))#error in kappa
        A_std_E = np.sqrt(np.abs(A0_E*var_E[1]))
        
        perc_diff = (np.max([kappa_N,kappa_E]) - np.min([kappa_N,kappa_E]))/np.min([kappa_N,kappa_E])
        
        kappa_GM = np.sqrt(kappa_N * kappa_E)
        kappa_std = (1/2)*np.sqrt(((kappa_std_N**2)*kappa_E/kappa_N) + ((kappa_std_E**2)*kappa_N/kappa_E))
    
        return(kappa_GM, kappa_std, perc_diff)
    else:
        return('','','')
    

def GM_average_kappa(stn, M, N_amp2, E_amp2, sigmaN2, sigmaE2, freq, freq_model, samprate, deg, delta_freq):
    
    import numpy as np
    import inv_fns as freqrange
    from kappa_inversion import run_kappa_inv
    
    # Get Geometric mean
    N_amp = np.sqrt(N_amp2)
    E_amp = np.sqrt(E_amp2)
    data_NE = np.sqrt(N_amp * E_amp)
    
    # Get standard deviation
    sigma_NE = (1/2)*np.sqrt((((sigmaN2/(2*N_amp2))**2)*N_amp*E_amp) + (((sigmaE2/(2*E_amp2))**2)*N_amp*E_amp))
    
    fe,fx,fc,_,_,_ = freqrange.get_freq_range(stn,freq,data_NE,samprate,M,model=freq_model,deg=deg)
    
    if fx-fe >=delta_freq:

        # Obtain frequencies and amplitudes in frequency range used to obtain Kappa
        freq_ind = np.where((freq >= fe) & (freq <= fx))[0]
        inv_freq = freq[freq_ind]
        inv_amp = (data_NE)[freq_ind] #data from log to lin
        std = sigma_NE[freq_ind]
        
        kappa, A0, var = run_kappa_inv(inv_freq, inv_amp, std)
        
        # Append standard deviations 
        kappa_std = np.sqrt(np.abs(var[0]))#error in kappa
        A_std = np.sqrt(np.abs(A0*var[1]))
    
        return(fe,fx,fc,data_NE,kappa,kappa_std,A0,A_std)

    else:
        return('','','','','NaN','','','')
        # return()

