#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 22:33:42 2022

@author: tnye
"""

def compute_spectra_GM(st,nfft):
    
    import numpy as np
    from mtspec import mtspec
    
    # Compute power spectra
    E_amp2, freq, E_jack, E_fstat, E_dof =  mtspec(st[0].data,delta=st[0].stats.delta,time_bandwidth=4,nfft=nfft,number_of_tapers=7,quadratic=True,statistics=True)
    N_amp2, freq, N_jack, N_fstat, N_dof =  mtspec(st[1].data,delta=st[0].stats.delta,time_bandwidth=4,nfft=nfft,number_of_tapers=7,quadratic=True,statistics=True)
    Z_amp2, freq, Z_jack, Z_fstat, Z_dof =  mtspec(st[2].data,delta=st[0].stats.delta,time_bandwidth=4,nfft=nfft,number_of_tapers=7,quadratic=True,statistics=True)
    
    # Get FAS (L2 norm used to get averaged horizontal)
    E_amp = np.sqrt(E_amp2)
    N_amp = np.sqrt(N_amp2)
    Z_amp = np.sqrt(Z_amp2)

    NE_amp = np.sqrt(E_amp * N_amp)
    
    # Get sigma
    sigma2N = (N_jack[:,1] - N_jack[:,0])/3.29
    sigma2E = (E_jack[:,1] - E_jack[:,0])/3.29
    sigmaNE = (1/2)*np.sqrt((((sigma2N/(2*N_amp2))**2)*N_amp*E_amp) + (((sigma2E/(2*E_amp2))**2)*N_amp*E_amp))
    
    sigma2Z = (Z_jack[:,1] - Z_jack[:,0])/3.29
    sigmaZ = 0.5 * sigma2Z / Z_amp
    
    return([NE_amp, sigmaNE], [Z_amp, sigmaZ], freq)


def compute_spectra_L2norm(st,nfft):
    
    import numpy as np
    from mtspec import mtspec
    
    # Compute power spectra
    E_amp2, freq, E_jack, E_fstat, E_dof =  mtspec(st[0].data,delta=st[0].stats.delta,time_bandwidth=4,nfft=nfft,number_of_tapers=7,quadratic=True,statistics=True)
    N_amp2, freq, N_jack, N_fstat, N_dof =  mtspec(st[1].data,delta=st[0].stats.delta,time_bandwidth=4,nfft=nfft,number_of_tapers=7,quadratic=True,statistics=True)
    Z_amp2, freq, Z_jack, Z_fstat, Z_dof =  mtspec(st[2].data,delta=st[0].stats.delta,time_bandwidth=4,nfft=nfft,number_of_tapers=7,quadratic=True,statistics=True)
    
    # Get FAS (L2 norm used to get averaged horizontal)
    data_NE_2 = N_amp2 + E_amp2
    NE_amp = np.sqrt(data_NE_2)
    
    Z_amp = np.sqrt(Z_amp2)
    
    # Get sigma
    sigma2N = (N_jack[:,1] - N_jack[:,0])/3.29
    sigma2E = (E_jack[:,1] - E_jack[:,0])/3.29
    sigmaNE = np.sqrt((N_amp2/NE_amp**2.)*sigma2N + ((E_amp2/NE_amp**2.)*sigma2E))
    
    sigma2Z = (Z_jack[:,1] - Z_jack[:,0])/3.29
    sigmaZ = 0.5 * sigma2Z / Z_amp
    
    return([NE_amp, sigmaNE], [Z_amp, sigmaZ], freq)


def compute_spectra_rotd50(st,nfft):
    
    import numpy as np
    from mtspec import mtspec
    import rotd50
    
    # Get rotated horizontal components
    rotd50_horiz = rotd50.compute_rotd50(st[0].data, st[1].data)
    
    # Compute power spectra
    NE_amp2, freq, NE_jack, NE_fstat, NE_dof =  mtspec(rotd50_horiz,delta=st[0].stats.delta,time_bandwidth=4,nfft=nfft,number_of_tapers=7,quadratic=True,statistics=True)
    Z_amp2, freq, Z_jack, Z_fstat, Z_dof =  mtspec(st[2].data,delta=st[0].stats.delta,time_bandwidth=4,nfft=nfft,number_of_tapers=7,quadratic=True,statistics=True)
    
    # Get FAS
    NE_amp = np.sqrt(NE_amp2)
    Z_amp = np.sqrt(Z_amp2)
    
    # Get sigma
    sigma2NE = (NE_jack[:,1] - NE_jack[:,0])/3.29
    sigmaNE = 0.5 * sigma2NE / NE_amp
    
    sigma2Z = (Z_jack[:,1] - Z_jack[:,0])/3.29
    sigmaZ = 0.5 * sigma2Z / Z_amp
    
    return([NE_amp, sigmaNE], [Z_amp, sigmaZ], freq)


def hann_filter_time(st):
    
    import numpy as np
    from numpy.fft import rfft
    from scipy.signal import get_window
    
    for i in range(len(st)):
        # window_time = np.hanning(len(component))
        # window_freq = rfft(window_time)   
        window = get_window('hamming', st[i].stats.npts)
        st[i].data = window * st[i].data
    
    return(st)


def hann_filter_freq(spec_array):
    
    import numpy as np
    from numpy.fft import fft, rfft
    from scipy.signal import get_window
    
    smooth_spec = []
    for i in range(len(spec_array)):  
        window = get_window('hamming', len(spec_array[i]))
        window_freq = rfft(window)
        amp = np.convolve(spec_array[i],window_freq,mode='same')
        smooth_spec.append(amp)
    
    return(smooth_spec)


def compute_HVSR(filename):
    """
    Computes HVSR for a stream and returns the computed f0 with uncertainty. 
    
    Inputs:
        filename: File to stream mseed object with all 3 components.
    
    Returns:
        f0(float): fundamental site frequency 
    """
    
    import hvsrpy
    from hvsrpy import utils
    
    ################################# Parameters ##################################

    # Window length in seconds. In general low frequency peaks require longer window lengths.
    # See the SESAME guidelines for specific window length recommendations.
    windowlength = 60
    
    # Boolean to control whether Butterworth filter is applied. 
    # Geopsy does not apply a bandpass filter.
    filter_bool = False        
    # Low-cut frequency for bandpass filter.
    filter_flow = 0.1                   
    # High-cut frequency for bandpass filter.
    filter_fhigh = 30                   
    # Filter order.
    filter_order = 5
    
    # Width of cosine taper {0. - 1.}. Geopsy default of 0.05 is equal to 0.1 -> 0.1 is recommended
    width = 0.1
    
    # Konno and Ohmachi smoothing constant. 40 is recommended.
    bandwidth = 40
    
    # Minimum frequency after resampling
    resample_fmin = 0.1  
    # Maximum frequency after resampling
    resample_fmax = 50
    # Number of frequencies after resampling
    resample_fnum = 200
    # Type of resampling {'log', 'linear'}
    resample_type = 'log'
    
    # Upper and lower frequency limits to restrict peak selection. To use the entire range use `None`.
    peak_f_lower = None
    peak_f_upper = None
    
    # Method for combining horizontal components {"squared-average", "geometric-mean", "single-azimuth"}.
    # Geopsy's default is "squared-average" -> "geometric-mean" is recommended.
    method = "geometric-mean"
    # If method="single-azimuth", set azimuth in degree clock-wise from north. If method!="single-azimuth", value is ignored.
    azimuth = 0
    
    # Boolean to control whether frequency domain rejection proposed by Cox et al. (2020) is applied.
    # Geopsy does not offer this functionality.
    rejection_bool = True
    # Number of standard deviations to consider during rejection. Smaller values will reject more windows -> 2 is recommended.
    n = 2
    # Maximum number of iterations to perform for rejection -> 50 is recommended.
    max_iterations = 50
    
    # Distribution of f0 {"lognormal", "normal"}. Geopsy default "normal" -> "lognormal" is recommended.
    distribution_f0 = "lognormal"
    # Distribution of mean curve {"lognormal", "normal"}. Geopsy default "lognormal" -> "lognormal" is recommended.
    distribution_mc = "lognormal"
    
    
    ################################ Calcualte H/V ################################
    
    sensor = hvsrpy.Sensor3c.from_mseed(filename)
    bp_filter = {"flag":filter_bool, "flow":filter_flow, "fhigh":filter_fhigh, "order":filter_order}
    resampling = {"minf":resample_fmin, "maxf":resample_fmax, "nf":resample_fnum, "res_type":resample_type}
    hv = sensor.hv(windowlength, bp_filter, width, bandwidth, resampling, method, f_low=peak_f_lower, f_high=peak_f_upper, azimuth=azimuth)
    
    f0 = hv.mean_f0_frq(distribution_f0)
    
    return(f0)



def plot_HVSR(hv):
    
    import matplotlibplyplot as plt
    
    # Manually set the ylimits of the HVSR figures. Default is None so limits will be set automatically.
    ymin, ymax = 0, 10
    fig = plt.figure(figsize=(6,6), dpi=150)
    gs = fig.add_gridspec(nrows=6,ncols=6)
    
    ax0 = fig.add_subplot(gs[0:2, 0:3])
    ax1 = fig.add_subplot(gs[2:4, 0:3])
    ax2 = fig.add_subplot(gs[4:6, 0:3])
    
    if rejection_bool:
        ax3 = fig.add_subplot(gs[0:3, 3:6])
        ax4 = fig.add_subplot(gs[3:6, 3:6])
    else:
        ax3 = fig.add_subplot(gs[0:3, 3:6])
        ax4 = False
    
    individual_width = 0.3
    median_width = 1.3
    for ax, title in zip([ax3, ax4], ["Before Rejection", "After Rejection"]):
        # Rejected Windows
        if title=="After Rejection":
            if len(hv.rejected_window_indices):
                label = "Rejected"
                for amp in hv.amp[hv.rejected_window_indices]:
                    ax.plot(hv.frq, amp, color='#00ffff', linewidth=individual_width, zorder=2, label=label)
                    label=None
                
        # Accepted Windows
        label="Accepted"
        for amp in hv.amp[hv.valid_window_indices]:
            ax.plot(hv.frq, amp, color='#888888', linewidth=individual_width,
                    label = label if title=="Before Rejection" else "")
            label=None
        
        # Window Peaks
        ax.plot(hv.peak_frq, hv.peak_amp, linestyle="", zorder=2,
                marker='o', markersize=2.5, markerfacecolor="#ffffff", markeredgewidth=0.5, markeredgecolor='k',
                label="" if title=="Before Rejection" and rejection_bool else r"$f_{0,i}$")
        
        # Peak Mean Curve
        ax.plot(hv.mc_peak_frq(distribution_mc), hv.mc_peak_amp(distribution_mc), linestyle="", zorder=4,
                marker='D', markersize=4, markerfacecolor='#66ff33', markeredgewidth=1, markeredgecolor='k', 
                label = "" if title=="Before Rejection" and rejection_bool else r"$f_{0,mc}$")
        
        # Mean Curve
        label = r"$LM_{curve}$" if distribution_mc=="lognormal" else "Mean"   
        ax.plot(hv.frq, hv.mean_curve(distribution_mc), color='k', linewidth=median_width,
                label="" if title=="Before Rejection" and rejection_bool else label)
        
        # Mean +/- Curve
        label = r"$LM_{curve}$"+" ± 1 STD" if distribution_mc=="lognormal" else "Mean ± 1 STD"
        ax.plot(hv.frq, hv.nstd_curve(-1, distribution_mc),
                color='k', linestyle='--', linewidth=median_width, zorder=3,
                label = "" if title=="Before Rejection" and rejection_bool else label)
        ax.plot(hv.frq, hv.nstd_curve(+1, distribution_mc),
                color='k', linestyle='--', linewidth=median_width, zorder=3)
    
        # f0 +/- STD
        if ymin is not None and ymax is not None:
            ax.set_ylim((ymin, ymax))
        label = r"$LM_{f0}$"+" ± 1 STD" if distribution_f0=="lognormal" else "Mean f0 ± 1 STD"    
        _ymin, _ymax = ax.get_ylim()
        ax.plot([hv.mean_f0_frq(distribution_f0)]*2, [_ymin, _ymax], linestyle="-.", color="#000000")
        ax.fill([hv.nstd_f0_frq(-1, distribution_f0)]*2 + [hv.nstd_f0_frq(+1, distribution_f0)]*2, [_ymin, _ymax, _ymax, _ymin], 
                color = "#ff8080",
                label="" if title=="Before Rejection" and rejection_bool else label)
        ax.set_ylim((_ymin, _ymax))
        
        ax.set_xscale('log')
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("HVSR Amplitude")
        if rejection_bool:
            if title=="Before Rejection":
                print("\nStatistics before rejection:")
                hv.print_stats(distribution_f0)
                c_iter = hv.reject_windows(n, max_iterations=max_iterations, 
                                           distribution_f0=distribution_f0, distribution_mc=distribution_mc)
            elif title=="After Rejection":
                fig.legend(ncol=4, loc='lower center', bbox_to_anchor=(0.51, 0), columnspacing=2)
    
                print("\nAnalysis summary:")  
                display(pd.DataFrame(columns=[""], index=["Window length", "No. of windows", "Number of iterations to convergence", "No. of rejected windows"], 
                        data=[f"{windowlength}s", str(sensor.ns.nseries), f"{c_iter} of {max_iterations} allowed", str(sum(hv.rejected_window_indices))]))            
                print("\nStatistics after rejection:")
                hv.print_stats(distribution_f0)
        else:
            display(pd.DataFrame(columns=[""], index=["Window length", "No. of windows"], 
                             data=[f"{windowlength}s", str(sensor.ns.nseries)]))
            hv.print_stats(distribution_f0)
            fig.legend(loc="upper center", bbox_to_anchor=(0.77, 0.4))
            break
        ax.set_title(title)
    
    norm_factor = sensor.normalization_factor
    for ax, timerecord, name in zip([ax0,ax1,ax2], [sensor.ns, sensor.ew, sensor.vt], ["NS", "EW", "VT"]):
        ctime = timerecord.time
        amp = timerecord.amp/norm_factor
        ax.plot(ctime.T, amp.T, linewidth=0.2, color='#888888')
        ax.set_title(f"Time Records ({name})")
        ax.set_yticks([-1, -0.5, 0, 0.5, 1])
        ax.set_xlim(0, windowlength*timerecord.nseries)
        ax.set_ylim(-1, 1)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Normalized Amplitude')
        ax.plot(ctime[hv.rejected_window_indices].T, amp[hv.rejected_window_indices].T, linewidth=0.2, color="cyan")
    
    if rejection_bool:
        axs = [ax0, ax3, ax1, ax4, ax2]
    else:
        axs = [ax0, ax3, ax1, ax2]
        
    for ax, letter in zip(axs, list("abcde")):    
        ax.text(0.97, 0.97, f"({letter})", ha="right", va="top", transform=ax.transAxes, fontsize=12)
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
    
    
    fig.tight_layout(h_pad=1, w_pad=2, rect=(0,0.08,1,1))
    plt.show()
    
