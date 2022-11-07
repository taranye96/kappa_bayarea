#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 14:15:25 2018

@author: tnye
"""

###############################################################################
# This module contains functions for computing Arias intensity.
###############################################################################

# Imports mports
import numpy as np
from scipy import integrate
import scipy.constants as sp
import matplotlib.pyplot as plt


def get_arias_intensity(acc, dt, starttime=0):
    """
    Calculates Arias Intensity.

    Args:
        acc (array): Array of acceleration values in m/s/s.
        dt (float): Time between each sample in s.
        starttime (float): Time in s, after start of record, to start
            integrating. Usually P-wave arrival time or 0 for full record.

    Returns:
        Ia (array): Array of Arias intensity values in m/s with respect to
            time.
        NIa (array): Array of normalized Arias intensity values with respect to
            time. 
    """

    g = sp.g
    npts = len(acc)
    t = np.linspace(0, (npts-1)*dt, num=npts)
    
    # Choose acceleration values starting at the specificed starttime
    acc2 = acc[t >= starttime]
    
    # Calculate Arias intensity 
    Int = integrate.cumtrapz(acc2*acc2, dx=dt)
    Ia = Int * np.pi/(2*g)

    # Calculate normalized Arias intensity
    # divide arias intensity by its max value
    NIa = Ia/np.amax(Ia)
    return(Ia, NIa)


def get_time_from_percent(NIa, p, dt):
    """
    Find the closest value to the desired %Arias intensity and
    calculate the duration time associated with the percent.

    Args:
        NIa (array): Array of normalized Arias intensity values with respect to
            time.
        p (float): Percent (0 to 1) of Arias Intensity.
        dt (float): Time in between each record in s.

    Returns:
        time (float): Duration time to reach specified percent of Arias
            intensity.
    """

    npts = len(NIa)
    t = np.linspace(0, (npts-1)*dt, num=npts)

    time = t[np.argmin(np.abs(p-NIa))]
    return(time)


def get_time_for_duration_parameters(NIa, dt, p1, p2):
    """
    Find cloest values to the desired %Arias intensity range and calc total
    time for the range.

    Args:
        NIa (array): Array of normalized Arias intensity values with respect to
            time.
        dt (float): Time in between each record in s.
        p1 (float): Lower boundary of %Arias itensity range.
        p2 (float): Upper boundary of %Arias intensity range.

    Returns:
        dur_time(float): Time duration for paramteters in s.
    """

    t1 = get_time_from_percent(NIa, p1, dt)
    t2 = get_time_from_percent(NIa, p2, dt)

    dur_time = t2 - t1
    return(dur_time)


def plot_arias_intensity(Ia, dt):
    """Plots Arias intensity.

    Args:
        Ia (array): Array of Arias intensity values with respect to time.
        dt (float): time in between each record in s.

    Returns:
    """

    npts = len(Ia)
    t = np.linspace(0, (npts-1)*dt, num=npts)

    # Plot Arias Intensity
    plt.plot(t, Ia, 'k-')
    plt.xlabel('Time (s)')
    plt.ylabel('Arias Intensity (m/s)')
    plt.savefig('Arias_Intensity.png')


def plot_durations(NIa, dt, durations, ax, xlab=True):
    """Plots duration on normalized Arias intensity graph.

    Args:
        NIa (array): Array of normalized Arias intensity values with respect to
            time.
        dt (float): Time in between each record in s.
        Durations (array): Array of lists of duration parameters.

    Returns:
    """

    npts = len(NIa)
    t = np.linspace(0, (npts-1)*dt, num=npts)
    ax.plot(t, NIa, 'k-')
    if xlab:
        ax.set_xlabel('Time (s)')
    ax.set_ylabel('Norm Arias Intensity (m/s)')
    for i in range(len(durations)):
        p1 = durations[i][0]
        p2 = durations[i][1]
        t1 = get_time_from_percent(NIa, p1, dt)
        t2 = get_time_from_percent(NIa, p2, dt)
        height = (1/(len(durations)+1) * i) + 1/(len(durations)+1)
        ax.plot(t1, p1, 'ok')
        ax.plot(t2, p2, 'ok')
        ax.annotate('', xy=(t1, height), xytext=(t2, height),
                    arrowprops=dict(arrowstyle='<->'))
        label = '$D_{%i{-}%i}$' % (100 * durations[i][0],
                                   100 * durations[i][1])
        ax.text(t2, height, label, style='italic',
                horizontalalignment='left',
                verticalalignment='center')
