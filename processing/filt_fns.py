#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 23:24:38 2022

@author: tnye
"""

###############################################################################
# This module contains functions for filtering data. They are called in 
# pre-process_data.py and parallel_pre-process_data.py.
###############################################################################

# Imports
import numpy as np
from obspy import read
from scipy.signal import butter, filtfilt

def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y
