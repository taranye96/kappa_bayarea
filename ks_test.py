#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 12:41:51 2022

@author: tnye
"""

###############################################################################
# Script used to run a Kolmogorov–Smirnov (K-S) test on the kappa fault zone 
# data. 
###############################################################################

# Imports 
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

fault = pd.read_csv('/Users/tnye/kappa/data/fault_zone_stns.csv')
kappa_list = pd.read_csv('/Users/tnye/kappa/traditional_method/models/test/updated_stns/updated_stns_kappa.out',delimiter='\t')[' kappa(s) ']

kappa_within = kappa_list[np.where(np.array(fault['Fault_Zone']=='within'))[0]]
kappa_out = kappa_list[np.where(np.array(fault['Fault_Zone']=='outside'))[0]]

ks_2samp(kappa_out,kappa_within)


c_alpha = 1.36 # for significance level of 0.05
crit = c_alpha* np.sqrt((len(kappa_within)+len(kappa_out))/(len(kappa_within)*len(kappa_out)))
