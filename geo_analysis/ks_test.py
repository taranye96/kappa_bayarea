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

model_name = 'model1'

fault = pd.read_csv(f'/Users/tnye/kappa/data/fault_zone_stns_{model_name}.csv')
kappa_list = pd.read_csv(f'/Users/tnye/kappa/traditional_method/models/{model_name}/{model_name}_kappa.out',delimiter='\t')[' kappa(s) ']

kappa_within = kappa_list[np.where(np.array(fault['Fault_Zone']=='within'))[0]]
kappa_out = kappa_list[np.where(np.array(fault['Fault_Zone']=='outside'))[0]]

stat, pval = ks_2samp(kappa_out,kappa_within)


c_alpha = 1.36 # for significance level of 0.05
crit = c_alpha* np.sqrt((len(kappa_within)+len(kappa_out))/(len(kappa_within)*len(kappa_out)))

print(f'P-value: {round(pval,3)}')
print(f'Statistics: {round(stat, 3)}')
print(f'Critical value: {round(crit,3)}')
