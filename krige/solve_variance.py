#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 14:03:14 2021

@author: tnye
"""

###############################################################################
# Script that calculates the variance of the site kappa estimates to determine
# what the nugget should be in the kirging. 
###############################################################################

# Imports
import numpy as np
import pandas as pd

model_name = 'model1'
home_dir = f'/Users/tnye/kappa/traditional_method/models/{model_name}'

std_dev = np.array(pd.read_csv(f'{home_dir}/{model_name}_kappa.out', delimiter='\t')[' kappa_std '])
kappa = np.array(pd.read_csv(f'{home_dir}/{model_name}_kappa.out', delimiter='\t')[' kappa(s) '])

perc_std = std_dev/kappa
avg_perc = np.mean(perc_std)*100

nugget = np.log10(1+(avg_perc/100))**2
print(f'nugget = {round(nugget,4)}')
