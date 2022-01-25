#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 17:27:51 2021

@author: tnye
"""

# Imports
import numpy as no
import pandas as pd

###############################################################################
# Script that makes a GMT text file with the kriged kappa data. 
###############################################################################

model_name = 'model1'

df = pd.read_csv(f'/Users/tnye/kappa/krige/{model_name}_krige_k0_lonlat.csv')

f = open(f'/Users/tnye/kappa/krige/{model_name}_krige_k0_lonlat.txt','w')
f.write('#Longitude,Latitude,kappa(s)\n')
for i in range(len(df)):
    f.write(f"{df['Longitude'].iloc[i]},{df['Latitude'].iloc[i]},{df['pred_k0'].iloc[i]}\n")
f.close()

f = open(f'/Users/tnye/kappa/krige/{model_name}_pred_stddev_lonlat.txt','w')
f.write('#Longitude,Latitude,k0_stddev\n')
for i in range(len(df)):
    f.write(f"{df['Longitude'].iloc[i]},{df['Latitude'].iloc[i]},{df['k0_stddev'].iloc[i]}\n")
f.close()
