#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 17:27:51 2021

@author: tnye
"""

###############################################################################
# Script that makes a GMT text file with the kriged kappa data. 
###############################################################################

# Imports
import numpy as np
import pandas as pd

model_name = 'model4.6.7'
semivariogram_sill = 0.0001946156

# Read in .csv with kriged predictions
df = pd.read_csv(f'/Users/tnye/kappa/krige/{model_name}/{model_name}_culled_krige_k0_linear.csv')

# Make file for k0 predictions
f = open(f'/Users/tnye/kappa/krige/{model_name}/{model_name}_culled_pred_k0.txt','w')
f.write('#Longitude,Latitude,kappa(s)\n')
for i in range(len(df)):
    f.write(f"{df['Longitude'].iloc[i]},{df['Latitude'].iloc[i]},{df['pred_k0'].iloc[i]}\n")
f.close()

# Make file for pred stddev (log 10)
f = open(f'/Users/tnye/kappa/krige/{model_name}/{model_name}_culled_pred_stddev.txt','w')
f.write('#Longitude,Latitude,Log10 k0 stddev\n')
for i in range(len(df)):
    f.write(f"{df['Longitude'].iloc[i]},{df['Latitude'].iloc[i]},{df['k0_stddev'].iloc[i]}\n")
f.close()


# Get regions where stddev is within sill
ind = np.where(df['k0_stddev']<=np.sqrt(semivariogram_sill))[0]
df_cut = df.loc[ind]

# Make file for k0 predictions
f = open(f'/Users/tnye/kappa/krige/{model_name}/{model_name}_culled_pred_k0_cut.txt','w')
f.write('#Longitude,Latitude,kappa(s)\n')
for i in range(len(df_cut)):
    f.write(f"{df_cut['Longitude'].iloc[i]},{df_cut['Latitude'].iloc[i]},{df_cut['pred_k0'].iloc[i]}\n")
f.close()

# Make file for pred stddev (log 10)
f = open(f'/Users/tnye/kappa/krige/{model_name}/{model_name}_culled_pred_stddev_cut.txt','w')
f.write('#Longitude,Latitude,Log10 k0 stddev\n')
for i in range(len(df_cut)):
    f.write(f"{df_cut['Longitude'].iloc[i]},{df_cut['Latitude'].iloc[i]},{df_cut['k0_stddev'].iloc[i]}\n")
f.close()
