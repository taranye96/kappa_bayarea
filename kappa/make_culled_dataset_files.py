#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 22:41:14 2022

@author: tnye
"""

###############################################################################
# This script makes files using only the culled dataset. 
###############################################################################

# Imports
import numpy as np
import pandas as pd
import csv

model = 'model4.6.7'

############ Make kappa coordiantes csv using culled stations only ############

df = pd.read_csv(f'/Users/tnye/kappa/traditional_method/models/{model}/{model}_kappa.out',delimiter='\t')
xstns = pd.read_csv(f'/Users/tnye/kappa/traditional_method/models/{model}/no_krige_stns.csv')['Station'].values

rm_ind = []
for stn in xstns:
    try:
        i = np.where(df['#site ']==stn)[0][0]
        rm_ind.append(i)
    except:
        continue

df = df.drop(index=rm_ind).reset_index(drop=True)

df.to_csv(f'/Users/tnye/kappa/traditional_method/models/{model}/{model}_kappa_culled.csv')


############## Make kappa file for GMT using culled stations only #############

df = pd.read_csv(f'/Users/tnye/kappa/traditional_method/models/{model}/{model}_kappa.out',delimiter='\t')
GMT_kappa_df = pd.read_csv(f'/Users/tnye/kappa/GMT/data/kappa/{model}.txt',delimiter=',')
xstns = pd.read_csv(f'/Users/tnye/kappa/traditional_method/models/{model}/no_krige_stns.csv')['Station'].values

rm_ind = []
for stn in xstns:
    try:
        i = np.where(df['#site ']==stn)[0][0]
        rm_ind.append(i)
    except:
        continue

GMT_kappa_df = GMT_kappa_df.drop(index=rm_ind).reset_index(drop=True)

GMT_kappa_df.to_csv(f'/Users/tnye/kappa/GMT/data/kappa/{model}_culled.csv',index=False)


########## Make kappa stddev file for GMT using culled stations only ##########

df = pd.read_csv(f'/Users/tnye/kappa/traditional_method/models/{model}/{model}_kappa.out',delimiter='\t')
GMT_kappa_df = pd.read_csv(f'/Users/tnye/kappa/GMT/data/kappa/{model}_std.txt',delimiter=',')
xstns = pd.read_csv(f'/Users/tnye/kappa/traditional_method/models/{model}/no_krige_stns.csv')['Station'].values

rm_ind = []
for stn in xstns:
    try:
        i = np.where(df['#site ']==stn)[0][0]
        rm_ind.append(i)
    except:
        continue

GMT_kappa_df = GMT_kappa_df.drop(index=rm_ind).reset_index(drop=True)

GMT_kappa_df.to_csv(f'/Users/tnye/kappa/GMT/data/kappa/{model}_std_culled.csv',index=False)


