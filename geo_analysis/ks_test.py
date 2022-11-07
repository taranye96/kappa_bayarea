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
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

model_name = 'model4.6.7'

fault = pd.read_csv(f'/Users/tnye/kappa/data/fault_zone_stns_{model_name}_culled.csv')
kappa_list = pd.read_csv(f'/Users/tnye/kappa/traditional_method/models/{model_name}/{model_name}_kappa_culled.out',delimiter=',')[' kappa(s) ']

kappa_within = kappa_list[np.where(np.array(fault['Fault_Zone']=='within'))[0]]
kappa_out = kappa_list[np.where(np.array(fault['Fault_Zone']=='outside'))[0]]

stat, pval = ks_2samp(kappa_out,kappa_within)


c_alpha = 1.36 # for significance level of 0.05
crit = c_alpha* np.sqrt((len(kappa_within)+len(kappa_out))/(len(kappa_within)*len(kappa_out)))

print(f'P-value: {round(pval,3)}')
print(f'Statistics: {round(stat, 3)}')
print(f'Critical value: {round(crit,3)}')

kappa = np.concatenate((kappa_out, kappa_within))

loc = np.array([])
loc = np.append(loc, np.full(len(kappa_out),'Outside'))
loc = np.append(loc, np.full(len(kappa_within),'Within'))

dataset_dict = {'$\kappa_0$ (s)':kappa,'Relation to Fault Zone':loc}
df = pd.DataFrame(data=dataset_dict)

fig, ax = plt.subplots(figsize=(8,6))
plt.grid(alpha=0.3)
sns.violinplot(x ='Relation to Fault Zone', y ='$\kappa_0$ (s)', data = df)
plt.savefig(f'/Users/tnye/kappa/plots/paper/{model_name}_fault_violin.png', dpi=300)


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'


fig, ax = plt.subplots(1,1, figsize=(8,6))
ax.hist(kappa_out, alpha=0.8, label='Outside', edgecolor='black', lw=0.5)
ax.hist(kappa_within, alpha=0.8, label='Within', edgecolor='black', lw=0.5)
plt.legend()
ax.set_xlabel('$\kappa_0$ (s)')
ax.set_ylabel('No. Stations')
ax.yaxis.grid()
ax.set_axisbelow(True)
plt.savefig(f'/Users/tnye/kappa/plots/paper/{model_name}_fault_hist.png', dpi=300)


