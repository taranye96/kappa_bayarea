#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 16:26:38 2021

@author: tnye
"""

###############################################################################
# Script that plots the overal pga residuals using the ASK14 gmm. 
###############################################################################

# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/ASK14_pga.csv')

obs = np.array(df['Obs_PGA(m/s/s)'])
ask_pga = np.array(df['ASK14_PGA(m/s/s)'])
res = np.array(df['ASK14_Residual(m/s/s)'])
logres = np.log10(obs) - np.log10(ask_pga)
rrup = np.array(df['Rrup(km)'])
logrrup = np.log10(rrup)

x = np.linspace(0, 250)
y = np.zeros(len(x))
plt.figure()
plt.scatter(rrup, res)
plt.plot(x,y,ls='--',c='black')
plt.xlabel('Rrup(km)')
plt.ylabel('PGA Residual (m/s/s)')
plt.xlim(0,250)
plt.title('Abrahmson et al. 2014 GMPE Residuals')
# plt.show()
plt.savefig('/Users/tnye/kappa/plots/gmm/gmpe_residuals2.png', dpi=300)
plt.close

x = np.linspace(0, 2.5)
y = np.zeros(len(x))
plt.scatter(logrrup, logres)
plt.plot(x,y,ls='--',c='black')
plt.xlabel('Log Rrup(km)')
plt.ylabel('Log PGA Residuals (m/s/s)')
plt.xlim(0,2.5)
