#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 14:50:41 2021

@author: tnye
"""

###############################################################################
# Script that plots histogram of time to reach 80% Arias Intensity. 
###############################################################################

# Imports
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt

energy_data = genfromtxt('/Users/tnye/kappa/data/energy.log',delimiter=",")[:,2]
# energy_data = genfromtxt('/Users/tnye/kappa/data/energy.log',dtype=None,delimiter=",")
# energy = energy_data[:,2]
# events = energy_data[:,0]
# stations = energy_data[:,1]


plt.figure()
plt.hist(energy_data,10)
plt.xlabel('Time (s) from S-arrival to reach 80% energy release')
plt.ylabel('Counts')
plt.show()
plt.savefig('/Users/tnye/kappa/plots/misc/acc_energy_hist.png', dpi=300)