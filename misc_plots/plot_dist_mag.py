#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 12:42:31 2021

@author: tnye
"""

###############################################################################
# Script that plots hypocentral distance vs earthquake magnitude for all the
# events and broadband stations. 
###############################################################################

# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

catalog = pd.read_csv('/Users/tnye/kappa/data/flatfiles/SNR_5_file.csv')

magnitude = np.array(catalog['Mag'])
rhyp = np.array(catalog['rhyp'])
depth = np.array(catalog['Qdep']/1000)

plt.figure()
sc = plt.scatter(rhyp, magnitude, c=depth, cmap='plasma_r', s=2)
# plt.ylim(2.0,6.0)
plt.grid(alpha=0.5)
cbar = plt.colorbar(sc, fraction=0.03, pad=0.01)
cbar.ax.invert_yaxis()
cbar.set_label('Depth (km)', rotation=270, labelpad=10)
plt.xlabel('Recrod Rhyp (km)')
plt.ylabel('Magnitdue')
plt.savefig('/Users/tnye/kappa/plots/stats/rhyp_mag.png', dpi=300)