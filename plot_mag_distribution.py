#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 21:11:56 2021

@author: tnye
"""

###############################################################################
# Makes historgrams of magnitude and fc distrivution for the Bay Area events.
###############################################################################

# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

############################# Plot Magnitudes #################################

# Read in events dataframe
events_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/SNR_5_file.csv')

# Select out unique events
events = np.unique(np.array(events_df['OrgT']))

# Get magnitude for events
mag = []
for event in events:
    m = events_df['Mag'][np.where(events_df==event)[0][0]]
    mag.append(m)

# Make histogram
plt.hist(mag, 15)
plt.xlabel('Magnitude')
plt.ylabel('Number of Events')
plt.title('Bay Area Events Magnitude Distribution')

plt.savefig('/Users/tnye/kappa/plots/misc/mag_dist.png', dpi=300)
plt.close()


############################## Calculate fc ###################################

B = 3100 #m/s
stressdrop = 5e6 #Pa

# If local magnitude is than 3, convert to moment magnitude
# Mag = []
# for m in mag:
#     if m < 3.0:
#         M = 0.884 + 0.754*m
#     else:
#         M = m
#     Mag.append(M)

Mag = []
for m in mag:
    if m >= 3.5:
        Mag.append(m)

# Convert Magnitude to corner frequency 
fc_list = []
for M in Mag:
    # Calcualte seismic moment
    M0 = 10.**((3./2.)*M + 9.1)
    # Calculate corner frequency
    fc = B*(stressdrop/(8.47*M0))**(1./3.)
    fc_list.append(fc)

# Make histogram
plt.hist(fc_list, 15)
plt.xlabel('Corner Frequency (Hz)')
plt.ylabel('Number of Events')
plt.title('Bay Area Events fc Distribution')

plt.savefig('/Users/tnye/kappa/plots/misc/fc_dist.png', dpi=300)
plt.close()
