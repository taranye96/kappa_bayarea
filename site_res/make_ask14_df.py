#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 00:23:50 2021

@author: tnye
"""

###############################################################################
# Script used to make a dataframe with PGA and PGV predicitions using the ASK14
# GMM. Takes as an input the gmm_init.csv dataframe. 
###############################################################################

# Imports
import numpy as np
import scipy.constants as sp
import pandas as pd
from obspy import read
from glob import glob
import sys
sys.path.append("/Users/tnye/kappa/code/")
import gmm_call as gmm
import signal_average_fns as avg
import tsueqs_main_fns as tmf

model_name = 'model1'

df = pd.read_csv(f'/Users/tnye/kappa/data/flatfiles/gmm_init_{model_name}.csv')

obs_pga = []
obs_pgv = []
ask14_pga = []
ask14_pga_std = []
ask14_pgv = []
ask14_pgv_std = []
boore_pga = []
boore_pgv_std = []
dists = []
mags = []

# Get observed PGA and PGV
for i in range(len(df)):
    
    stn = df['Station Name'].iloc[i]
    org = df['OrgT'].iloc[i]
    yyyy,mth,dd = org.split(' ')[0].split('-')
    hh,mm,ss = org.split(' ')[1].split(':')
    event = f'Event_{yyyy}_{mth}_{dd}_{hh}_{mm}_{ss}'
    
    records = glob(f'/Users/tnye/kappa/data/waveforms/acc_all/Vs3.1/cut2_20/{event}/*_{stn}*')
    
    acc_euc_norm = avg.get_eucl_norm_2comp(read(records[0])[0].data,read(records[1])[0].data)
    obs_pga.append(np.max(np.abs(acc_euc_norm)))
    
    # Convert to velocity 
    E_vel_unfilt = tmf.accel_to_veloc(read(records[0]))
    N_vel_unfilt = tmf.accel_to_veloc(read(records[1]))
    
    # Highpass filter vleocity
    stsamprate = read(records[0])[0].stats.sampling_rate
    fcorner = 1/15.                        
    order = 2
    
    E_vel = tmf.highpass(E_vel_unfilt,fcorner,stsamprate,order,zerophase=True)
    N_vel = tmf.highpass(N_vel_unfilt,fcorner,stsamprate,order,zerophase=True)
    
    # Get PGV residual 
    vel_euc_norm = avg.get_eucl_norm_2comp(E_vel[0].data,N_vel[0].data)
    obs_pgv.append(np.max(np.abs(vel_euc_norm)))
 
# Get ASK14 PGA
for i in range(len(df)):
    
    # Get parameters
    M = df['Magnitude'].iloc[i]
    Rrup = df['Rrup(km)'].iloc[i]
    vs30 = df['Vs30(m/s)'].iloc[i]
    ztor = df['Ztor(km)'].iloc[i]
    if pd.notna(df['Rake'].iloc[i]):
        rake = df['Rake'].iloc[i]
    else:
        rake = 0.0
    if pd.notna(df['Dip'].iloc[i]):
        dip = df['Dip'].iloc[i]
    else:
        dip = 90.0
    if pd.notna(df['Width(km)'].iloc[i]):
        width = df['Width(km)'].iloc[i]
    else:
        width = 10.0
    z1pt0 = df['Z1pt0(m)'].iloc[i]
    
    # Call ASK14 GMM
    pga, pga_std = gmm.ask14('PGA',M,Rrup,vs30,ztor,rake,dip,width,z1pt0)
    pga = np.exp(pga)*sp.g #ln %g to m/s/s
    ask14_pga.append(np.abs(pga[0]))
    ask14_pga_std.append(pga_std[0][0])
    
    pgv, pgv_std = gmm.ask14('PGV',M,Rrup,vs30,ztor,rake,dip,width,z1pt0)
    pgv = np.exp(pgv)*sp.g #ln %g to m/s/s
    ask14_pgv.append(np.abs(pgv[0]))
    ask14_pgv_std.append(pgv_std[0][0])
    
    dists.append(Rrup)
    mags.append(M)

# Get PGA and PGV residual
pga_res = np.log(np.array(obs_pga)) - np.log(np.array(ask14_pga))
pgv_res = np.log(np.array(obs_pgv)) - np.log(np.array(ask14_pgv))
    
df['Obs_PGA(m/s/s)'] = obs_pga
df['Obs_PGV(m/s)'] = obs_pgv
df['ASK14_PGA(m/s/s)'] = ask14_pga
df['ASK14_PGA_std'] = ask14_pga_std
df['lnASK14_PGA_Res'] = pga_res
df['ASK14_PGV(m/s)'] = ask14_pgv
df['ASK14_PGV_std'] = ask14_pgv_std
df['lnASK14_PGV_Res'] = pgv_res
# df['Boore14_PGA'] = boore_pga
# df['Boore14_PGA_std'] = boore_std
# df['Boore14_Residual'] = boore_res

df.to_csv(f'/Users/tnye/kappa/data/flatfiles/ASK14_pga_{model_name}.csv')


# Plot
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

# Set up colormap
viridis = plt.get_cmap('viridis_r') 
cNorm  = colors.Normalize(vmin=3.5, vmax=5.5)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=viridis)
colors = []
for i in range(len(mags)):
    colors.append(scalarMap.to_rgba(mags[i]))

fig = plt.figure(figsize=(10,8))
ax = plt.gca()
fig.colorbar(scalarMap, label='Magnitude', ax=ax)
# ax.scatter(dists[i],boore_pga[i],s=5,c=colorVal,alpha=0.5,lw=.8,label='Boore14')
ax.scatter(dists,obs_pga,marker='x',s=12,color=colors,alpha=0.7,lw=.8,label='Observed')
ax.scatter(dists,ask14_pga,marker='.',s=12,color=colors,alpha=0.7,lw=.8,label='ASK14')
ax.set_xlabel('Rrup(km)')
ax.set_ylabel('PGA(m/s/s)')
ax.set_yscale('log')
ax.set_xscale('log')
plt.legend()
plt.savefig(f'/Users/tnye/kappa/plots/gmm/obs_v_ask14_{model_name}.png', dpi=300)