#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 00:23:50 2021

@author: tnye
"""

###############################################################################
# Script used to make a dataframe with PGA and PGV predicitions using the ASK14
# and BA14 GMMs. Takes as an input the gmm_init.csv dataframe from
# initialize_gmm_df.py. 
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
import tsueqs_main_fns as tmf
from rotd50 import compute_rotd50

model_name = 'model4.6.7'

df = pd.read_csv(f'/Users/tnye/kappa/data/flatfiles/gmm_init_{model_name}_culled.csv')

wf_dir = '/Users/tnye/kappa/data/waveforms/bay_acc_2000-2021/instrument_corrected/cut_SNR_3_80%'

obs_pga = []
obs_pgv = []
ask14_pga = []
ask14_pga_std = []
ask14_pga_Vs30ref = []
ask14_pga_Vs30ref_std = []
ask14_pgv = []
ask14_pgv_std = []
ask14_pgv_Vs30ref = []
ask14_pgv_Vs30ref_std = []
boore_pga = []
boore_pga_std = []
boore_pga_Vs30ref = []
boore_pga_Vs30ref_std = []
boore_pgv = []
boore_pgv_std = []
boore_pgv_Vs30ref = []
boore_pgv_Vs30ref_std = []
dists = []
mags = []

ex_ind = []

# Get observed PGA and PGV
for i in range(len(df)):
    
    stn = df['Station Name'].iloc[i]
    org = df['OrgT'].iloc[i]
    yyyy,mth,dd = org.split(' ')[0].split('-')
    hh,mm,ss = org.split(' ')[1].split(':')
    event = f'Event_{yyyy}_{mth}_{dd}_{hh}_{mm}_{ss}'
    
    records = glob(f'{wf_dir}/{event}/*_{stn}_H*')
    
    st_E = read(records[0])
    st_N = read(records[1])
    
    
    if len(st_E[0].data) != len(st_N[0].data):
        
        # Check for gaps:
        st = st_E.append(st_N[0])
        gaps = st.get_gaps()
        
        if gaps == []:
            # Minimum record length between the components
            npts = np.min([len(st_E[0].data),len(st_N[0].data)])
            
            # Shorten the components to be the same length 
            st_N[0].data = st_N[0].data[:npts]
            st_E[0].data = st_E[0].data[:npts]
        
    
    if len(st_E[0].data) == len(st_N[0].data):
        # acc_euc_norm = avg.get_eucl_norm_2comp(read(records[0])[0].data,read(records[1])[0].data)
        rotd50_acc = compute_rotd50(st_E[0].data,st_N[0].data)
        obs_pga.append(np.max(np.abs(rotd50_acc)))
    
        # Convert to velocity 
        E_vel_unfilt = tmf.accel_to_veloc(st_E)
        N_vel_unfilt = tmf.accel_to_veloc(st_N)
        
        # Highpass filter vleocity
        stsamprate = read(records[0])[0].stats.sampling_rate
        fcorner = 1/15.                        
        order = 2
        
        E_vel = tmf.highpass(E_vel_unfilt,fcorner,stsamprate,order,zerophase=True)
        N_vel = tmf.highpass(N_vel_unfilt,fcorner,stsamprate,order,zerophase=True)
    
        # Get PGV residual 
        rotd50_vel = compute_rotd50(E_vel[0].data,N_vel[0].data)
        obs_pgv.append(np.max(np.abs(rotd50_vel)))

        # Get parameters
        M = df['Magnitude'].iloc[i]
        Rrup = df['Rrup(km)'].iloc[i]
        vs30 = df['Vs30(m/s)'].iloc[i]
        vs30_meas = df['Vs30_meas(m/s)'].iloc[i]
        if vs30_meas == 0.0:
            vs30_meas = False
        else:
            vs30_meas = True
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
        pga, pga_std = gmm.ask14('PGA',M,Rrup,vs30,vs30_meas,ztor,rake,dip,width,z1pt0)
        pga = np.exp(pga)*sp.g #ln %g to m/s/s
        ask14_pga.append(np.abs(pga[0]))
        ask14_pga_std.append(pga_std[0][0])
        
        pga, pga_std = gmm.ask14('PGA',M,Rrup,760.0,False,ztor,rake,dip,width,z1pt0)
        pga = np.exp(pga)*sp.g #ln %g to m/s/s
        ask14_pga_Vs30ref.append(np.abs(pga[0]))
        ask14_pga_Vs30ref_std.append(pga_std[0][0])
        
        pgv, pgv_std = gmm.ask14('PGV',M,Rrup,vs30,vs30_meas,ztor,rake,dip,width,z1pt0)
        pgv = np.exp(pgv)*sp.g #ln %g to m/s/s
        ask14_pgv.append(np.abs(pgv[0]))
        ask14_pgv_std.append(pgv_std[0][0])
        
        pgv, pgv_std = gmm.ask14('PGV',M,Rrup,760.0,False,ztor,rake,dip,width,z1pt0)
        pgv = np.exp(pgv)*sp.g #ln %g to m/s/s
        ask14_pgv_Vs30ref.append(np.abs(pgv[0]))
        ask14_pgv_Vs30ref_std.append(pgv_std[0][0])
        
        # Call Boore14 GMM
        pga, pga_std = gmm.boore2014('PGA',M,Rrup,vs30,vs30_meas,ztor=ztor,rake=rake)
        pga = np.exp(pga)*sp.g #ln %g to m/s/s
        boore_pga.append(np.abs(pga[0]))
        boore_pga_std.append(pga_std[0][0])
        
        pga, pga_std = gmm.boore2014('PGA',M,Rrup,760.0,False,ztor,rake)
        pga = np.exp(pga)*sp.g #ln %g to m/s/s
        boore_pga_Vs30ref.append(np.abs(pga[0]))
        boore_pga_Vs30ref_std.append(pga_std[0][0])
        
        pgv, pgv_std = gmm.boore2014('PGV',M,Rrup,vs30,vs30_meas,ztor,rake)
        pgv = np.exp(pgv)*sp.g #ln %g to m/s/s
        boore_pgv.append(np.abs(pgv[0]))
        boore_pgv_std.append(pgv_std[0][0])
        
        pgv, pgv_std = gmm.boore2014('PGV',M,Rrup,760.0,False,ztor,rake)
        pgv = np.exp(pgv)*sp.g #ln %g to m/s/s
        boore_pgv_Vs30ref.append(np.abs(pgv[0]))
        boore_pgv_Vs30ref_std.append(pgv_std[0][0])
        
        dists.append(Rrup)
        mags.append(M)

    else:
        ex_ind.append(i)

df = df.drop(index=ex_ind)

# Get PGA and PGV residual
ask14_pga_res = np.log(np.array(obs_pga)) - np.log(np.array(ask14_pga))
ask14_Vs30ref_pga_res = np.log(np.array(obs_pga)) - np.log(np.array(ask14_pga_Vs30ref))
ask14_pgv_res = np.log(np.array(obs_pgv)) - np.log(np.array(ask14_pgv))
ask14_Vs30ref_pgv_res = np.log(np.array(obs_pgv)) - np.log(np.array(ask14_pgv_Vs30ref))

boore_pga_res = np.log(np.array(obs_pga)) - np.log(np.array(boore_pga))
boore_Vs30ref_pga_res = np.log(np.array(obs_pga)) - np.log(np.array(boore_pga_Vs30ref))
boore_pgv_res = np.log(np.array(obs_pgv)) - np.log(np.array(boore_pgv))
boore_Vs30ref_pgv_res = np.log(np.array(obs_pgv)) - np.log(np.array(boore_pgv_Vs30ref))
    
df['Obs_PGA(m/s/s)'] = obs_pga
df['Obs_PGV(m/s)'] = obs_pgv
df['ASK14_PGA(m/s/s)'] = ask14_pga
df['ASK14_PGA_std'] = ask14_pga_std
df['lnASK14_PGA_Res'] = ask14_pga_res
df['ASK14_PGA_Vs30ref(m/s/s)'] = ask14_pga_Vs30ref
df['ASK14_PGA_Vs30ref_std'] = ask14_pga_Vs30ref_std
df['lnASK14_PGA_Vs30ref_Res'] = ask14_Vs30ref_pga_res
df['ASK14_PGV(m/s)'] = ask14_pgv
df['ASK14_PGV_std'] = ask14_pgv_std
df['lnASK14_PGV_Res'] = ask14_pgv_res
df['ASK14_PGV_Vs30ref(m/s)'] = ask14_pgv_Vs30ref
df['ASK14_PGV_Vs30ref_std'] = ask14_pgv_Vs30ref_std
df['lnASK14_PGV_Vs30ref_Res'] = ask14_Vs30ref_pgv_res
df['Boore14_PGA(m/s/s)'] = boore_pga
df['Boore14_PGA_std'] = boore_pga_std
df['lnBoore14_PGA_Res'] = boore_pga_res
df['Boore14_PGA(m/s/s)'] = boore_pga_Vs30ref
df['Boore14_PGA_std'] = boore_pga_Vs30ref_std
df['lnBoore14_PGA_Vs30ref_Res'] = boore_Vs30ref_pga_res
df['Boore14_PGV(m/s)'] = boore_pgv
df['Boore14_PGV_std'] = boore_pgv_std
df['lnBoore14_PGV_Res'] = boore_pgv_res
df['Boore14_PGV_Vs30ref(m/s)'] = boore_pgv_Vs30ref
df['Boore14_PGV_Vs30ref_std'] = boore_pgv_Vs30ref_std
df['lnBoore14_PGV_Vs30ref_Res'] = boore_Vs30ref_pgv_res

df.to_csv(f'/Users/tnye/kappa/data/flatfiles/GMM_{model_name}_culled.csv')


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
ax.scatter(dists[i],boore_pga[i],s=5,c=colorVal,alpha=0.5,lw=.8,label='Boore14')
ax.scatter(dists,obs_pga,marker='x',s=12,color=colors,alpha=0.7,lw=.8,label='Observed')
ax.scatter(dists,ask14_pga,marker='.',s=12,color=colors,alpha=0.7,lw=.8,label='ASK14')
ax.set_xlabel('Rrup(km)')
ax.set_ylabel('PGA(m/s/s)')
ax.set_yscale('log')
ax.set_xscale('log')
plt.legend()
plt.savefig(f'/Users/tnye/kappa/plots/gmm/obs_v_gmm_{model_name}_culled.png', dpi=300)