#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 00:58:34 2022

@author: tnye
"""

###############################################################################
# Script that runs the kappa calculations.
# Outputs are:
    # model.log with run information 
    # event_density.log with number of events recorded at each station
    # kappa.out file with k0 estiamtes for each station
    # .csv files for each station with record kappa estimations 
    # plots of the station records with frequency range lines
    # plots of the station kappa regression
    # GMT text file with station k0 values
###############################################################################

# Imports 
import numpy as np
import pandas as pd
from glob import glob
from os import path, makedirs
from obspy import read
from mtspec import mtspec
from plot_stn_spectra import plot_spectra
from plot_kappa_regression import plot_kappa
from make_gmt_kappa import make_gmt
import tsueqs_main_fns as tmf
import kappa_methods

import warnings
warnings.filterwarnings("ignore")

################################ Parameters ###################################

model_name='model4.6.7'
home_dir = f'/Users/tnye/kappa/traditional_method/models/{model_name}'

wf_dir = '/Users/tnye/kappa/data/waveforms/bay_acc_2000-2021/instrument_corrected/cut_SNR_3_80%'
mseeds = sorted(glob(f'{wf_dir}/*/*'))

outfilename = f'{home_dir}/{model_name}_kappa.out'

startfreq = 'min of 1st derivative'
endfreq = '-8e-5 Diff1 threshold'
deg = 15              # polyfit degree
max_dist = 150        # max epicentral distance in km
min_mag = 3.5         # min magnitude
max_mag = 5.5         # max magnitude
delta_freq = 10       # min frequency ranged needed to solve for k0
min_events = 8        # min number of events recorded at a station to solve for k0

freq_model = 'm4'     # desired model in inv_fns for selecting fx

bb = 0
sm = 0

# Make directories 
if not path.exists(f'{home_dir}'):
    makedirs(f'{home_dir}')
        
if not path.exists(f'{home_dir}/stn_flatfiles'):
    makedirs(f'{home_dir}/stn_flatfiles')
    
# Make log file with number of events recorded on each station
f_event_density = open(f'{home_dir}/{model_name}.log','w')

# Make log file with criteria
f_criteria = open(f'{home_dir}/{model_name}.log','w')
f_criteria.write('L2_Norm kappa with 25% diff\n')
f_criteria.write('Geometric mean HVSR\n')
f_criteria.write(f'Minimum magnitude: {min_mag}\n')
f_criteria.write(f'Maximum Repi (km): {max_dist}\n')
f_criteria.write(f'Minimum events recorded per station: {min_events}\n')
f_criteria.write(f'Frequency range model: {freq_model}\n')
f_criteria.write(f'Min df (Hz): {delta_freq} (Hz)\n')
f_criteria.write(f'Frequency range (Hz): {startfreq}-{endfreq}\n')
f_criteria.write(f'Polyfit degree: {deg}\n')

# Read in dataframe 
event_file = '/Users/tnye/kappa/data/flatfiles/data_download/bay_2000-2021_SNR_irfilt.csv'


######################### Run Kappa Calculations ##############################

# Read in event df and make sure all names are strings
lst_str_cols = ['Name','Channel']
dict_dtypes = {x : 'str'  for x in lst_str_cols}
event_df = pd.read_csv(event_file, dtype=dict_dtypes)

# Get list of stations
stations = []
for file in mseeds:
    stn = file.split('/')[-1].split('_')[1]
    stations.append(stn)
stations = np.unique(stations)

n_stations = 0
n_records = 0

stns = []
stn_kappa_list = []
stn_kappa_std_list = []
r2_list = []
mag_list = np.array([])
depth_list = np.array([])
# rhyp_list = np.array([])
repi_list = np.array([])
stn_list = np.array([])
final_events = np.array([])

events = []

# Loop through stations
for stn in stations:
    
    stn_events = []
    stn_files = glob(f'{wf_dir}/*/*{stn}_H*E_*')
    for file in stn_files:
        ev = file.split('/')[-2]
        yyyy,mth,dd,hh,mm,ss = ev.split('_')[1:]
        stn_events.append(f'{yyyy}-{mth}-{dd} {hh}:{mm}:{ss}')

    # Select out only data from event_df for current station
    stn_df0 = event_df[(event_df['Name'] == stn)].reset_index(drop=True)
    stn_df0 = stn_df0.sort_values(by=['OrgT'],ascending=False)
    
    # Remove events from station df that we do not have records for
    stn_df1 = stn_df0.copy().reset_index(drop=True)
    ex_ind = []
    for i in range(len(stn_df0)):
        org = stn_df0['OrgT'].iloc[i]
        if org not in stn_events:
            ex_ind.append(i)
    stn_df1 = stn_df1.drop(index=ex_ind).reset_index(drop=True)
            
    # Remove events from station df that do not fit the magnitude and dist criteria
    ex = []
    ex_ind = []
    stn_df2 = stn_df1.copy()
    stn_df2['repi'] = ""
    for i in range(len(stn_df1)):
        Slon = stn_df1['Slon'].iloc[i]
        Slat = stn_df1['Slat'].iloc[i]
        Qlon = stn_df1['Qlon'].iloc[i]
        Qlat = stn_df1['Qlat'].iloc[i]
        mag = stn_df1['Mag'].iloc[i]
        repi = tmf.compute_repi(Slon,Slat,Qlon,Qlat)
        stn_df2['repi'].iloc[i] = repi
        if repi > max_dist or mag < min_mag or mag > max_mag:
            ex_ind.append(i)
    stn_df2 = stn_df2.drop(index=ex_ind).reset_index(drop=True)
    
    f_event_density.write(f'{stn}: {len(stn_df2)} events\n')
    
    ## Our event list should now match the event dataframe

    ################################# Tstar (Kappa) ###############################
    
    # Number of events
    I = len(stn_df2)
    
    if I == 0:
        print(f'{stn} recorded no events')
    else:
        print(f'{stn} recorded {I} events')
    
        # Initialize lists for kappa 
        fe_list = []
        fx_list = []
        fc_list = []
        freq_list = []
        amp_list = []
        kappa_list = []
        A0_list = []
        L2norm = []
        kappa_std = []
        A_std = []
        
        temp_stns = []
        temp_mag = []
        temp_repi = []
        temp_events = []
        temp_depth = []
        
        ex_ind = []
        # Loop through events
        for i in range(I):
            
            ############################ Set up inversion #############################
            
            event = stn_df2['OrgT'].iloc[i]
            yyyy,mth,dd = event.split(' ')[0].split('-')
            hh,mm,ss = event.split(' ')[1].split(':')
            M = stn_df2['Mag'].iloc[i]
            
            # Gather horizontal component waveforms for this station/event pair
            N_record = glob(f'{wf_dir}/Event_{yyyy}_{mth}_{dd}_{hh}_{mm}_{ss}/*{stn}_H*N_*')[0]
            E_record = glob(f'{wf_dir}/Event_{yyyy}_{mth}_{dd}_{hh}_{mm}_{ss}/*{stn}_H*E_*')[0]
            
            # Read in records
            st_N = read(N_record)
            st_E = read(E_record)
            
            # Get power spectra
            if st_N[0].stats.npts == st_E[0].stats.npts:
                
                # Get stats
                dt = st_N[0].stats.delta
                samprate = st_N[0].stats.sampling_rate
                nfft = st_N[0].stats.npts
                
                # Compute power spectra
                N_amp2, freq, N_jack, N_fstat, N_dof =  mtspec(st_N[0].data*100,delta=dt,nfft=nfft,time_bandwidth=4,number_of_tapers=7,quadratic=True,statistics=True)
                E_amp2, freq, E_jack, E_fstat, E_dof =  mtspec(st_E[0].data*100,delta=dt,nfft=nfft,time_bandwidth=4,number_of_tapers=7,quadratic=True,statistics=True)
                
                # Get standard deviation of power spectra
                sigmaN2 = (N_jack[:,1] - N_jack[:,0])/3.29
                sigmaE2 = (E_jack[:,1] - E_jack[:,0])/3.29
                
                try:
                    fe,fx,fc,data_NE,kappa,k_std,A0,A0_std = kappa_methods.L2_norm_kappa(stn, M, N_amp2, E_amp2, sigmaN2, sigmaE2, freq, freq_model, samprate, deg, delta_freq)
                
                    fe_list.append(fe)
                    fx_list.append(fx)
                    fc_list.append(fc)
                    freq_list.append(freq)
                    amp_list.append(data_NE)
                    
                    # Append tstar(kappa) and A0 to lists
                    kappa_list.append(float(kappa))
                    A0_list.append(A0)
                    
                    # Append standard deviations 
                    kappa_std.append(k_std)#error in kappa
                    A_std.append(A0_std)
                    
                    temp_stns.append(stn)
                    temp_mag.append(M)
                    temp_repi.append(repi)
                    temp_depth.append(stn_df2['Qdep'].iloc[i])
                    temp_events.append(event)
                except:
                    ex_ind.append(i)

            else:
                ex_ind.append(i)
        
        stn_df2 = stn_df2.drop(index=ex_ind)
        stn_df2 = stn_df2.reset_index(drop=True)
        
        # Check how many records are left now
        if len(kappa_list) >= min_events:
            
            stn_df2['fc'] = fc_list
            stn_df2['fe'] = fe_list
            stn_df2['fx'] = fx_list
            stn_df2['A0'] = A0_list
            stn_df2['A0 std'] = A_std
            stn_df2['Kappa(s)'] = kappa_list
            stn_df2['Kappa std dev'] = kappa_std
            
            stn_df2.to_csv(f'{home_dir}/stn_flatfiles/{stn}.csv')
            
            n_stations += 1
            
            if stn_df2['Channel'].iloc[0] == 'HH':
                bb+=1
            elif stn_df2['Channel'].iloc[0] == 'HN':
                sm+=1
        
            # Plot spectra
            if not path.exists(f'{home_dir}/plots'):
                makedirs(f'{home_dir}/plots')
            figpath = f'{home_dir}/plots/{stn}_record_spectra.png'
            
            # Plot spectra with frequency ranges indicated
            plot_spectra(freq_list,amp_list,fe_list,fx_list,mag,kappa_list,A0_list,stn_df2,st_N[0].stats.sampling_rate,figpath)
        
            # Plot kappa
            kappa,kappa_std,r_squared = plot_kappa(model_name,home_dir,stn_df2,min_events)
            stn_kappa_list.append(kappa)
            stn_kappa_std_list.append(kappa_std)
            r2_list.append('')
            stns.append(stn)
            
            n_records += I
            
            stn_list = np.append(stn_list, temp_stns)
            mag_list = np.append(mag_list, temp_mag)
            repi_list = np.append(repi_list, temp_repi)
            depth_list = np.append(depth_list, temp_depth)
            final_events = np.append(final_events, temp_events)

data = {'Station':stn_list,'Event':final_events,'Mag':mag_list,'Depth':depth_list,'Repi':repi_list}
final_df = pd.DataFrame(data)
final_df.to_csv(f'/Users/tnye/kappa/data/flatfiles/final_dataset_{model_name}.csv')
           
f_criteria.write(f'Average Polyfit R2: {np.mean(r2_list)}\n')
f_criteria.write('\n')
f_criteria.write(f'Total Number of Events: {len(np.unique(final_events))}\n')
f_criteria.write(f'Total Number of Useable Stations: {n_stations}\n')
f_criteria.write(f'Total Broadband: {bb}\n')
f_criteria.write(f'Total Strong Motion of Useable Stations: {sm}\n')
f_criteria.write(f'Total Records: {n_records}')
    
# Save kappa to text file
outfile = open(outfilename, 'w')
out = np.column_stack((list(stns), stn_kappa_list, stn_kappa_std_list, r2_list))

# Write value to file and save
outfile.write('#site \t kappa(s) \t kappa_std \t R-squared\n')
np.savetxt(outfile, out, fmt='%s', delimiter='\t')
outfile.close()
    
f_criteria.close()
f_event_density.close()

# Make file for plotting in GMT
make_gmt(home_dir,model_name,outfilename,event_df)
        
