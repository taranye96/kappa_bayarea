#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 17:51:00 2021

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
from numpy.polynomial import polynomial as P
from scipy.signal import argrelextrema
from sklearn.metrics import r2_score
import inv_fns as freqrange
from kappa_inversion import run_kappa_inv
from plot_stn_spectra import plot_spectra
from plot_kappa_regression import plot_kappa
from make_gmt_kappa import make_gmt

################################ Parameters ###################################

model_name='model1'
home_dir = f'/Users/tnye/kappa/traditional_method/models/{model_name}'

spectra_dir = '/Users/tnye/kappa/data/spectra/acc'

outfilename = f'{home_dir}/{model_name}_kappa.out'

startfreq = 'min of 1st derivative'
endfreq = '2e-6 threshold'
# endfreq = 'Zero 2nd derivative'
max_dist = 250
min_mag = 3.5
max_mag = 5.5

min_events = 10

freq_model = 'm1'

bb = 0
sm = 0

# Make directories 
if not path.exists(f'{home_dir}'):
    makedirs(f'{home_dir}')
        
if not path.exists(f'{home_dir}/stn_flatfiles'):
    makedirs(f'{home_dir}/stn_flatfiles')

# Make log file with criteria
f_criteria = open(f'{home_dir}/{model_name}.log','w')
f_criteria.write('Including stations BJOB and N001')
f_criteria.write(f'Minimum magnitude: {min_mag}\n')
f_criteria.write(f'Maximum distance (km): {max_dist}\n')
f_criteria.write(f'Minimum events recorded per station: {min_events}\n')
f_criteria.write(f'Frequency range (Hz): {startfreq}-{endfreq}\n')

# Read in dataframe 
event_file = '/Users/tnye/kappa/data/flatfiles/full_bay_flatfile.csv'

r2_list = []


######################### Run Kappa Calcualtions ##############################

# Read in event df and make sure all names are strings
lst_str_cols = ['Name','Channel']
dict_dtypes = {x : 'str'  for x in lst_str_cols}
event_df = pd.read_csv(event_file, dtype=dict_dtypes)

# Get list of stations
stations = np.unique(np.array(event_df['Name']))

# Stations Valerie and I decided to remove
# stations = np.delete(stations, np.where(stations=='BJOB'))
# stations = np.delete(stations, np.where(stations=='N001'))

# Start log file for number of events record on station
f_event_density = open(f'{home_dir}/event_density.log','w')

n_stations = 0
n_records = 0

stns = []
stn_kappa_list = []
stn_kappa_std_list = []
r2_list = []
mag_list = np.array([])
depth_list = np.array([])
rhyp_list = np.array([])
rrup_list = np.array([])
stn_list = np.array([])
final_events = np.array([])

events = []

# Loop through stations
for stn in stations:
    
    temp_stns = []
    temp_mag = []
    temp_rhyp = []
    temp_rrup = []
    temp_events = []
    temp_depth = []

    # Select out only data from event_df for current station
    stn_df0 = event_df[(event_df['Name'] == stn)].reset_index(drop=True)
    stn_df0 = stn_df0.sort_values(by=['OrgT'],ascending=False)
    
    # Get list of event directories for which we have spectra 
    event_list = glob(f'{spectra_dir}/*/*_{stn}_*')
    event_list = sorted(event_list,reverse=True)
    
    # Make a list of the event origin times using event directory names
    record_origins = []
    for i in range(len(event_list)):
        yyyy,mth,dd,hr,mm,sec = event_list[i].split('/')[-1].split('_')[4:10]
        sec = sec.split('.')[0]
        org = f'{yyyy}-{mth}-{dd} {hr}:{mm}:{sec}' 
        record_origins.append(org)
    
    # Remove duplicate origin times for separate events
        # I'm assuming the 2nd is the one we have data for since it likely 
        # overwrote the first one
    stn_df1 = stn_df0.copy().reset_index(drop=True)
    ex_ind = [idx for idx, val in enumerate(np.array(stn_df0['OrgT'])) if val in np.array(stn_df0['OrgT'])[:idx]]
    stn_df1 = stn_df1.drop(ex_ind)
    
    # Remove events from event files that are not in the dataframe
        # Perhaps removing strong motion data for station that I also have broadband data)
    missing = []
    missing_ind = []
    for i in range(len(record_origins)):
        if record_origins[i] not in np.array(stn_df1['OrgT']):
            missing.append(record_origins[i])
            missing_ind.append(i)
    for index in sorted(missing_ind, reverse=True):
        del event_list[index]
    
    # Remove events from event df that are not in the event files
    missing = []
    missing_ind = []
    stn_df2 = stn_df1.copy().reset_index()
    for i in range(len(stn_df1)):
        if stn_df1['OrgT'].iloc[i] not in record_origins:
            a=stn_df1['OrgT'].iloc[i]
            missing.append(stn_df1['OrgT'].iloc[i])
            missing_ind.append(i)
            stn_df2 = stn_df2.drop([i])
    
   
    # Remove events from event df that do not fit the magnitude and dist criteria
    ex = []
    ex_ind = []
    stn_df3 = stn_df2.copy()
    for i in range(len(stn_df2)):
        rhyp = stn_df2['rhyp'].iloc[i]
        mag = stn_df2['Mag'].iloc[i]
        org = stn_df2['OrgT'].iloc[i]
        depth = stn_df2['Qdep'].iloc[i]/1000
        rrup = np.tan(np.arccos(depth/rhyp))*depth
        if rhyp >= max_dist or mag <= min_mag or mag >= max_mag:
            ex_ind.append(i)
            ex.append(org)
            stn_df3 = stn_df3.drop(index=stn_df3.index[np.where(stn_df3['OrgT']==org)[0][0]])
    
    # Remove events from event list that do not fit the magnitude and dist criteria
    for index in sorted(ex_ind, reverse=True):
        del event_list[index]
    
    f_event_density.write(f'{stn}: {len(event_list)} events\n')
    
    ## Our event list should now match the event dataframe
    
    # Exclude these stations from kappa inversion because of site response
        # But keep for dataset stats
    if (stn != 'BJOB') & (stn != 'N001'):

        ################################# Tstar (Kappa) ###############################
        
        # Number of events
        I = len(event_list)
        print(f'{stn} recorded {I} events')
        
        if I == 0:
            print(f'{stn} recorded no events')
        else:
                
            # Number of frequency bins
            J = 150 
        
            # Initialize lists for kappa (tstar)
            fe_list = []
            fx_list = []
            fc_list = []
            kappa_list = []
            A0_list = []
            L2norm = []
            kappa_std = []
            A_std = []
            
            # Loop through events
            ex_ind = []
            for i in range(I):
                
                ############################ Set up inversion #############################
                
                event = event_list[i].split('/')[-2]
                yyyy,mth,dd,hh,mm,sec = event.split('_')[1:]
                M = event_df['Mag'].iloc[np.where(event_df['OrgT']==f'{yyyy}-{mth}-{dd} {hh}:{mm}:{sec}')[0][0]]
                
                # Do not include any M6 events if those show up
                if M < 6:
                
                    # Read in spectra data for event i at current station
                    data = np.genfromtxt(event_list[i], comments = '#', dtype = float)
                    freq = data.T[0]
                    amp = data.T[1]
                    
                    # Get theoretical corner frequency and frequency ranges (fe, fx)
                    fe,fx,fc,_,_,_ = freqrange.get_freq_range(freq,amp,M,model=freq_model)
                    
                    # Make sure frequency range is at least 10 Hz
                    if fx-fe >=10:
                        
                        fe_list.append(fe)
                        fx_list.append(fx)
                        fc_list.append(fc)
                        
                        # Obtain frequencies and amplitudes in frequency range used to obtain Kappa
                        freq_ind = np.where((data.T[0] >= fe) & (data.T[0] <= fx))[0]
                        inv_freq = data.T[0][freq_ind]
                        amplitude = (data.T[1][freq_ind]) #data from log to lin
                        std = data.T[2][freq_ind]
                        
                        kappa, A0, var = run_kappa_inv(inv_freq, amplitude, std)
                        
                        # Append tstar(kappa) and A0 to lists
                        kappa_list.append(float(kappa))
                        A0_list.append(A0)
                        
                        # Append standard deviations 
                        kappa_std.append(np.sqrt(np.abs(var[0])))#error in kappa
                        A_std.append(np.sqrt(np.abs(A0*var[1])))
                        
                        temp_stns.append(stn)
                        temp_mag.append(M)
                        temp_rhyp.append(stn_df3['rhyp'].iloc[i])
                        temp_depth.append(stn_df3['Qdep'].iloc[i]/1000) #m to km
                        temp_events.append(event)
                        # temp_rrup.append()
                
                    else:
                        ex_ind.append(i)
                        I = I-1
            
            # Check how many records are left now
            if I >= min_events:
                
                # Remove records from df if freq range was too small (< 10 Hz)
                if len(ex_ind) > 0:
                    stn_df3 = stn_df3.reset_index(drop=True)
                    stn_df3 = stn_df3.drop(ex_ind,axis=0)
                
                stn_df3['fc'] = fc_list
                stn_df3['fe'] = fe_list
                stn_df3['fx'] = fx_list
                stn_df3['A0'] = A0_list
                stn_df3['A0 std'] = A_std
                stn_df3['Kappa(s)'] = kappa_list
                stn_df3['Kappa std dev'] = kappa_std
                
                stn_df3.to_csv(f'{home_dir}/stn_flatfiles/{stn}.csv')
                
                n_stations += 1
                
                if stn_df3['Channel'].iloc[0] == 'HH':
                    bb+=1
                elif stn_df3['Channel'].iloc[0] == 'HN':
                    sm+=1
            
                # Plot spectra
                if not path.exists(f'{home_dir}/plots'):
                    makedirs(f'{home_dir}/plots')
                figpath = f'{home_dir}/plots/{stn}_spectra.png'
                
                # Plot spectra with frequency ranges indicated
                plot_spectra(stn,home_dir,min_mag,min_events,event_df,spectra_dir,figpath)
            
                # Plot kappa
                kappa,kappa_std,r_squared = plot_kappa(model_name,home_dir,stn_df3,min_events)
                stn_kappa_list.append(kappa)
                stn_kappa_std_list.append(kappa_std)
                r2_list.append(r_squared)
                stns.append(stn)
                
                n_records += I
                
                stn_list = np.append(stn_list, temp_stns)
                mag_list = np.append(mag_list, temp_mag)
                rhyp_list = np.append(rhyp_list, temp_rhyp)
                depth_list = np.append(depth_list, temp_depth)
                final_events = np.append(final_events, temp_events)

    else:
        I = len(event_list)
        for i in range(I):
            event = event_list[i].split('/')[-2]
            yyyy,mth,dd,hh,mm,sec = event.split('_')[1:]
            M = event_df['Mag'].iloc[np.where(event_df['OrgT']==f'{yyyy}-{mth}-{dd} {hh}:{mm}:{sec}')[0][0]]
            stn_list = np.append(stn_list, stn)
            mag_list = np.append(mag_list, M)
            depth_list = np.append(depth_list, stn_df3['Qdep'].iloc[i]/1000)
            rhyp_list = np.append(rhyp_list, stn_df3['rhyp'].iloc[i])
            final_events = np.append(final_events, event)

data = {'Station':stn_list,'Event':final_events,'Mag':mag_list,'Depth':depth_list,'Rhyp':rhyp_list}
final_df = pd.DataFrame(data)
final_df.to_csv(f'/Users/tnye/kappa/data/flatfiles/final_dataset_{model_name}.csv')
           
f_criteria.write(f'Average Polyfit R2: {np.mean(r2_list)}\n')
f_criteria.write('\n')
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
        
