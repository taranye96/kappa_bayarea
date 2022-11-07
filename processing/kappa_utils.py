#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 15:46:51 2019

@author: eking

This code was downloded from https://github.com/vSahakian/kappabay.
"""

import numpy as np
import pandas as pd
from obspy.core.utcdatetime import UTCDateTime
import obspy
import datetime
from datetime import datetime, timedelta
import glob
import os
import os.path as path
import time

#%%
def cut_swave(infile, cutfile, Sarriv, cutstart, cutend):
    #reads in a sac file
    #full path to output file
    #cuts at t1 sec before s arrival and t2 seconds after
    from obspy import read    
    import datetime
    from obspy.core.utcdatetime import UTCDateTime
    import matplotlib.pyplot as plt
    
   
    # Read in obspy stream 
    stream = read(infile)
    tr = stream[0]
    
    ## Get start time:
    streamstart = tr.stats['starttime']
    
    
    start = UTCDateTime(Sarriv - datetime.timedelta(0,cutstart))
    end = UTCDateTime(Sarriv + datetime.timedelta(0,cutend))
    # start = UTCDateTime(orig + datetime.timedelta(0,stime-cutstart))
    # end = UTCDateTime(orig + datetime.timedelta(0,stime+cutend))
  
#    print(start, end)
#    print(cuttime - datetime.timedelta(seconds = t1),cuttime + datetime.timedelta(seconds = t2))
#    print(tr.slice(cuttime - datetime.timedelta(seconds = t1),cuttime + datetime.timedelta(seconds = t2), nearest_sample=True))
#    cut_trace = tr.slice(cuttime - datetime.timedelta(seconds = t1),cuttime + datetime.timedelta(seconds = t2), nearest_sample=True)#make a new trace by slicing
    cut_trace = tr.slice(start, end, nearest_sample=True)#make a new trace by slicing
    cut_trace.write(cutfile, format = 'mseed')    

#%%
def check_diff(con_dir1, con_dir2):
    
    data_dir = '/home/eking/Documents/internship/data/events/SNR_5'
    
    fig_dir = '/home/eking/Documents/internship/data/figs/SNR_5'
    
    con_stn1 = glob.glob(data_dir + '/'+ con_dir1 + '/[!2]*.out')
    
    con_stn2 = glob.glob(data_dir + '/'+ con_dir2 + '/[!2]*.out')
    
    for i in range(len(con_stn1)):#for each station
    #make each row into an array
        stn1 = np.genfromtxt(con_stn1[i])
        stn2 = np.genfromtxt(con_stn2[i])
        stnid = path.basename(con_stn1[i]).split('.')[0]
        freq = stn1.T[0]
        amp1 = stn1.T[1]
        std1 = stn1.T[2]
        amp2 = stn2.T[1]
        std2 = stn2.T[2]
        diff = np.abs(amp1 - amp2)
        stdiff = np.abs(std1 -std2)
        fig = plt.figure(figsize = (12,10))
        plt.axes().set_xscale("log")
        plt.axes().set_yscale("log")
        plt.ylabel('diffrence ', fontsize = 16)
        plt.xlim(0.5,70)
        plt.errorbar(freq , diff, yerr=stdiff, fmt = 'b')
        plt.grid()
        #plt.loglog(freq, Brune_list[ind], color = 'blue', label = 'Brune spectra')
        #plt.legend(loc = 'lower right', fontsize = 16)
        plt.xlabel('Frequency (Hz)', fontsize = 16)
        plt.title(stnid, fontsize = 16)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.tick_params(axis='both', which='both', length = 5, width = 1)
       # plt.text(0.7, .1, 'Median log(diff) 1-32.7 Hz (demeaned): ' + str(round(sum_list[ind],3)), fontsize = 16)
        plt.savefig(fig_dir + '/difflog/' + stnid + '.png')
        plt.show()
    
#%%
#
#def plot_raw_spectra(spectrafile):
#    spec = np.genfromtxt(spectrafile)
#    freq = spec.T[0]
#    amp = spec.T[1]
#    std = spec.T[2]
#    fig = plt.figure(figsize = (12,10))
#    plt.axes().set_xscale("log")
#    plt.axes().set_yscale("log")
#    plt.ylabel('Velocity amplitude (m)', fontsize = 16)
#    plt.xlabel('Frequency (Hz)', fontsize = 16)
#    #plt.xlim(0.5,50)
#    #plt.ylim(10e-6,3)
#    plt.minorticks_on
#    plt.grid(True, which='both')
#    plt.tick_params(axis='both', which='major', labelsize=15)
#    plt.tick_params(axis='both', which='both', length = 5, width = 1)
#    plt.title(path.basename(spectrafile), fontsize = 16)
#    plt_ind = np.argmin(np.abs(freq - 35))
#        
#    plt.errorbar(freq[0:plt_ind] , amp[0:plt_ind], yerr=std[0:plt_ind])#
#    plt.savefig('/home/eking/Documents/internship/data/Kappa/raw_spectra/'  + path.basename(spectrafile) + '.png')
#    plt.show()
#brune_dir = '/home/eking/Documents/internship/data/Kappa/SNR_3/75_bins/15/cut_record_spectra_15/Event_2010_06_28_14_47_04/'
#spectra_files = glob.glob(brune_dir + '*.out')
#print(len(spectra_files))
#
#for i in spectra_files:
#    plot_raw_spectra(i)

#plot_raw_spectra()
#%%
#combines parallelized flatfiles after data DL
def flatfile_comb(path0,path1,path2,path3,path4,path5,savepath):
    f1 = pd.read_csv(path1)
    f2 = pd.read_csv(path2)
    f3 = pd.read_csv(path3)
    f4 = pd.read_csv(path4)
    f5 = pd.read_csv(path5)
    f0 = pd.read_csv(path0)
    
    flatfile = pd.concat([f0,f1,f2,f3,f4,f5],axis=0,ignore_index=True)
    flatfile.to_csv(savepath)


#flatfile_comb('/home/eking/Documents/internship/data/collected_mpi_flatfiles/collected_mpi_flatfilescollected_flatfile+0.csv','/home/eking/Documents/internship/data/collected_mpi_flatfiles/collected_mpi_flatfilescollected_flatfile+1.csv','/home/eking/Documents/internship/data/collected_mpi_flatfiles/collected_mpi_flatfilescollected_flatfile+2.csv','/home/eking/Documents/internship/data/collected_mpi_flatfiles/collected_mpi_flatfilescollected_flatfile+3.csv','/home/eking/Documents/internship/data/collected_mpi_flatfiles/collected_mpi_flatfilescollected_flatfile+4.csv','/home/eking/Documents/internship/data/collected_mpi_flatfiles/collected_mpi_flatfilescollected_flatfile+5.csv','/home/eking/Documents/internship/data/collected_mpi_flatfiles/PhatPhile.csv')
#
#bigfile = pd.read_csv('/home/eking/Documents/internship/data/collected_mpi_flatfiles/PhatPhile.csv')
#print(len(bigfile))
#%%
def comp_SNR(stream, noisedur, sigstart, sigdur):
    
    time = stream[0].times()
    data = stream[0].data
    noisestart = 1
    noiseend = noisestart + noisedur
    sigend = sigstart + sigdur
    # find indecies
    noise_ind = np.where((time<=noiseend) & (time>=noisestart))[0]
    noise_amp = data[noise_ind]
    noise_avg = float(np.std(noise_amp))
    
    sig_ind = np.where((time<=sigend) & (time>=sigstart))[0]
    sig_amp = data[sig_ind]
    sig_avg = float(np.std(sig_amp))
    
    if noise_avg == 0:
        noise_avg = 1E-10
    if sig_avg == 0:
        sig_avg = 1E-10
        
    SNR=sig_avg/noise_avg
    

    
    
    return SNR
    
    
#%%
#calculates the number of stations each event was recorded at
def num_stas(flatfile):
    alldata = pd.read_csv(flatfile)
    allstas = alldata['Name']
    allevs = alldata['Quake#']
    uniq_stas = np.unique(allstas)
    uniq_evs = np.unique(allevs)
    num_sta = []
    for i in range(len(uniq_evs)):
        i_ev = uniq_evs[i]
        i_ev_ind = np.where(alldata['Quake#']==i_ev)
        i_ev_stas = alldata.loc[i_ev_ind]['Name']
        i_numstas = len(i_ev_stas)
        num_sta.append(i_numstas)
    return num_sta, uniq_stas
#%%
# calculates the number of events each station recorded
def num_evs(flatfile):
    alldata = pd.read_csv(flatfile)
    allstas = alldata['Name']
    allevs = alldata['Quake#']
    uniq_stas = np.unique(allstas)
    uniq_evs = np.unique(allevs)
    num_ev = []
    for i in range(len(uniq_stas)):
        i_st = uniq_stas[i]
        i_st_ind = np.where(alldata['Name']==i_st)
        i_sta_evs = alldata.loc[i_st_ind]['Quake#']
        i_numevs = len(i_sta_evs)
        num_ev.append(i_numevs)
    return num_ev, uniq_evs

#%%

        


