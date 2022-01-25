#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 15:36:55 2021

@author: tnye
"""

###############################################################################
# Script used to combine broadband and strong motion station data flatfiles
# into one main flatfile. 
###############################################################################

import pandas as pd
import glob

path = '/Users/tnye/kappa/data/flatfiles/sm_acc/collected_mpi_flatfiles' # use your path
all_files = sorted(glob.glob(path + "/*.csv"))

li = []

df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/SNR_5_file.csv', index_col=None, header=0)
li.append(df)
df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/cesmd_flatfile.csv', index_col=None, header=0)
li.append(df)

frame = pd.concat(li, axis=0, ignore_index=True)

frame.to_csv('/Users/tnye/kappa/data/flatfiles/full_bay_flatfile.csv')

