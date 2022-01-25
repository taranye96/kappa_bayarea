#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 10:59:07 2021

@author: tnye
"""

###############################################################################
# Script that makes text file with station paths to be used in GMT. 
###############################################################################

# Imports 
import pandas as pd
import numpy as np

df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/full_bay_flatfile.csv')
outfile='/Users/tnye/kappa/data/paths/paths_250.txt'

idx = np.where((np.array(df['Mag']>=3.5)) & (np.array(df['rhyp']<=250)))
mag = df['Mag'][idx[0]]
rhyp = df['rhyp'][idx[0]]
newdf = df[idx[0]]
for i in {2..536100}
do
  rhyp=$(awk -v i="$i" -F, 'NR==i {print $13}' $seg_file)
  mag=$(awk -v i="$i" -F, 'NR==i {print $8}' $seg_file)
  if [[ $(echo "$mag >= 3.5" | bc) == 1 ]]; then
    if [[ $(echo "$rhyp <= 250" | bc) == 1 ]]; then
    # Get station lon, lat
    echo ">-Z1.515476448" >>/Users/tnye/kappa/data/paths/paths_250.txt
    awk -v i="$i" -F, 'NR==i {print $5,$4}' $seg_file>>$outfile
    awk -v i="$i" -F, 'NR==i {print $10,$9}' $seg_file>>$outfile
    fi
  fi

done