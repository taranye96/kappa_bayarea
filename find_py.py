#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 13:27:27 2022

@author: tnye
"""

from glob import glob

# Glob all python files. recursive=True allows us to loop through all subdirectories, not just the first 
py_files = glob('/Users/tnye/**/*.py', recursive=True)

# Loop through python files to find the file that makes the datamap figure
for file in py_files:
    with open(file) as f:
        if '/supplemental/' in f.read():
            print(file)
            # break
        # else:
        #     continue
