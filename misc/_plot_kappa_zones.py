#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 15:32:50 2022

@author: tnye
"""

# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('/Users/tnye/kappa/data/kappa_zones.csv')

delta_k = np.array(df['Kappa'])[np.where(np.array(df['Kappa_Zone'])=='River Delta')[0]]
san_k = np.array(df['Kappa'])[np.where(np.array(df['Kappa_Zone'])=='San Andreas')[0]]
hay_k = np.array(df['Kappa'])[np.where(np.array(df['Kappa_Zone'])=='Hayward')[0]]

cat = [delta_k, san_k, hay_k]
labels = ['Delta', 'San Andreas', 'Hayward']

plt.figure(figsize=(12,8))
for i, a in enumerate(cat):
# kappa = kappa_list[np.where(np.array(fault['Fault_Zone']==cat))[0]]
# mean_kappa = np.mean(kappa)
# kappa_units.append(mean_kappa)
    n, bins, patches = plt.hist(a, bins=23, range=(0.005,0.115), alpha=0.5, label=labels[i])
    c=patches[0].get_facecolor()
    plt.vlines(np.mean(a), 0, 50, color=c, ls='--')

plt.ylim(0,30)
plt.xlabel('Kappa (s)')
plt.ylabel('Number of stations')
plt.legend()

plt.title('Kappa Zone')


sns.set(style = 'whitegrid')
sns.violinplot(x ="Kappa_Zone", y="Kappa", data=df)


q3, q1 = np.percentile(delta_k, [75 ,25])
print(f'Delta: {round(q1,3)}-{round(q3,3)}')
q3, q1 = np.percentile(hay_k, [75 ,25])
print(f'Haywawrd: {round(q1,3)}-{round(q3,3)}')
q3, q1 = np.percentile(san_k, [75 ,25])
print(f'San Andreas: {round(q1,3)}-{round(q3,3)}')



# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('/Users/tnye/kappa/data/kappa_zones2.csv')

west_k = np.array(df['Kappa'])[np.where(np.array(df['Kappa_Zone'])=='West')[0]]
east_k = np.array(df['Kappa'])[np.where(np.array(df['Kappa_Zone'])=='East')[0]]

cat = [west_k, east_k]
labels = ['West Bay', 'East Bay']

plt.figure(figsize=(12,8))
for i, a in enumerate(cat):
# kappa = kappa_list[np.where(np.array(fault['Fault_Zone']==cat))[0]]
# mean_kappa = np.mean(kappa)
# kappa_units.append(mean_kappa)
    n, bins, patches = plt.hist(a, bins=23, range=(0.005,0.115), alpha=0.5, label=labels[i])
    c=patches[0].get_facecolor()
    plt.vlines(np.mean(a), 0, 50, color=c, ls='--')

plt.ylim(0,30)
plt.xlabel('Kappa (s)')
plt.ylabel('Number of stations')
plt.legend()
plt.title('Kappa Zone')
plt.show()

plt.figure()
sns.set(style = 'whitegrid')
sns.violinplot(x ="Kappa_Zone", y="Kappa", data=df)
plt.show()


q3, q1 = np.percentile(west_k, [75 ,25])
print(f'West: {round(q1,3)}-{round(q3,3)}')
q3, q1 = np.percentile(east_k, [75 ,25])
print(f'East: {round(q1,3)}-{round(q3,3)}')
