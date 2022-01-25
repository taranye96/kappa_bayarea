#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 11 02:39:18 2021

@author: tnye
"""

###############################################################################
# That makes histogram of kappa for the different geologyic units or
# relation to fault zone. 
###############################################################################

# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

geo = pd.read_csv('/Users/tnye/kappa/data/stn_geology2.csv')
kappa_list = pd.read_csv('/Users/tnye/kappa/traditional_method/models/test/updated_stns/updated_stns_kappa.out',delimiter='\t')[' kappa(s) ']

# categories = ['Quaternary', 'Sedimentary', 'Franciscan', 'Great Valley', 'Metamorphic', 'Volcanic' ]
# categories = ['Alluvium', 'Artificial Fill', 'Mud', 'Sand', 'Misc. Quaternary', 'Sedimentary', 'Volcanic', 'Metamorphic' ]
categories = ['Alluvium', 'Sedimentary', 'Volcanic', 'Sand', 'Melange', 'Misc. Quaternary', 'Artificial Fill', 'Metamorphic' ]
kappa_units = []
plt.figure(figsize=(12,8))
a=0
for i, cat in enumerate(categories):
    kappa = kappa_list[np.where(np.array(geo['Geology']==cat))[0]]
    mean_kappa = np.mean(kappa)
    print(f'{cat}: {round(mean_kappa,3)}')
    # med_kappa = np.median(kappa)
    kappa_units.append(mean_kappa)
    # kappa_units.append(med_kappa)
    
   
    n, bins, patches = plt.hist(kappa, alpha=0.5, label=cat)
    c=patches[0].get_facecolor()
    plt.vlines(mean_kappa, 0, 50, color=c, ls='--')
    # plt.annotate(labels[i], (mean_kappa, 28-a))
    plt.ylim(0,30)
    plt.xlabel('Kappa (s)')
    plt.ylabel('Number of stations')
    plt.legend()
    a+=2

plt.show()
plt.savefig('/Users/tnye/kappa/plots/kappa_vs_geology.png', dpi=300)
plt.close()


###############################################################################

# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

fault = pd.read_csv('/Users/tnye/kappa/data/fault_zone_stns.csv')
kappa_list = pd.read_csv('/Users/tnye/kappa/traditional_method/models/test/updated_stns/updated_stns_kappa.out',delimiter='\t')[' kappa(s) ']

# categories = ['Quaternary', 'Sedimentary', 'Franciscan', 'Great Valley', 'Metamorphic', 'Volcanic' ]
# categories = ['Alluvium', 'Artificial Fill', 'Mud', 'Sand', 'Misc. Quaternary', 'Sedimentary', 'Volcanic', 'Metamorphic' ]
categories = ['outside', 'within']
kappa_units = []
plt.figure(figsize=(12,8))
a=0
for i, cat in enumerate(categories):
    kappa = kappa_list[np.where(np.array(fault['Fault_Zone']==cat))[0]]
    mean_kappa = np.mean(kappa)
    kappa_units.append(mean_kappa)
    # med_kappa = np.median(kappa)
    # kappa_units.append(med_kappa)
    
   
    n, bins, patches = plt.hist(kappa, bins=23, range=(0.005,0.115), alpha=0.5, label=cat)
    c=patches[0].get_facecolor()
    plt.vlines(mean_kappa, 0, 50, color=c, ls='--')
    # plt.vlines(med_kappa, 0, 50, color=c, ls='--')
    # plt.annotate(labels[i], (mean_kappa, 28-a))
    plt.ylim(0,30)
    plt.xlabel('Kappa (s)')
    plt.ylabel('Number of stations')
    plt.legend()
    a+=2

plt.title('Fault Damage Zone')
plt.show()
plt.savefig('/Users/tnye/kappa/plots/kappa_vs_faults.png', dpi=300)
plt.close()



# geo = pd.read_csv('/Users/tnye/kappa/data/stn_geology.csv')
# for i in range(len(geo)):
#     unit = geo['Unit'].iloc[i]
#     description = geo.iloc[i].Description
   
#     if unit != float('nan'):
#         if [(unit[0] == 'a') | (unit[0] == 'Q')]:
#             if [(geo['Latitude'].iloc[i] >= 37.45) & (geo['Longitude'] >= -122.38)]:
#                 geo['Geology'].iloc[i] = 'East Quaternary'
#             elif [(geo['Latitude'].iloc[i] <= 37.45) and (geo['Longitude'] <= -122.38)]:
#                 geo['Geology'].iloc[i] = 'West Quaternary'
#             else:
#                 geo['Geology'] = ['Other Quaternary']
            
#         elif 'sediment' in description.lower():
#             if geo['Latitude'].iloc[i] >= 37.45 and geo['Longitude'] >= -122.38:
#                 geo['Geology'].iloc[i] = 'East Sedimentary'
#             elif geo['Latitude'].iloc[i] <= 37.45 and geo['Longitude'] <= -122.38:
#                 geo['Geology'].iloc[i] = 'West Sedimentary'
#             else:
#                 geo['Geology'] = ['Other Quaternary']
           
#             if 'complex' in description.lower():
#                 geo['Geology'].iloc[i] = 'Sedimentary Complex'
#         elif 'meta' in description.lower() or 'ite' in description.lower():
#             geo['Geology'].iloc[i] = 'Metamorphic'
#             if 'complex' in description.lower():
#                 geo['Geology'].iloc[i] = 'Metamorphic Complex'
#         elif 'volcanic' in description.lower() or 'granite' in description.lower():
#             geo['Geology'].iloc[i] = 'Volcanic'
#             if 'complex' in description.lower():
#                 geo['Geology'].iloc[i] = 'Volcanic Complex'
#         else:
#             geo['Geology'].iloc[i] = 'unknown'
        
# geo.to_csv('/Users/tnye/kappa/data/stn_geology.csv')
    
    