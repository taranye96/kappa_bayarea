#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 16:09:00 2022

@author: tnye
"""

###############################################################################
# Script used to make the kappa vs site residual figures for the manuscript. 
###############################################################################

# Imports 
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import seaborn as sns
from scipy import stats
from statsmodels.stats.power import tt_ind_solve_power
from scipy.stats import pearsonr
from plotting_fns import place_column_text
from matplotlib.ticker import MultipleLocator, ScalarFormatter


class Format:
    end = '\033[0m'
    underline = '\033[4m'

model_name = 'model1'
home_dir = f'/Users/tnye/kappa/traditional_method/models/{model_name}'

site_res_file = pd.read_csv(f'/Users/tnye/kappa/residual_decomp/output/Bay_dS_meml_{model_name}.txt',sep='\t')
stn_file = pd.read_csv(f'/Users/tnye/kappa/traditional_method/models/{model_name}/{model_name}_kappa.out',sep='\t')
site_ind = []
pga_ds = []
pgv_ds = []
pga_ds_ms = []
pga_ds_est = []
pgv_ds_ms = []
pgv_ds_est = []
vs30_ms = []
vs30_est = []
vs30_idx = []
k_ms = []
k_est = []
kappa = []
stns = np.unique(stn_file['#site '])

vs30_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/USGS_vs30.csv',encoding='ISO-8859-1')
vs30_df["VS30 Measurement 1 (m/s)"] = pd.to_numeric(vs30_df["VS30 Measurement 1 (m/s)"], downcast="float", errors="ignore")
vs30_list = np.array(vs30_df['VS30 Measurement 1 (m/s)'])
stn_list = np.array(vs30_df['Station Number'])
nga_stns_df = pd.read_csv('/Users/tnye/kappa/data/flatfiles/NGA_West2_SiteDatabase_V032.csv')

vs30_file = np.genfromtxt(f'/Users/tnye/kappa/data/Vs30/{model_name}_vs30.txt') 

kappa_file = pd.read_csv(f'{home_dir}/{model_name}_kappa.out',sep='\t')
kappa_list = kappa_file[' kappa(s) ']
kappa_std = kappa_file[' kappa_std ']

for i, stn in enumerate(stns):
    ind = np.where(site_res_file['Station Name']==stn)[0][0]
    pga_ds.append(site_res_file['lnASK14_PGA_Res dS'].iloc[ind])
    pgv_ds.append(site_res_file['lnASK14_PGV_Res dS'].iloc[ind])
    kappa.append(kappa_list[i])
    
    if stn in stn_list:
        stn_idx = np.where(stn_list==stn)[0][0]
        v = vs30_list[stn_idx]
        if v != 0.0:
            # vs30.append(v)
            vs30_ms.append(v)
            vs30_idx.append('#F1830E')
            k_ms.append(kappa_list[i])
            pga_ds_ms.append(site_res_file['lnASK14_PGA_Res dS'].iloc[ind])
            pgv_ds_ms.append(site_res_file['lnASK14_PGV_Res dS'].iloc[ind])
        # Check if in NGA flatfile
        elif stn in np.array(nga_stns_df['Station ID']):
            ind = np.where(nga_stns_df['Station ID']==stn)[0][0]
            vs30_est.append(nga_stns_df['Vs30 for Analysis (m/s)'].iloc[ind])
            k_est.append(kappa_list[i])
            pga_ds_est.append(site_res_file['lnASK14_PGA_Res dS'].iloc[ind])
            pgv_ds_est.append(site_res_file['lnASK14_PGV_Res dS'].iloc[ind])
        else:
            # vs30.append(vs30_file[i])
            vs30_est.append(vs30_file[i])
            vs30_idx.append('#0EA2F1')
            k_est.append(kappa_list[i])
            pga_ds_est.append(site_res_file['lnASK14_PGA_Res dS'].iloc[ind])
            pgv_ds_est.append(site_res_file['lnASK14_PGV_Res dS'].iloc[ind])
    else:
        print(f'no data for stn {stn}')


# dataset_dict = {'Kappa':kappa, 'Vs30':vs30, 'PGA_site_res':pga_ds, 'PGV_site_res':pgv_ds}
# df = pd.DataFrame(data=dataset_dict)
# df['color']=vs30_idx

ms_dict = {'Kappa':k_ms, 'Vs30':vs30_ms, 'PGA_site_res':pga_ds_ms, 'PGV_site_res':pgv_ds_ms}
est_dict = {'Kappa':k_est, 'Vs30':vs30_est, 'PGA_site_res':pga_ds_est, 'PGV_site_res':pgv_ds_est}
ms_df = pd.DataFrame(data=ms_dict)
est_df = pd.DataFrame(data=est_dict)
kappa_dict = {'Kappa':kappa, 'PGA_site_res':pga_ds, 'PGV_site_res':pgv_ds}
kappa_df = pd.DataFrame(kappa_dict)

#%%
################################ Kappa vs Vs30 ################################

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'

fig, ax = plt.subplots(figsize=(6,5))
ax.set_xlim(np.log10(np.array([125,1000])))
ax.set_ylim(np.log10(np.array([0.003,0.12])))
sns.regplot(ax=ax, x=np.log10(est_df["Vs30"]), y=np.log10(est_df["Kappa"]), robust=False, truncate=False, label='$V_{S30,est}$', order=2)
sns.regplot(ax=ax, x=np.log10(ms_df["Vs30"]), y=np.log10(ms_df["Kappa"]), robust=False, truncate=False, label='$V_{S30,meas}$', order=2)
ax.set_ylabel('$\kappa_0$ (s)',fontsize=14)
ax.set_xlabel('$V_{S30}$ (m/s)',fontsize=14)
ax.tick_params(which='minor', top=True, right=True)
ax.grid(linestyle='--',alpha=0.5,which='both')
ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.25), ncol=2)
ax.collections[1].set_label('95% Confidence interval')
plt.subplots_adjust(top=0.96, bottom=0.2, left = 0.125, right = 0.96)

xticks = np.log10(np.array([1000]))
yticks = np.log10(np.array([0.01,0.1]))
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_yticks(np.log10(np.array([0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.11,0.12])), minor=True)
ax.set_xticks(np.log10(np.array([200,300,400,500,600,700,800,900])), minor=True)
formatter = lambda x, pos: f'{10 ** x:g}'
ax.get_xaxis().set_major_formatter(formatter)
ax.get_yaxis().set_major_formatter(formatter)

    
plt.show()
plt.savefig(f'/Users/tnye/kappa/plots/paper/kappa_vs30_2nd_loglog.png', dpi=300)


#%%
####################### Kappa vs Vs30: Linear in loglog #######################

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'

fig, ax = plt.subplots(figsize=(6,5))
ax.set_xlim(np.log10(np.array([125,1000])))
ax.set_ylim(np.log10(np.array([0.003,0.12])))
sns.regplot(ax=ax, x=np.log10(est_df["Vs30"]), y=np.log10(est_df["Kappa"]), robust=False, truncate=False, label='$V_{S30,est}$')
sns.regplot(ax=ax, x=np.log10(ms_df["Vs30"]), y=np.log10(ms_df["Kappa"]), robust=False, truncate=False, label='$V_{S30,meas}$')
ax.set_ylabel('$\kappa_0$ (s)',fontsize=14)
ax.set_xlabel('$V_{S30}$ (m/s)',fontsize=14)
ax.tick_params(direction="out",labelright=False,top=True,right=True,length=5,labelsize=12)
ax.tick_params(which='minor', top=True, right=True)
ax.grid(linestyle='--',alpha=0.5,which='both')
ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.25), ncol=2)
ax.collections[1].set_label('95% Confidence interval')
plt.subplots_adjust(top=0.96, bottom=0.2, left = 0.125, right = 0.96)

xticks = np.log10(np.array([1000]))
yticks = np.log10(np.array([0.01,0.1]))
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_yticks(np.log10(np.array([0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.11,0.12])), minor=True)
ax.set_xticks(np.log10(np.array([200,300,400,500,600,700,800,900])), minor=True)
formatter = lambda x, pos: f'{10 ** x:g}'
ax.get_xaxis().set_major_formatter(formatter)
ax.get_yaxis().set_major_formatter(formatter)

pop_size = len(vs30_est)
r,pval = pearsonr(np.log10(k_est),np.log10(vs30_est))
pval_text = f'{round(pval,2)}'
if round(pval,2) == 0:
    pval_text = '< 0.001'
power = tt_ind_solve_power(effect_size=r,nobs1=pop_size,alpha=0.05)
textstr = '\n'.join((
    '$V_{S30,est}$',
    f'Pearson R: {round(r,2)}',
    f'P-value: {pval_text}',
    f'Power: {round(power,2)}'))
props = dict(boxstyle='round', facecolor='white', alpha=0.7)
ax.text(0.05, 0.05, textstr, transform=ax.transAxes, fontsize=12, va='bottom', ha='left', bbox=props)

pop_size = len(vs30_ms)
r,pval = pearsonr(np.log10(k_ms),np.log10(vs30_ms))
pval_text = f'{round(pval,2)}'
if round(pval,2) == 0:
    pval_text = '< 0.001'
power = tt_ind_solve_power(effect_size=r,nobs1=pop_size,alpha=0.05)
textstr = '\n'.join((
    '$V_{S30,meas}$',
    f'Pearson R: {round(r,2)}',
    f'P-value: {pval_text}',
    f'Power: {round(power,2)}'))
props = dict(boxstyle='round', facecolor='white', alpha=0.7)
# ax.text(0.76, 0.975, textstr, transform=ax.transAxes, fontsize=10, va='top', ha='left', bbox=props)
ax.text(0.345, 0.05, textstr, transform=ax.transAxes, fontsize=12, va='bottom', ha='left', bbox=props)
    # 
# plt.show()
plt.savefig(f'/Users/tnye/kappa/plots/paper/kappa_vs30_linear_loglog.png', dpi=300)



