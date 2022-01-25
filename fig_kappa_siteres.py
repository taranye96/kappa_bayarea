#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 14:42:42 2021

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
        else:
            # vs30.append(vs30_file[i])
            vs30_est.append(vs30_file[i])
            vs30_idx.append('#0EA2F1')
            k_est.append(kappa_list[i])
            pga_ds_est.append(site_res_file['lnASK14_PGA_Res dS'].iloc[ind])
            pgv_ds_est.append(site_res_file['lnASK14_PGV_Res dS'].iloc[ind])


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

fig, ax = plt.subplots(figsize=(6.5,5))
# sns.regplot(ax=ax, x="Kappa", y="Vs30", data=df, robust=False, scatter=False)
# ax.scatter(kappa, vs30, c=vs30_idx)
ax.set_xlim(0,0.12)
ax.set_ylim(100,850)
sns.regplot(ax=ax, x="Kappa", y="Vs30", data=est_df, robust=False, truncate=False, label='$V_{S30,est}$')
sns.regplot(ax=ax, x="Kappa", y="Vs30", data=ms_df, robust=False, truncate=False, label='$V_{S30,meas}$')
ax.set_ylabel('$V_{S30}$ (m/s)',fontsize=12)
ax.set_xlabel('$\kappa_0$ (s)',fontsize=12)
ax.tick_params(direction="out",labelright=False,top=True,right=True,labelsize=12)
ax.yaxis.set_minor_locator(MultipleLocator(100))
ax.yaxis.set_major_locator(MultipleLocator(200))
ax.xaxis.set_minor_locator(MultipleLocator(0.01))
ax.xaxis.set_major_locator(MultipleLocator(0.02))
ax.tick_params(which='minor', top=True, right=True)
ax.grid(linestyle='--',alpha=0.5)
# ax.legend(loc='upper left')
ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.25), ncol=2)
ax.collections[1].set_label('95% Confidence interval')

pop_size = len(vs30_ms)
r,pval = pearsonr(k_ms,vs30_ms)
pval_text = f'{round(pval,2)}'
if round(pval,2) == 0:
    pval_text = '< 0.001'
power = tt_ind_solve_power(effect_size=r,nobs1=pop_size,alpha=0.05)
textstr = '\n'.join((
    '$V_{S30,meas}$',
    f'Pearson R: {round(r,3)}',
    f'P-value: {pval_text}',
    f'Power: {round(power,2)}'))
props = dict(boxstyle='round', facecolor='white', alpha=0.7)
ax.text(0.76, 0.975, textstr, transform=ax.transAxes, fontsize=10, va='top', ha='left', bbox=props)

pop_size = len(vs30_est)
r,pval = pearsonr(k_est,vs30_est)
pval_text = f'{round(pval,2)}'
if round(pval,2) == 0:
    pval_text = '< 0.001'
power = tt_ind_solve_power(effect_size=r,nobs1=pop_size,alpha=0.05)
textstr = '\n'.join((
    '$V_{S30,est}$',
    f'Pearson R: {round(r,3)}',
    f'P-value: {pval_text}',
    f'Power: {round(power,2)}'))
props = dict(boxstyle='round', facecolor='white', alpha=0.7)
ax.text(0.76, 0.78, textstr, transform=ax.transAxes, fontsize=10, va='top', ha='left', bbox=props)
plt.subplots_adjust(top=0.96, bottom=0.2, left = 0.105, right = 0.96)
    
plt.show()
plt.savefig(f'/Users/tnye/kappa/plots/paper/kappa_vs30_{model_name}.png', dpi=300)


################################ Set up figure ################################
#%%
# Stats parameters
kappa_metrics = np.array([[kappa,pga_ds],[kappa,pgv_ds]])
ms_metrics = np.array([[vs30_ms,pga_ds_ms],[vs30_ms,pgv_ds_ms]])
est_metrics = np.array([[vs30_est,pga_ds_est],[vs30_est,pgv_ds_est]])
sig_level = 0.05

# Figure parameters
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'

fig = plt.figure(figsize=(10,14))
ax1 = fig.add_subplot(2,2,1)
ax2 = fig.add_subplot(2,2,2, sharey = ax1)
ax3 = fig.add_subplot(2,2,3, sharex = ax1)
ax4 = fig.add_subplot(2,2,4, sharex=ax2, sharey=ax3)

# Kappa vs PGA
ax1.grid(linestyle='--',alpha=0.5)
# sns.regplot(ax=ax1, x="Kappa", y="PGA_site_res", data=df, robust=True, line_kws={'label': 'Regression line'})
ax1.set_ylim(-2.25,1.25)
ax1.set_xlim(0,0.12)
sns.regplot(ax=ax1, x="Kappa", y="PGA_site_res", data=kappa_df, robust=False, truncate=False,)
# ax1.set_xlabel('')
ax1.set_xlabel('$\kappa_0$ (s)',fontsize=12)
ax1.set_ylabel('ln(PGA $\delta S_j$)',fontsize=12)
ax1.tick_params(direction="out",labelbottom=True,labelright=False,top=True,right=True,labelsize=12)
ax1.yaxis.set_minor_locator(MultipleLocator(0.5))
ax1.yaxis.set_major_locator(MultipleLocator(1))
ax1.xaxis.set_minor_locator(MultipleLocator(0.01))
ax1.xaxis.set_major_locator(MultipleLocator(0.02))
ax1.tick_params(which='minor', top=True, right=True)
ax1.collections[1].set_label('95% Confidence interval')
ax1.text(-0.175,1.0,'(a)',transform=ax1.transAxes,fontsize=12,va='top',ha='right',weight='bold')

# Vs30 vs PGA
ax2.grid(linestyle='--',alpha=0.5)
ax2.set_ylim(-2.25,1.25)
ax2.set_xlim(100,850)
sns.regplot(ax=ax2, x="Vs30", y="PGA_site_res", data=est_df, robust=False, truncate=False, label='$V_{S30,est}$')
sns.regplot(ax=ax2, x="Vs30", y="PGA_site_res", data=ms_df, robust=False, truncate=False, label='$V_{S30,meas}$')
ax2.set_xlabel('$V_{S30}$ (m/s)',fontsize=12)
ax2.set_ylabel('ln(PGA $\delta S_j$)',fontsize=12)
ax2.collections[1].set_label('95% Confidence interval')
ax2.collections[3].set_label('95% Confidence interval')
ax2.tick_params(direction="out",labelleft=True,labelbottom=True,labelright=False,top=True,right=True,labelsize=12)
ax2.yaxis.set_minor_locator(MultipleLocator(0.5))
ax2.yaxis.set_major_locator(MultipleLocator(1))
ax2.xaxis.set_minor_locator(MultipleLocator(100))
ax2.xaxis.set_major_locator(MultipleLocator(200))
ax2.tick_params(which='minor', top=True, right=True)
ax2.text(-0.175,1.0,'(b)',transform=ax2.transAxes,fontsize=12,va='top',ha='right',weight='bold')

# Kappa vs PGV
ax3.grid(linestyle='--',alpha=0.5)
ax3.set_ylim(-2.25,1.25)
ax3.set_xlim(0.0,0.12)
sns.regplot(ax=ax3, x="Kappa", y="PGV_site_res", data=kappa_df, truncate=False, robust=False)
ax3.set_xlabel('$\kappa_0$ (s)',fontsize=12)
ax3.set_ylabel('ln(PGV $\delta S_j$)',fontsize=12)
ax3.tick_params(direction="out",top=True,right=True)
ax3.yaxis.set_minor_locator(MultipleLocator(0.5))
ax3.yaxis.set_major_locator(MultipleLocator(1))
ax3.xaxis.set_minor_locator(MultipleLocator(0.01))
ax3.xaxis.set_major_locator(MultipleLocator(0.02))
ax3.tick_params(which='minor', top=True, right=True)
ax3.collections[1].set_label('95% Confidence interval')
ax3.text(-0.175,1.0,'(c)',transform=ax3.transAxes,fontsize=12,va='top',ha='right',weight='bold')
# ax3.legend(loc='upper right')
ax3.legend(loc="lower center", bbox_to_anchor=(0.5, -0.35))


# Vs30 vs PGV
ax4.grid(linestyle='--',alpha=0.5)
ax4.set_ylim(-2.25,1.25)
ax4.set_xlim(100,850)
sns.regplot(ax=ax4, x="Vs30", y="PGV_site_res", data=est_df, robust=False, truncate=False, label='$V_{S30,est}$')
sns.regplot(ax=ax4, x="Vs30", y="PGV_site_res", data=ms_df, robust=False, truncate=False, label='$V_{S30,meas}$')
ax4.set_xlabel('$V_{S30}$ (m/s)',fontsize=12)
# ax4.set_ylabel('')
ax4.set_ylabel('ln(PGV $\delta S_j$)',fontsize=12)
ax4.collections[1].set_label('95% Confidence interval')
ax4.collections[3].set_label('95% Confidence interval')
ax4.tick_params(direction="out",labelleft=True,top=True,right=True,labelsize=12)
ax4.yaxis.set_minor_locator(MultipleLocator(0.5))
ax4.yaxis.set_major_locator(MultipleLocator(1))
ax4.xaxis.set_minor_locator(MultipleLocator(100))
ax4.xaxis.set_major_locator(MultipleLocator(200))
ax4.tick_params(which='minor', top=True, right=True)
ax4.text(-0.175,1.0,'(d)',transform=ax4.transAxes,fontsize=12,va='top',ha='right',weight='bold')
# ax4.legend(loc='upper right')
ax4.legend(loc="lower center", bbox_to_anchor=(0.5, -0.45), ncol=2)

# Add stats to subplots
# Kappa metrics
for i, ax in enumerate([ax1,ax3]):
    pop_size = len(kappa_metrics[0][0])
    r,pval = pearsonr(kappa_metrics[i][0],kappa_metrics[i][1])
    pval_text = f'{round(pval,2)}'
    if round(pval,2) == 0:
        pval_text = '< 0.001'
    power = tt_ind_solve_power(effect_size=r,nobs1=pop_size,alpha=sig_level)
    textstr = '\n'.join((
        f'Pearson R: {round(r,2)}',
        f'P-value: {pval_text}',
        f'Statistical Power: {round(power,2)}'))
    props = dict(boxstyle='round', facecolor='white', alpha=0.7)
    ax.text(0.05, 0.05, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='bottom', bbox=props)

for i, ax in enumerate([ax2,ax4]):
    
    est_pop_size = len(est_metrics[0][0])
    est_r,est_pval = pearsonr(est_metrics[i][0],est_metrics[i][1])
    est_pval_text = f'{round(est_pval,2)}'
    if round(est_pval,2) == 0:
        est_pval_text = '< 0.001'
    est_power = tt_ind_solve_power(effect_size=est_r,nobs1=est_pop_size,alpha=sig_level)
    est_txtstr = '\n'.join((
        '$V_{S30,est}$',
        f'Pearson R: {round(est_r,2)}',
        f'P-value: {est_pval_text}',
        f'Statistical Power: {round(est_power,2)}'))
    
    ms_pop_size = len(ms_metrics[0][0])
    ms_r,ms_pval = pearsonr(ms_metrics[i][0],ms_metrics[i][1])
    ms_pval_text = f'{round(ms_pval,2)}'
    if round(ms_pval,3) == 0:
        ms_pval_text = '< 0.001'
    ms_power = tt_ind_solve_power(effect_size=ms_r,nobs1=ms_pop_size,alpha=sig_level)
    ms_txtstr = '\n'.join((
        '$V_{S30,meas}}$',
        f'Pearson R: {round(ms_r,2)}',
        f'P-value: {ms_pval_text}',
        f'Statistical Power: {round(ms_power,2)}'))
    
    r'$\underline{Measured V_{S30}}$'
    
    props = dict(boxstyle='round', facecolor='white', alpha=0.7)
    ax.text(0.05, 0.05, est_txtstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='bottom', bbox=props)
    # ax.text(0.05, 0.95, ms_txtstr, transform=ax.transAxes, fontsize=10,
    #     va='top', bbox=props)
    ax.text(0.6, 0.05, ms_txtstr, transform=ax.transAxes, fontsize=10,
        va='bottom', bbox=props)
    
plt.subplots_adjust(wspace=0.35, hspace=0.3, left=0.1, right=0.97, top=0.98, bottom=0.17)
plt.show()
plt.savefig(f'/Users/tnye/kappa/plots/paper/site_param_correlation_{model_name}.png', dpi=300)
