# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
from scipy.stats import linregress
from matplotlib import gridspec
from lib.io import load_CARRA_data, load_promice_old
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
station= 'QAS_U'

# path_aws = 'output/KAN_M_100_layers_3'
# path_carra = 'output/KAN_M_100_layers_4'

path_aws = 'C:/Users/bav/OneDrive - Geological survey of Denmark and Greenland/Code/GEUS model/Output/QAS_U_0_IWC_CL_100_layers_4'
path_carra = 'output/QAS_U_100_layers_2'

df_out_aws = xr.open_dataset(path_aws+'/'+station+'_surface.nc').to_dataframe()
df_out_aws.index = df_out_aws.index.round('H')
df_in_aws = load_promice_old("QAS_U_CARRA.txt")
df_out_aws [ df_in_aws.columns] = df_in_aws.values
df_out_aws.index = df_out_aws.index - pd.Timedelta('1D')
del df_in_aws

df_out_carra = xr.open_dataset(path_carra+'/'+station+'_surface.nc').to_dataframe()
df_out_carra.index = df_out_carra.index.round('H')
df_in_carra = load_CARRA_data("./input/weather data/CARRA_at_AWS.nc", station)
df_out_carra [ df_in_carra.columns] = df_in_carra
del df_in_carra

df_out_aws['melt_cumul'] = df_out_aws.melt_mweq.cumsum()
df_out_carra['melt_cumul'] = df_out_carra.melt_mweq.cumsum()
# df_out_aws['runoff_cumul'] = df_out_aws.zrogl.cumsum()
df_out_carra['runoff_cumul'] = df_out_carra.zrogl.cumsum()
df_out_aws['SR_net'] = df_out_aws.SRin - df_out_aws.SRout
df_out_carra['SR_net'] = df_out_carra.SRin - df_out_carra.SRout

df_out_carra.loc[df_out_carra.Albedo<0.2, 'Albedo'] = np.nan
df_out_aws.loc[df_out_aws.Albedo<0.2, 'Albedo'] = np.nan

df_out_aws = df_out_aws[~df_out_aws.index.duplicated(keep='first')]
df_out_carra = df_out_carra[~df_out_carra.index.duplicated(keep='first')]
# df_out_carra['Tsurf']  = (df_out_aws.LRout_mdl / 0.98 / 5.67e-08 -(1 - 0.98) * df_out_carra.LRin)**(1/4)  
df_out_carra['Tsurf']  = (df_out_aws.LRout_mdl / 5.67e-08)**(1/4)  
common_idx = df_out_aws.index.intersection(df_out_carra.index)
df_out_aws = df_out_aws.loc[common_idx, :]
df_out_carra = df_out_carra.loc[common_idx, :]
plt.close('all')
#%%
plt.close('all')
for var in ['LRin', 'Tsurf', 'melt_cumul', 'melt_mweq']:
    ME = np.mean(df_out_carra[var] - df_out_aws[var])
    RMSE = np.sqrt(np.mean((df_out_carra[var] - df_out_aws[var])**2))

    fig = plt.figure(figsize=(12, 4))
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    
    # first plot
    df_out_aws[var].plot(ax=ax1, label='AWS')
    df_out_carra[var].plot(ax=ax1, label='CARRA')
    ax1.set_ylabel(var)
    ax1.set_title(station)
    ax1.legend()

    # second plot
    ax2.plot(df_out_aws[var], df_out_carra[var], marker='.',ls='None')
    ax2.set_xlabel('AWS')
    ax2.set_ylabel('CARRA')
    ax2.set_title(var)
    slope, intercept, r_value, p_value, std_err = linregress(df_out_aws[var], df_out_carra[var])
    max_val = max(df_out_aws[var].max(), df_out_carra[var].max())
    min_val = min(df_out_aws[var].min(), df_out_carra[var].min())
    ax2.plot([min_val, max_val], [min_val, max_val], 'k-', label='1:1 Line')
    regression_line = slope * df_out_aws[var] + intercept
    ax2.plot(df_out_aws[var], regression_line, 'r-', label='Linear Regression')
    ax2.legend(loc='lower right')

    
    # Annotate with RMSE and ME
    ax2.annotate(f'RMSE: {RMSE:.2f}\nME: {ME:.2f}', 
                 xy=(0.05, 0.95), xycoords='axes fraction', 
                 horizontalalignment='left', verticalalignment='top',
                 fontsize=10, bbox=dict(boxstyle="round,pad=0.3",
                                        edgecolor='black', facecolor='white'))



# %% Producing old format
# df_old = pd.read_csv('data_KAN_M_combined_hour.txt',sep='\t')
# df_old = pd.DataFrame(columns=df_old.columns)

# df_old['Year'] = df_in.index.year.values
# df_old['MonthOfYear'] = df_in.index.month.values
# df_old['HourOfDayUTC'] = df_in.index.hour.values
# df_old['DayOfYear'] = df_in.index.dayofyear.values
# for v in df_in.columns:
#     if v in df_old.columns:
#         df_old[v] = df_in[v].values
        
# df_old.replace({np.nan: -999}).to_csv('QAS_U_CARRA.txt', sep ='\t', index=None)