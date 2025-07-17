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
import os
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
station= 'QAS_U'

if __name__ == "__main__":
    os.chdir('..')


station_list =  [s.replace('.nc', '') for s in os.listdir("./input/weather data/CARRA_at_AWS/")]
for station in station_list:
    print(station)
    path_1 = f'output/HH precipitation/{station}_100_layers_3h'
    path_2 = f'output/new/{station}_100_layers_3h'
    name_1 = 'HH snow/rain transition'
    name_2 = 'snow/rain transition at 0 degC'
    if not os.path.isfile(path_1+'/'+station+'_surface.nc'):
        continue
    if not os.path.isfile(path_2+'/'+station+'_surface.nc'):
        continue
    df_output_1 = xr.open_dataset(path_1+'/'+station+'_surface.nc').to_dataframe()

    # df_in_aws = load_promice_old("QAS_U_CARRA.txt")
    # df_output_1 [ df_in_aws.columns] = df_in_aws.values
    # df_output_1.index = df_output_1.index - pd.Timedelta('1D')
    # del df_in_aws

    df_output_2 = xr.open_dataset(path_2+'/'+station+'_surface.nc').to_dataframe()
    # df_output_2.index = df_output_2.index.round('H')
    # df_in_carra = load_CARRA_data("./input/weather data/CARRA_at_AWS.nc", station)
    # df_output_2 [ df_in_carra.columns] = df_in_carra
    # del df_in_carra

    df_output_1['melt_cumul'] = df_output_1.melt_mweq.cumsum()
    df_output_2['melt_cumul'] = df_output_2.melt_mweq.cumsum()
    # df_output_1['runoff_cumul'] = df_output_1.zrogl.cumsum()
    # df_output_2['runoff_cumul'] = df_output_2.zrogl.cumsum()
    df_output_1['SR_net'] = df_output_1.SRin - df_output_1.SRout
    df_output_2['SR_net'] = df_output_2.SRin - df_output_2.SRout

    # df_output_2.loc[df_output_2.Albedo<0.2, 'Albedo'] = np.nan
    # df_output_1.loc[df_output_1.Albedo<0.2, 'Albedo'] = np.nan

    df_output_1 = df_output_1[~df_output_1.index.duplicated(keep='first')]
    df_output_2 = df_output_2[~df_output_2.index.duplicated(keep='first')]
    # df_output_2['Tsurf']  = (df_output_1.LRout_mdl / 0.98 / 5.67e-08 -(1 - 0.98) * df_output_2.LRin)**(1/4)
    df_output_2['Tsurf']  = (df_output_1.LRout_mdl / 5.67e-08)**(1/4)
    common_idx = df_output_1.index.intersection(df_output_2.index)
    df_output_1 = df_output_1.loc[common_idx, :]
    df_output_2 = df_output_2.loc[common_idx, :]

    #%%


    var_list = ['snowfall_mweq', 'smb_mweq']

    fig, axes = plt.subplots(len(var_list),2,  figsize=(12, 4 * len(var_list)))

    for i, var in enumerate(var_list):
        ME = np.mean(df_output_2[var] - df_output_1[var])
        RMSE = np.sqrt(np.mean((df_output_2[var] - df_output_1[var]) ** 2))

        ax1 = axes[i, 0] if len(var_list) > 1 else axes[0]
        ax2 = axes[i, 1] if len(var_list) > 1 else axes[1]

        # First plot
        df_output_1[var].cumsum().plot(ax=ax1, label=name_1)
        df_output_2[var].cumsum().plot(ax=ax1, label=name_2)
        ax1.set_ylabel(var)
        ax1.legend()

        # Second plot
        ax2.plot(df_output_1[var], df_output_2[var], marker='.', ls='None')
        ax2.set_xlabel(name_1)
        ax2.set_ylabel(name_2)
        ax2.set_title(var)

        # slope, intercept, r_value, p_value, std_err = linregress(df_output_1[var], df_output_2[var])
        max_val = max(df_output_1[var].max(), df_output_2[var].max())
        min_val = min(df_output_1[var].min(), df_output_2[var].min())

        ax2.plot([min_val, max_val], [min_val, max_val], 'k-', label='1:1 Line')
        # regression_line = slope * df_output_1[var] + intercept
        # ax2.plot(df_output_1[var], regression_line, 'r-', label='Linear Regression')
        ax2.legend(loc='lower right')

        # Annotate with RMSE and ME
        ax2.annotate(f'RMSE: {RMSE:.2f}\nME: {ME:.2f}',
                     xy=(0.05, 0.95), xycoords='axes fraction',
                     horizontalalignment='left', verticalalignment='top',
                     fontsize=10, bbox=dict(boxstyle="round,pad=0.3",
                                            edgecolor='black', facecolor='white'))

    plt.suptitle(station)
    plt.tight_layout()
    plt.show()
    fig.savefig(f'{station}_snowfall_scheme.png',dpi=200)
