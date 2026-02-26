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

def load_data(path, station):
    df_output = xr.open_dataset(path+'/'+station+'_surface.nc').to_dataframe()
    slwc = xr.open_dataset(path+'/'+station+'_slwc.nc')
    df_output['total_slwc'] = slwc.sum(dim='level').slwc.to_series()
    percolation_depth = slwc.depth.where(slwc.slwc>0).min("level").where((slwc.slwc>0).any("level"))
    df_output['percolation_depth'] = percolation_depth.to_series()
    total_snowc = xr.open_dataset(path+'/'+station+'_snowc.nc')
    total_snic = xr.open_dataset(path+'/'+station+'_snic.nc')
    df_output['total_snowc'] = total_snowc.sum(dim='level').snowc.to_series()
    df_output['snowc_1'] = total_snowc.isel(level=0).snowc.to_series()
    df_output['total_snic'] = total_snic.sum(dim='level').snic.to_series()

    # df_in_aws = load_promice_old("QAS_U_CARRA.txt")
    # df_output_1 [ df_in_aws.columns] = df_in_aws.values
    # df_output_1.index = df_output_1.index - pd.Timedelta('1D')
    # del df_in_aws
    df_output['melt_cumul'] = df_output.melt_mweq.cumsum()
    df_output['SR_net'] = df_output.SRin - df_output.SRout

    df_output = df_output[~df_output.index.duplicated(keep='first')]

    return df_output

station_list =  [s.replace('.nc', '') for s in os.listdir("./input/weather data/CARRA_at_AWS/")]
# for station in station_list:
for station in ['DY2']:
    print(station)
    path_1 = f'output/{station}_100_layers_3h_old'
    path_2 = f'output/{station}_100_layers_3h_copy'
    name_1 = 'without copies of arrays'
    name_2 = 'with copies of arrays'
    df_output_1=  load_data(path_1, station)
    df_output_2=  load_data(path_2, station)


    # df_output_2['Tsurf']  = (df_output_1.LRout_mdl / 0.98 / 5.67e-08 -(1 - 0.98) * df_output_2.LRin)**(1/4)
    df_output_2['Tsurf']  = (df_output_1.LRout_mdl / 5.67e-08)**(1/4)
    common_idx = df_output_1.index.intersection(df_output_2.index)
    df_output_1 = df_output_1.loc[common_idx, :]
    df_output_2 = df_output_2.loc[common_idx, :]

    #%%


    var_list = ['smb_mweq', 'refreezing_mweq', 'runoff_mweq', 'total_slwc', 'percolation_depth', 'total_snowc', 'total_snic']

    fig, axes = plt.subplots(len(var_list),2,  figsize=(12, 4 * len(var_list)))
    fig.subplots_adjust(top=0.97, bottom= 0.04, hspace=0.15)

    for i, var in enumerate(var_list):
        ME = np.mean(df_output_2[var] - df_output_1[var])
        RMSE = np.sqrt(np.mean((df_output_2[var] - df_output_1[var]) ** 2))

        ax1 = axes[i, 0] if len(var_list) > 1 else axes[0]
        ax2 = axes[i, 1] if len(var_list) > 1 else axes[1]

        # First plot
        if var in ['smb_mweq', 'refreezing_mweq', 'runoff_mweq',]:
            tmp_1 = df_output_1[var].cumsum()
            tmp_2 = df_output_1[var].cumsum()
            label = var + ' cumulated'
        else:
            tmp_1 = df_output_1[var]
            tmp_2 = df_output_1[var]
            label = var
        tmp_1.plot(ax=ax1, label=name_1)
        tmp_2.plot(ax=ax1, label=name_2)
        ax1.set_ylabel(label)
        ax1.legend(title=var)
        ax1.grid()
        # Second plot
        ax2.plot(df_output_1[var], df_output_2[var], marker='.', ls='None', label=var)
        ax2.set_xlabel(var + ' '+ name_1)
        ax2.set_ylabel(var + ' '+ name_2)
        ax2.grid()

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
    fig.savefig(f'{station}_compare.png',dpi=120)

# %%
    plt.figure()
    df_output_2.snowc_1.plot(marker='o', label='snowc first layer')
    df_output_2.snowbkt.plot(marker='o', label='snow bucket')
    df_output_2.melt_mweq.plot(marker='o', label='melt')
    plt.legend()
