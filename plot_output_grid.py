# -*- coding: utf-8 -*-
"""
Created on %(date)s
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import netCDF4 as nc
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import matplotlib
import os
from lib.initialization import Struct

def plot_pixel_multifile(pixel=120096, var='LHF'):
    output_folder = 'C:/Users/bav/data_save/output firn model/'
    list_folder = [folder for folder in os.listdir(output_folder) if 'grid' in folder]

    # if isintance(var, str) :
    #     var=[var]
    if var in ['L', 'LHF', 'SHF', 'theta_2m', 'q_2m', 'ws_10m', 'Re', 'melt_mweq', 'sublimation_mweq', 'SRin', 'SRout', 'LRin', 'LRout_mdl', 'zrogl', 'snowbkt', 'snowthick']:
        file_name = 'surface'
    
    
    list_path = []
    for folder in list_folder:
        year = folder.split('_')[1]
        month = folder.split('_')[2]
        list_path.append(
            output_folder+folder+'/pixel_'+str(pixel) +'_' + year + '_' + month +'_100_layers_3H' \
                +'/pixel_'+str(pixel)+'_'+year+'_'+month+'_'+file_name+'.nc')
        
    ds_out = xr.open_mfdataset(list_path)
    
    def new_fig():
        fig, ax = plt.subplots(1, 1, sharex=True, figsize=(15, 10))
        plt.subplots_adjust(
            left=0.1, right=0.9, top=0.95, bottom=0.1, wspace=0.2, hspace=0.05
        )
        return fig, ax

    fig, ax = new_fig()
    ds_out[var].plot(ax=ax)
    fig.suptitle('pixel '+str(pixel))
    

def plot_surface_var_pixel(pixel=120096, filetag="summary", var_list=None):
    output_folder = 'C:/Users/bav/data_save/output firn model/'
    list_folder = sorted([folder for folder in os.listdir(output_folder) if 'grid' in folder])
    
    list_path = []
    for folder in list_folder:
        year = folder.split('_')[1]
        month = folder.split('_')[2]
        list_path.append(
            output_folder+folder+'/pixel_'+str(pixel) +'_' + year + '_' + month +'_100_layers_3H' \
                +'/pixel_'+str(pixel)+'_'+year+'_'+month+'_surface.nc')

    tmp =pd.read_csv(
        output_folder+folder+'/pixel_'+str(pixel) +'_' + year + '_' + month +'_100_layers_3H' \
            +'/constants.csv', dtype={'key':str})
    tmp['value_num'] = pd.to_numeric(tmp.value, errors='coerce')
    tmp.loc[tmp.value_num.notnull(),'value'] = tmp.loc[tmp.value_num.notnull(),'value_num']
    tmp = tmp.set_index('key')[['value']]
    c = Struct(**tmp.to_dict()['value'] )
    c.plot_path = output_folder+'plots_pixel_'+str(pixel)\
        +list_folder[0].replace('grid','')+list_folder[-1].replace('grid','')
    os.makedirs(c.plot_path, exist_ok=True)
        
    df = xr.open_mfdataset(list_path).to_dataframe()
    
    def new_fig():
        fig, ax = plt.subplots(8, 1, sharex=True, figsize=(15, 10))
        plt.subplots_adjust(
            left=0.1, right=0.9, top=0.97, bottom=0.1, wspace=0.2, hspace=0.05
        )
        return fig, ax

    if not var_list:
        var_list = df.columns
        df.columns = df.columns.astype(str)

    fig, ax = new_fig()
    count = 0
    count_fig = 0

    for i, var in enumerate(var_list):
        var = str(var)           
        if "_origin" in var.lower():
            continue
        if var + "_Origin" in df.columns:
            df[var].plot(ax=ax[count], color="k", label="_no_legend_")
            for k in df[var + "_Origin"].unique():
                tmp = df.copy()
                tmp.loc[df[var + "_Origin"] != k, var] = np.nan
                tmp[var].plot(ax=ax[count], label="origin: " + str(int(k)))
                ax[count].legend()
        else:
            df[var].plot(ax=ax[count])

        ax[count].set_ylabel(var)
        ax[count].grid()
        ax[count].set_xlim((df.index[0], df.index[-1]))
        
        if var == "L": ax[count].set_ylim((-30000, 30000))

        count = count + 1

        if (count == len(ax)) & (var != var_list[-1]):
            ax[0].set_title(c.station)
            plt.savefig(
                c.plot_path + "/" + c.station + "_summary_" + str(count_fig),
                bbox_inches="tight",
            )
            count_fig = count_fig + 1
            fig, ax = new_fig()
            count = 0
    if count < len(ax)-2:
        count = count - 1
        ax[count].xaxis.set_tick_params(which="both", labelbottom=True)

        for k in range(count + 1, len(ax)):
            ax[k].set_axis_off()
    ax[0].set_title(c.station)
    
    plt.savefig(
        c.plot_path + "/" + c.station + "_summary_" + str(count_fig),
        bbox_inches="tight",
    )
if __name__ == "__main__":
    plot_surface_var_pixel()
