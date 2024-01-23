# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import lib.plot as lpl
from lib.initialization import Struct
import pandas as pd
import lib.io as io

def main(output_path= 'C:/Users/bav/data_save/output firn model/',
         run_name = 'UPE_U_100_layers_SU',
          station= 'UPE_U'):
    # %%
    run_name = 'UPE_U_100_layers_SU'
    station= 'UPE_U'

    c = Struct(**pd.read_csv(output_path+run_name+'/constants.csv',
                             dtype={'key':str})
               .set_index('key').to_dict()['value'] )
    c.RunName=run_name
    df_in = io.load_surface_input_data(c)
    try:
        df_out = xr.open_dataset(c.output_path+run_name+'/'+station+'_surface.nc').to_dataframe()
        df_in = df_in.loc[df_out.index[0]:df_out.index[-1],:]
        lpl.plot_summary(df_out, c, 'SEB_output')
    except Exception as e:
        print(e)
    for var in ['T_ice','density_bulk','slwc']:
        lpl.plot_var(c.station, c.RunName, var, zero_surf=False)
    #%%
    lpl.plot_var_start_end(c.station, c.RunName, 'T_ice')
    # extracting surface height
    filename = c.output_path + run_name + "/" + station + "_T_ice.nc"
    ds = xr.open_dataset(filename).transpose()
    df_out['surface_height'] = ds.depth.isel(level=-1) - ds.depth.isel(level=-1).isel(time=0)
    del ds
        
    #%%
    fig = plt.figure()
    ax=plt.gca()
    df_out.melt_mweq.cumsum().plot(ax=ax, label='melt')
    df_in.Snowfallmweq.cumsum().plot(ax=ax, label='Snowfall')
    df_out.zrogl.cumsum().plot(ax=ax, label='runoff')
    try:
        df_out.zrfrz_sum.cumsum().plot(ax=ax, label='refreezing')
    except:
        pass
    df_out.snowthick.plot(ax=ax, label='snowthickness')
    plt.legend()
    fig.savefig(c.output_path+c.RunName+'/SMB.png', dpi=120)
    
    # %% Surface height evaluation
    # if 'SurfaceHeightm' in df_in.columns:
    #     plt.figure()
    #     plt.plot(df_out.index, -depth_act[-1,0] + depth_act[-1,:])
    #     plt.plot(df_in.index, df_in.SurfaceHeightm)
    # else:
    path_aws_l4 = '../PROMICE/PROMICE-AWS-toolbox/out/L4/'
    df_obs = pd.read_csv(path_aws_l4+c.station+'_L4.csv')
    df_obs.time= pd.to_datetime(df_obs.time)
    df_obs = df_obs.set_index('time')
    
    fig = plt.figure()
    plt.plot(df_out.index, df_out.surface_height, label='model')
    plt.plot(df_obs.index, df_obs.z_surf_combined, label='AWS')
    plt.legend()
    plt.ylabel('Surface height (m)')
    plt.title(c.station)
    fig.savefig(c.output_path+c.RunName+'/surface_height.png', dpi=120)
    
    # %% 
    lpl.plot_movie(station, run_name, 'T_ice')
    lpl.plot_movie(station, run_name, 'density_bulk')

    
    # %% 
    if c.station in ['KAN_M', 'QAS_M', 'QAS_U','TAS_A','THU_U2']:
        file = '../../Data/SUMup/data/SMB data/to add/SnowFox_GEUS/SF_'+c.station+'.txt'
    
        df_sf = pd.read_csv(file,delim_whitespace=True)
        df_sf[df_sf==-999] = np.nan
        df_sf['time'] = pd.to_datetime(df_sf[['Year','Month','Day']])
        df_sf = df_sf.set_index('time')
        df_sf['SWE_mweq'] =df_sf['SWE(cmWeq)']/100
    
        fig = plt.figure()
        ax=plt.gca()
        df_sf.SWE_mweq.plot(ax=ax, marker='o')
        (df_in.loc['2018-08-12':'2019-05-01'].Snowfallmweq).cumsum().plot(ax=ax, label='Snowfall')
        (df_in.loc['2019-09-01':'2020-05-01'].Snowfallmweq).cumsum().plot(ax=ax, label='Snowfall')
        plt.title(c.station)
        plt.ylabel('Snow accumulation (m w.e.)')
        fig.savefig(c.output_path+c.RunName+'/snowfox_eval.png', dpi=120)

# %%
        
if __name__ == "__main__":
    main()
