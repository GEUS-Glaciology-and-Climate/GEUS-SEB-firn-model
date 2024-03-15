# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import xarray as xr
import lib.plot as lpl
from lib.initialization import Struct
import pandas as pd
import os
import lib.io as io

station_list_list = [
                # ['KAN_L','KAN_M','KAN_U'],
                # ['KPC_U','KPC_L'],
                # ['QAS_L','QAS_M','QAS_A','QAS_U'],
                # ['UPE_L','UPE_U'],
                # ['THU_L','THU_U'],
                ['SWC','JAR','JAR2','JAR3', 'SMS1','SMS2','SMS3','SMS4']]

# plt.close('all')
for station_list in station_list_list:
    fig = plt.figure()
    
    for stid in station_list:
        output_path = './output/'
        run_name = stid+'_100_layers_3H'
        
        print(run_name)
        try: 
            tmp =pd.read_csv(output_path+'/'+ run_name+'/constants.csv', dtype={'key':str})
        except:
            continue
        tmp['value_num'] = pd.to_numeric(tmp.value, errors='coerce')
        tmp.loc[tmp.value_num.notnull(),'value'] = tmp.loc[tmp.value_num.notnull(),'value_num']
        tmp = tmp.set_index('key')[['value']]
        c = Struct(**tmp.to_dict()['value'] )
        c.RunName=run_name
            
        #  loading surface variables
        df_out = xr.open_dataset(c.output_path+run_name+'/'+c.station+'_surface.nc').to_dataframe()

    
        if 'surface_height' not in df_out.columns:
            filename = c.output_path + run_name + "/" + c.station + "_T_ice.nc"
            ds = xr.open_dataset(filename).transpose()
            ds['surface_height'] = (ds.depth.isel(level=-1)
                    -ds.depth.isel(level=-1).isel(time=0)
                    -(ds.depth.isel(level=-1).diff(dim='time')
                      .where(ds.depth.isel(level=-1)
                             .diff(dim='time')>6,0).cumsum()))
            ds['surface_height'].values[0] = 0
            df_out['surface_height'] = ds['surface_height'].values
            del ds
    
        
        path_aws_l4 = '../PROMICE/PROMICE-AWS-toolbox/out/L4/'
        if os.path.isfile(path_aws_l4+c.station+'_L4_ext.csv'):
            df_obs = pd.read_csv(path_aws_l4+c.station+'_L4_ext.csv')
        elif os.path.isfile(path_aws_l4+c.station+'_L4.csv'):
            df_obs = pd.read_csv(path_aws_l4+c.station+'_L4.csv')
        else:
            path_aws_l4 = '../PROMICE/GC-Net-Level-1-data-processing/L1/daily/'
            if os.path.isfile(path_aws_l4+c.station+'_daily.csv'):
                import nead
                df_obs = nead.read(path_aws_l4+c.station.replace(' ','')+'_daily.csv').to_dataframe()
                if 'HS_combined' not in df_obs.columns:
                    df_obs['HS_combined'] = df_obs.HS1
                df_obs = df_obs.rename(columns={'timestamp':'time',
                                                'HS_combined':'z_surf_combined',
                                                'T10m': 't_i_10m',
                                                'LHF':'dlhf_u',
                                                'SHF': 'dshf_u',
                                                'OLWR':'LRout',
                                                'Tsurf':'t_surf',
                                                })
            else:
                print('No observation for', c.station)
                
        df_obs.time= pd.to_datetime(df_obs.time).dt.tz_convert(None)
        df_obs = df_obs.set_index('time')
        # df_obs = df_obs.resample('M').asfreq()
        
        tmp = (df_obs.z_surf_combined -df_out.surface_height).mean()
        p = plt.plot(df_obs.index, df_obs.z_surf_combined-tmp, 
                  alpha=0.7, lw=4, label=stid)
        plt.plot(df_out.index, df_out.surface_height,color=p[-1].get_color(),
                 label='__nolegend__')
    p = plt.plot(df_obs.index, df_obs.z_surf_combined*np.nan, 
              alpha=0.7, lw=4, c='k', label='observation')
    plt.plot(df_out.index, df_out.surface_height*np.nan,color=p[-1].get_color(),
             label='GEUS SEB firn model')
    plt.legend()
    plt.ylabel('Surface height (m)')
