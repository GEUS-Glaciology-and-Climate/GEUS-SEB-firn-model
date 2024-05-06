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
import os
os.chdir('..')
import lib.plot as lpl
from lib.initialization import Struct
import pandas as pd
import os
import lib.io as io

def name_alias(stid):
    rename = {'South Dome':'SDM', 'Saddle':'SDL', 'NASA-U': 'NAU',
                'NASA-E': 'NAE', 'NEEM': 'NEM', 'EastGRIP': 'EGP',
                'DYE-2': 'DY2', 'Tunu-N':'TUN', 'Humboldt':'HUM', 
                # ['Summit', 'DMI'],
                # ['Summit', 'NOAA']
}
    if stid in rename.keys():
        return rename[stid]
    else:
        return stid

# output_path= 'C:/Users/bav/data_save/output firn model/spin up 3H/'
output_path = './output/'
for station in ['Humboldt','EastGRIP', 'KAN_U', 'CEN2', 'CP1','South Dome', 'Tunu-N','NASA-E','NASA-SE','NASA-U','Saddle']:
    run_name = station+'_100_layers_3H'
    
    print(run_name)
    tmp =pd.read_csv(output_path+'/'+ run_name+'/constants.csv', dtype={'key':str})
    tmp['value_num'] = pd.to_numeric(tmp.value, errors='coerce')
    tmp.loc[tmp.value_num.notnull(),'value'] = tmp.loc[tmp.value_num.notnull(),'value_num']
    tmp = tmp.set_index('key')[['value']]
    c = Struct(**tmp.to_dict()['value'] )
    c.RunName=run_name
    if c.surface_input_driver=='CARRA' and c.zdtime == 3600:
        resample=True
    else:
        resample=False
        
    df_in, c = io.load_surface_input_data(c, resample=resample)
    if output_path != c.output_path:
        print('Warning: Output has been moved from',c.output_path,'to',output_path)
        c.output_path = output_path
        
    # try:
    lpl.find_summer_surface_depths(c)

    
    # %% Surface height evaluation
    # extracting surface height
output_path = './output/'
for station in ['DYE-2']:  #['Humboldt','EastGRIP', 'KAN_U', 'CEN2', 'CP1','South Dome', 'Tunu-N','NASA-E','NASA-SE','NASA-U','Saddle']:
    run_name = station+'_100_layers_3H'
    
    print(run_name)
    tmp =pd.read_csv(output_path+'/'+ run_name+'/constants.csv', dtype={'key':str})
    tmp['value_num'] = pd.to_numeric(tmp.value, errors='coerce')
    tmp.loc[tmp.value_num.notnull(),'value'] = tmp.loc[tmp.value_num.notnull(),'value_num']
    tmp = tmp.set_index('key')[['value']]
    c = Struct(**tmp.to_dict()['value'] )
    c.RunName=run_name
    path_aws_l4 = '../PROMICE/PROMICE-AWS-toolbox/out/L4/'
    if os.path.isfile(path_aws_l4+name_alias(c.station)+'_L4_ext.csv'):
        df_obs = pd.read_csv(path_aws_l4+name_alias(c.station)+'_L4_ext.csv')
        obs_avail = True
    elif os.path.isfile(path_aws_l4+name_alias(c.station)+'_L4.csv'):
        df_obs = pd.read_csv(path_aws_l4+name_alias(c.station)+'_L4.csv')
        obs_avail = True
    else:
        path_aws_l4 = '../PROMICE/GC-Net-Level-1-data-processing/L1/daily/'
        if os.path.isfile(path_aws_l4+c.station.replace(' ','')+'_daily.csv'):
            import nead
            df_obs = nead.read(path_aws_l4+c.station.replace(' ','')+'_daily.csv').to_dataframe()
            df_obs = df_obs.rename(columns={'timestamp':'time',
                                            'HS_combined':'z_surf_combined',
                                            'T10m': 't_i_10m',
                                            'LHF':'dlhf_u',
                                            'SHF': 'dshf_u',
                                            'OLWR':'LRout',
                                            'Tsurf':'t_surf',
                                            })
            obs_avail = True
        else:
            print('No observation for', c.station)
            obs_avail = False

    df_out = xr.open_dataset(c.output_path+run_name+'/'+c.station+'_surface.nc').to_dataframe()

    
    from datetime import datetime as dt
    import time
    from scipy import stats
    
    def toYearFraction(date):
        def sinceEpoch(date): # returns seconds since epoch
            return time.mktime(date.timetuple())
        s = sinceEpoch
    
        year = date.year
        startOfThisYear = dt(year=year, month=1, day=1)
        startOfNextYear = dt(year=year+1, month=1, day=1)
    
        yearElapsed = s(date) - s(startOfThisYear)
        yearDuration = s(startOfNextYear) - s(startOfThisYear)
        fraction = yearElapsed/yearDuration
    
        return date.year + fraction
    

    if obs_avail:
        df_obs.time= pd.to_datetime(df_obs.time).dt.tz_convert(None)
        df_obs = df_obs.set_index('time')
        # df_obs = df_obs.resample(pd.infer_freq(df_out.index)).mean()
        df_obs = df_obs[['z_surf_combined']].dropna()
        fig = plt.figure()
        
        slope_obs, _, _, _, _ = stats.linregress(
            (df_obs.index-df_obs.index[0]).total_seconds(), 
                                                 df_obs['z_surf_combined'])
        slope_mod, _, _, _, _ = stats.linregress(
            (df_out.index-df_out.index[0]).total_seconds(), 
                                                 df_out['surface_height'])
        slope_obs = slope_obs*365*24*60*60
        slope_mod = slope_mod*365*24*60*60

        tmp = (df_obs.z_surf_combined -df_out.surface_height).mean()
        plt.plot(df_obs.index, df_obs.z_surf_combined-tmp, 
                  marker='.',ls='None', label='AWS: %0.2f m a-1'%slope_obs)
        plt.plot(df_out.index, df_out.surface_height,color='tab:red',
                 label='Model: %0.2f m a-1'%slope_mod)
        plt.legend()
        plt.grid()
        plt.ylabel('Surface height (m)')
        plt.title(c.station)
        fig.savefig(c.output_path+c.RunName+'/'+c.station+'_surface_height.png', dpi=120)