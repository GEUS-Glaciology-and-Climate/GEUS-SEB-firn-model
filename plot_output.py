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

run_name = 'QAS_U_100_layers_2'
station= 'QAS_U'

c = Struct(**pd.read_csv('output/'+run_name+'/constants.csv',
                         dtype={'key':str})
           .set_index('key').to_dict()['value'] )
c.RunName=run_name
df_in = io.load_surface_input_data(c)
df_out = xr.open_dataset('output/'+run_name+'/'+station+'_surface.nc').to_dataframe()

df_in = df_in.loc[df_out.index[0]:df_out.index[-1],:]


lpl.plot_summary(df_out, c, 'SEB_output')
for var in ['slwc','T_ice','density_bulk']:
    lpl.plot_var(c.station, c.RunName, var, zero_surf=False)

# extracting surface height
filename = "Output/" + run_name + "/" + station + "_T_ice.nc"
ds = xr.open_dataset(filename).transpose()
df_out['surface_height'] = ds.depth.isel(level=-1) - ds.depth.isel(level=-1).isel(time=0)
del ds
    
#%%
plt.figure()
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

plt.figure()
plt.plot(df_out.index, df_out.surface_height, label='model')
plt.plot(df_obs.index, df_obs.z_surf_combined, label='AWS')
plt.legend()
plt.ylabel('Surface height (m)')
plt.title(c.station)

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
    