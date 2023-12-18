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
import pandas as pd
import lib.io as io
import lib.plot as lpl
from lib.initialization import ImportConst
from lib.seb_smb_model import HHsubsurf
from os import mkdir
import time

# def run_SEB_firn():

start_time = time.time()
print('start processing')

# importing standard values for constants
c = ImportConst()

c.station = 'KAN_M'
# c.surface_input_path = "./input/weather data/data_KAN_U_2009-2019.txt"
# c.station = 'KAN_M'

# c.surface_input_path = "./input/weather data/data_"+c.station+"_combined_hour.txt"
# c.surface_input_driver = "AWS_old" 

c.surface_input_path = "./input/weather data/CARRA_at_AWS.nc"
c.surface_input_driver = "CARRA" 

# assigning constants specific to this simulation
c.snowthick_ini = 0.1
c.z_max = 50
c.num_lay = 100
# c.lim_new_lay = c.accum_AWS/c.new_lay_frac;

df_in = io.load_surface_input_data(c)

df_in = df_in.loc['2007-09-01':,:]

print('start/end of input file', df_in.index[0], df_in.index[-1])
# DataFrame for the surface is created, indexed with time from df_aws
df_out = pd.DataFrame()
df_out["time"] = df_in.index
df_out = df_out.set_index("time")
# %% 
print('reading inputs took %0.03f sec'%(time.time() -start_time))
start_time = time.time()
# The surface values are received 
(
    df_out["L"],
    df_out["LHF"],
    df_out["SHF"],
    df_out["theta_2m"],
    df_out["q_2m"],
    df_out["ws_10m"],
    df_out["Re"],
    df_out["melt_mweq"],
    df_out["sublimation_mweq"],
    df_out["SRin"],
    df_out["SRout"],
    df_out["LRin"],
    df_out["LRout_mdl"],
    snowc,
    snic,
    slwc,
    T_ice,
    zrfrz,
    rhofirn,
    zsupimp,
    dgrain,
    df_out['zrogl'],
    Tsurf,
    grndc,
    grndd,
    pgrndcapc,
    pgrndhflx,
    dH_comp,
    snowbkt,
    compaction,
    df_out['snowthick']
) = HHsubsurf(df_in, c)

print('\nHHsubsurf took %0.03f sec'%(time.time() -start_time))
start_time = time.time()

thickness_act = snowc * (c.rho_water / rhofirn) + snic * (c.rho_water / c.rho_ice)
depth_act = np.cumsum(thickness_act, 0)
density_bulk = (snowc + snic) / (snowc / rhofirn + snic / c.rho_ice)

# Writing output
c.RunName = c.station + "_" + str(c.num_lay) + "_layers"
i = 0
succeeded = 0
while succeeded == 0:
    try:
        mkdir(c.output_path + c.RunName)
        succeeded = 1
    except:
        if i == 0:
            c.RunName = c.RunName + "_" + str(i)
        else:
            c.RunName = c.RunName[: -len(str(i - 1))] + str(i)
        i = i + 1

# io.write_2d_netcdf(snic, 'snic', depth_act, df_in.index, c)
io.write_2d_netcdf(slwc, 'slwc', depth_act, df_in.index, c)
#io.write_2d_netcdf(rhofirn, 'rhofirn', depth_act, df_in.index, c)
io.write_1d_netcdf(df_out, c)
io.write_2d_netcdf(density_bulk, 'density_bulk', depth_act, df_in.index, c)
io.write_2d_netcdf(T_ice, 'T_ice', depth_act, df_in.index, c)
# io.write_2d_netcdf(rhofirn, 'rho_firn_only', depth_act, df_in.index, RunName)
# io.write_2d_netcdf(rfrz, 'rfrz', depth_act, df_in.index, RunName)
# io.write_2d_netcdf(dgrain, 'dgrain', depth_act, df_in.index, RunName)
# io.write_2d_netcdf(compaction, 'compaction', depth_act, df_in.index, RunName)
df_out['zrfrz_sum'] = zrfrz.sum(axis=0)

print('writing output files took %0.03f sec'%(time.time() -start_time))
start_time = time.time()

# Plot output
plt.close("all")
#lpl.plot_summary(df_in, c, 'input_summary', var_list = ['RelativeHumidity1','RelativeHumidity2'])
# %%
lpl.plot_summary(df_out, c, 'SEB_output')
for var in ['slwc','T_ice','density_bulk']:
    lpl.plot_var(c.station, c.RunName, var, ylim=(10, -5), zero_surf=False)

# LR modelled vs obs
plt.figure()
plt.plot(df_out["LRout_mdl"], 
         df_in.LongwaveRadiationUpWm2,
         marker='.',ls='None')


#%%
plt.figure()
ax=plt.gca()
df_out.melt_mweq.cumsum().plot(ax=ax, label='melt')
df_in.Snowfallmweq.cumsum().plot(ax=ax, label='Snowfall')
df_out.zrogl.cumsum().plot(ax=ax, label='runoff')
df_out.zrfrz_sum.cumsum().plot(ax=ax, label='refreezing')
df_out.snowthick.plot(ax=ax, label='snowthickness')
plt.legend()

# %% Surface height evaluation
# if 'SurfaceHeightm' in df_in.columns:
#     plt.figure()
#     plt.plot(df_out.index, -depth_act[-1,0] + depth_act[-1,:])
#     plt.plot(df_in.index, df_in.SurfaceHeightm)
# else:
path_aws_l4 = 'C:/Users/bav/OneDrive - Geological survey of Denmark and Greenland/Code/PROMICE/PROMICE-AWS-toolbox/out/L4/'
df_obs = pd.read_csv(path_aws_l4+c.station+'_L4.csv')
df_obs.time= pd.to_datetime(df_obs.time)
df_obs = df_obs.set_index('time')

plt.figure()
plt.plot(df_out.index, -depth_act[-1,0] + depth_act[-1,:])
plt.plot(df_obs.index, df_obs.z_surf_combined)
plt.ylabel('Surface height (m)')
plt.title(c.station)

# %% 
if c.station in ['KAN_M']:
    file = 'C:/Users/bav/OneDrive - Geological survey of Denmark and Greenland/Data/SUMup/data/SMB data/to add/SnowFox_GEUS/SF_KAN_M.txt'

    df_sf = pd.read_csv(file,delim_whitespace=True)
    df_sf[df_sf==-999] = np.nan
    df_sf['time'] = pd.to_datetime(df_sf[['Year','Month','Day']])
    df_sf = df_sf.set_index('time')
    df_sf['SWE_mweq'] =df_sf['SWE(cmWeq)']/100

    plt.figure()
    ax=plt.gca()
    df_sf.SWE_mweq.plot(ax=ax, marker='o')
    (df_in.loc['2018-08-12':'2019-05-01'].Snowfallmweq).cumsum().plot(ax=ax, label='Snowfall')
    (df_in.loc['2019-09-01':'2020-05-01'].Snowfallmweq).cumsum().plot(ax=ax, label='Snowfall')

#%%
print('plotting took %0.03f sec'%(time.time() -start_time))
start_time = time.time()


# if __name__ == "__main__":
#     run_SEB_firn()


    

