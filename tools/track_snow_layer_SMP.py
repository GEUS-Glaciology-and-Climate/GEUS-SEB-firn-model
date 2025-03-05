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

import pandas as pd
import os
from matplotlib import cm

os.chdir('..')
import lib.plot as lpl
from lib.initialization import Struct
import lib.io as io

def name_alias(stid):
    rename = {'South Dome':'SDM', 'Saddle':'SDL', 'NASA-U': 'NAU',
                'NASA-E': 'NAE', 'NEEM': 'NEM', 'EastGRIP': 'EGP',
                'DYE-2': 'DY2', 'Tunu-N':'TUN', 'CEN1':'CEN', 'CEN2':'CEN',
                'JAR1':'JAR', 'NASA-SE':'NSE','GITS':'CEN','Humboldt':'HUM',
                # ['Summit', 'DMI'],
                # ['Summit', 'NOAA']
}
    if stid in rename.keys():
        return rename[stid]
    else:
        return stid
# output_path= 'C:/Users/bav/data_save/output firn model/spin up 3H/'
output_path = './output/200_layers/'
run_name = 'EastGRIP_200_layers_3h'

print(run_name)
tmp =pd.read_csv(output_path+'/'+ run_name+'/constants.csv', dtype={'key':str})
# converting all numerical fields to numeric, except station
tmp['value_num'] = pd.to_numeric(tmp.value, errors='coerce')
msk = (tmp.value_num.notnull() & (tmp.key!='station'))
tmp.loc[msk,'value'] = tmp.loc[msk,'value_num']
tmp = tmp.set_index('key')[['value']]
# making it a structure
c = Struct(**tmp.to_dict()['value'] )
c.RunName=run_name

if c.surface_input_driver=='CARRA' and c.zdtime == 3600:
    print('resample')
    resample=True
else:
    resample=False

df_in, c = io.load_surface_input_data(c, resample=resample)
if output_path != c.output_path:
    print('Warning: Output has been moved from',c.output_path,'to',output_path)
    c.output_path = output_path

#  loading surface variables
try:
    df_out = xr.open_dataset(c.output_path+run_name+'/'+c.station+'_surface.nc', decode_cf=True).to_dataframe()
    df_in = df_in.loc[df_out.index[0]:df_out.index[-1],:]
except Exception as e:
    print(c.RunName, e)

# %% Surface height evaluation
# extracting surface height

# path_aws_l4 = '../thredds-data/level_3_sites/csv/hour/'
path_aws_l4 = 'C:/Users/bav/GitHub/PROMICE data/thredds/level_3_sites/hour/'
if os.path.isfile(path_aws_l4+name_alias(c.station)+'_hour.csv'):
    df_obs = pd.read_csv(path_aws_l4+name_alias(c.station)+'_hour.csv')
    obs_avail = True
else:
    path_aws_l4 = '../GC-Net-Level-1-data-processing/L1/hour/'
    if os.path.isfile(path_aws_l4+c.station.replace(' ','')+'.csv'):
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
        print(c.RunName, ': no weather observation was found')
        obs_avail = False

        # return []
    # else:
    #     tmp = pd.DataFrame()
# if len(df_obs)>0:
#     df_obs.time= pd.to_datetime(df_obs.time).dt.tz_convert(None)
#     if len(tmp)>0:
if obs_avail:
    df_obs.time= pd.to_datetime(df_obs.time)
    df_obs = df_obs.set_index('time')
    df_obs = df_obs.resample(pd.infer_freq(df_out.index)).mean()

    fig = plt.figure()
    tmp = (df_obs.z_surf_combined -df_out.surface_height).mean()
    plt.plot(df_obs.index, df_obs.z_surf_combined-tmp,
              marker='.',ls='None', label='AWS')
    plt.plot(df_out.index, df_out.surface_height,color='tab:red',
             label='model')
    plt.legend()
    plt.ylabel('Surface height (m)')
    plt.title(c.station)
    fig.savefig(c.output_path+c.RunName+'/'+c.station+'_surface_height.png', dpi=120)

# %% Tracking layer
site = c.station
output_path = c.output_path
run_name = c.RunName
filename = output_path+"/" + run_name + "/" + site + "_compaction.nc"
ds = xr.open_dataset(filename, decode_cf=True)

filename = output_path+"/" + run_name + "/" + site + "_surface.nc"
ds_surf= xr.open_dataset(filename, decode_cf=True)

compaction = ds["compaction"].data
time = ds["time"].data
depth_act = ds["depth"].data
H_surf = (ds.depth.isel(level=-1)
        -ds.depth.isel(level=-1).isel(time=0)
        -(ds.depth.isel(level=-1).diff(dim='time')
          .where(np.abs(ds.depth.isel(level=-1)
                  .diff(dim='time'))>1,0).cumsum())).data
H_surf = np.insert(H_surf,0,0)
df_ssd = pd.DataFrame(index=pd.to_datetime(time))
df_ssd['H_surf'] = H_surf
# df_ssd = df_ssd.loc[slice(None,'2024-06-01')]
df_ssd['H_surf'] = df_ssd['H_surf'] - df_ssd.loc['2024-06-01', 'H_surf'].mean()
df_obs['z_surf_combined'] = df_obs['z_surf_combined'] - df_obs.loc['2024-06-01', 'z_surf_combined'].mean(
    df_obs['z_surf_combined']=df_obs['z_surf_combined']-0.2
# %%

date_list = [pd.to_datetime(f'2024-01-{day}') for day in range(22,25,1)]
date_list = [pd.to_datetime(d) for d in ['2022-07-17','2023-01-02','2023-01-26',
                                         '2023-05-21',
                                         '2023-06-29','2023-09-02',
                                         '2023-09-20','2023-12-01','2024-02-24']]
for date_start in date_list:
    print(np.datetime64(date_start))
    df_ssd[str(date_start)] = lpl.track_horizon(time, H_surf, depth_act,
                                            compaction,
                                            np.datetime64(date_start),
                                            0, step=1)
    print(df_ssd[str(date_start)].iloc[-1])

# %%

filename = output_path+"/" + run_name + "/" + site + "_snowc.nc"
snowc = xr.open_dataset(filename).transpose()
filename = output_path+"/" + run_name + "/" + site + "_snic.nc"
snic = xr.open_dataset(filename).transpose()
filename = output_path+"/" + run_name + "/" + site + "_rhofirn.nc"
rhofirn = xr.open_dataset(filename).transpose()
ds['density'] = (snowc.snowc + snic.snic) / (snowc.snowc / rhofirn.rhofirn + snic.snic / 900)
ds['SWE'] = (snowc.snowc + snic.snic)


# %%
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


plt.close('all')
fig = plt.figure(figsize=(10,6))
gs = gridspec.GridSpec(1, 3, width_ratios=[1/2, 1/4, 1/4], wspace=0.01)

ax1 = fig.add_subplot(gs[0])
# first panel
# ax1.plot(df_obs.t_u.resample('D').mean())
# ax1.set_ylim(-140,10)
(ds_surf.melt_mweq.resample(time='D').sum()*1000).plot(ax=ax1,c='k',lw=2)
# ax1.set_ylim(-1,3)
ax1.set_ylabel('Daily melt (mm w.e.)                  ',ha='right', va='bottom')

ax=plt.twinx(ax1)

cmap = cm.get_cmap('tab10')

ax.plot(df_ssd.H_surf, label="Surface (model)")
ax.plot(df_obs.z_surf_combined, lw=3, alpha=1, label="Surface (observation)")

for i, d in  enumerate(date_list):
    ax.plot(df_ssd.index,- df_ssd[str(d)].values + df_ssd.H_surf,
            color='k',
            ls='--',
            alpha=0.6,
            label="_no_legend_")
ax.plot(np.nan, np.nan,'k--', alpha=0.6, label="tracked layers")
ax.legend(loc='lower center', bbox_to_anchor=(0.5,0.99))
ax.grid()
ax.set_ylabel("Height (m)")
ax.set_ylim(-3,0)
ax.set_xlim([pd.to_datetime(t) for t in ['2022-06-01', '2024-06-01']])

ax2 = fig.add_subplot(gs[1])
ds_sel = ds.sel(time='2024-06-01T00:00:00', method='nearest')

depth = -ds_sel.depth.values  # Ensure depth is negative (increasing downwards)
density = ds_sel.density.values

# Prepend depth=0 for the first step
depth_steps = np.insert(depth, 0, 0)
density_steps = np.insert(density, 0, density[0])

ax2.step(density_steps, depth_steps, where='pre')

ax2.set_ylim(-3, 0)
ax2.set_xlim(300, 400)
ax2.set_xlabel('Density (kg m⁻³)')
ax2.set_ylabel('')
ax2.grid(True, alpha=0.5)

ax3 = fig.add_subplot(gs[2])
ds_sel = ds.sel(time='2024-06-01T00:00:00', method='nearest')

depth = -ds_sel.depth.values  # Ensure depth is negative (increasing downwards)
cumulative_swe = (ds_sel.SWE * 1000).cumsum(dim='level').values  # Convert to mm

# Prepend depth=0 for the first step
depth_steps = np.insert(depth, 0, 0)
cumulative_swe_steps = np.insert(cumulative_swe, 0, cumulative_swe[0])

ax3.step(cumulative_swe_steps, depth_steps, where='pre')

ax3.set_xlabel('Cumulative SWE (mm)')
ax3.set_ylabel('Height relative to 2022-06-01 surface (m)')
ax3.set_ylim(-3, 0)
ax3.set_xlim(0, 1100)
ax3.grid(True, alpha=0.5)
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position("right")
ax.set_yticklabels([])
ax2.set_yticklabels([])

fig.suptitle(site)
fig.savefig('side analysis/layer_burial_eastGRIP.png', dpi=300)
