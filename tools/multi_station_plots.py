# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import os
os.chdir('..')

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import xarray as xr
import lib.plot as lpl
from lib.initialization import Struct
import pandas as pd
import lib.io as io
from lib.plot import interpolate_temperature_fast,load_sumup
import matplotlib.dates as mdate
locator = mdate.YearLocator(10)

#%% Surface height at ablation sites

station_list_list = [
                ['KAN_L','KAN_M','KAN_U'],
                ['SWC','JAR','JAR2','JAR3', 'SMS1','SMS2','SMS3','SMS4'],
                ['UPE_L','UPE_U'],
                ['THU_L','THU_U2'],
                ['KPC_U','KPC_L'],
                ['SCO_L','SCO_U'],
                ['TAS_L','TAS_U','TAS_A'],
                ['QAS_L','QAS_M','QAS_A','QAS_U'],
                ['NUK_L','NUK_U','NUK_A'],
                ]

# plt.close('all')
fig, ax_list = plt.subplots(3,3, sharex=True, sharey=True, figsize=(10,10))
plt.subplots_adjust(top=0.98,hspace=0, wspace=0, right=0.98)
ax_list=ax_list.flatten()
for station_list, ax in zip(station_list_list, ax_list):
    for stid in station_list:
        output_path = './output/HH precipitation/'
        run_name = stid+'_100_layers_3h'
        
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
        df_out = df_out.resample('W').first()
    
        
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
        df_obs = df_obs[~df_obs.index.duplicated(keep='first')]
        df_obs = df_obs.resample('W').first()
        
        tmp = (df_obs.z_surf_combined -df_out.surface_height).mean()
        p = ax.plot(df_obs.index, df_obs.z_surf_combined-tmp, 
                  alpha=0.7, lw=4, label=stid)
        ax.plot(df_out.index, df_out.surface_height,color=p[-1].get_color(),
                 label='__nolegend__')
    # p = ax.plot(df_obs.index, df_obs.z_surf_combined*np.nan, 
    #           alpha=0.7, lw=4, c='k', label='observation')
    # ax.plot(df_out.index, df_out.surface_height*np.nan,color=p[-1].get_color(),
    #          label='GEUS SEB firn model')
    ax.legend()
    ax.xaxis.set_major_locator(locator)
    ax.grid('both')

fig.text(0.5, 0.04, 'Years', ha='center')
fig.text(0.04, 0.5, 'Surface height (m)', va='center', rotation='vertical')
fig.savefig('side analysis/surface_height_SMB_vs_obs_ablation_HH.png', dpi=300)

#%% Surface height at accumulation sites

station_list = ['DYE-2','CP1','NASA-U','CEN2','Humboldt','EastGRIP',
                     'Tunu-N','NASA-E','NASA-SE','Saddle','South Dome',
                     'Summit']
aliases = {'DYE-2':'DY2','CP1':'CP1','NASA-U':'NAU','CEN2':'CEN2','Humboldt':'HUM','EastGRIP':'EGP',
                     'Tunu-N':'TUN','NASA-E':'NAE','NASA-SE':'NSE','Saddle':'SDL','South Dome':'SDM',
                     'Summit':'Summit'}

# plt.close('all')
fig, ax_list = plt.subplots(3,4, sharex=True, sharey=True, figsize=(10,10))
plt.subplots_adjust(top=0.98,hspace=0, wspace=0, right=0.98)
ax_list=ax_list.flatten()
for stid, ax in zip(station_list, ax_list):
    output_path = './output/new/'
    run_name = stid+'_100_layers_3h'
    
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
    df_out = df_out.resample('W').first()

    try: 
        path_aws_l4 = '../PROMICE/PROMICE-AWS-toolbox/out/L4/'
        if os.path.isfile(path_aws_l4+aliases[c.station]+'_L4_ext.csv'):
            df_obs = pd.read_csv(path_aws_l4+aliases[c.station]+'_L4_ext.csv')
        elif os.path.isfile(path_aws_l4+aliases[c.station]+'_L4.csv'):
            df_obs = pd.read_csv(path_aws_l4+aliases[c.station]+'_L4.csv')
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
        df_obs = df_obs[~df_obs.index.duplicated(keep='first')]
        df_obs = df_obs.resample('W').first()
        
        tmp = (df_obs.z_surf_combined -df_out.surface_height).mean()
        p = ax.plot(df_obs.index, df_obs.z_surf_combined-tmp, 
                  alpha=0.7, lw=4, label='observed')
        ax.plot(df_out.index, df_out.surface_height,color=p[-1].get_color(),
                 label='modelled')
    
        ax.legend(title=stid, loc='upper left')
        ax.xaxis.set_major_locator(locator)
        ax.grid('both')
    except:
        pass

fig.text(0.5, 0.04, 'Years', ha='center')
fig.text(0.04, 0.5, 'Surface height (m)', va='center', rotation='vertical')
fig.savefig('side analysis/surface_height_SMB_vs_obs_accumulation_new.png', dpi=300)

#%% SMB at all sites
compute_smb_evaluation = True

if compute_smb_evaluation:
    with xr.open_dataset("./input/weather data/CARRA_at_AWS.nc") as ds:
        ds_aws = ds.load()
        station_list = ds_aws.stid.load().values
    unwanted= ['NUK_K', 'MIT', 'ZAK_A', 'ZAK_L', 'ZAK_Lv3', 'ZAK_U', 'ZAK_Uv3', # local glaciers
                'LYN_L', 'LYN_T', 'FRE',  # local glaciers
                'KAN_B', 'NUK_B','WEG_B', # off-ice AWS
                'DY2','NSE','SDL','NAU','NAE','SDM','TUN','HUM','SWC', 'JAR', 'NEM',  # redundant
                'T2_08','S10','Swiss Camp 10m','SW2','SW4','Crawford Point 1', 'G1','EGP', # redundant
                'CEN1','THU_L2'
                ] 
    station_list = [s for s in station_list if s not in unwanted]
    station_list = [s for s in station_list if 'v3' not in s]
        
    ds_aws=ds_aws.where(~ds_aws.stid.isin(unwanted), drop=True)
    
    output_path= './output/HH precipitation/'

    
    df_all_smb = pd.DataFrame()
    for stid in station_list:
        run_name = stid + '_100_layers_3h'
        if not os.path.isfile(output_path+'/'+ run_name+'/constants.csv'):
            print('!!!!',stid,'run not finsished !!!!')
            continue

        tmp =pd.read_csv(output_path+'/'+ run_name+'/constants.csv', dtype={'key':str})
        tmp['value_num'] = pd.to_numeric(tmp.value, errors='coerce')
        tmp.loc[tmp.value_num.notnull(),'value'] = tmp.loc[tmp.value_num.notnull(),'value_num']
        tmp = tmp.set_index('key')[['value']]
        c = Struct(**tmp.to_dict()['value'] )
        print(c.station)
        df_out = xr.open_dataset(c.output_path+c.RunName+'/'+c.station+'_surface.nc').to_dataframe()
        # evaluate_smb_sumup(df_out, c)
        
        df_sumup_selec, df_meta_selec = load_sumup(var='SMB', name_var='name', c=c)
        if len(df_sumup_selec) == 0: continue

        msk = df_sumup_selec.start_date.isnull()
        df_sumup_selec.loc[msk, 'start_date'] = pd.to_datetime(df_sumup_selec.loc[msk, 'start_year'].astype(int).astype(str)+'-01-01')
        msk = df_sumup_selec.end_date.isnull()
        df_sumup_selec.loc[msk, 'end_date'] = pd.to_datetime(df_sumup_selec.loc[msk, 'end_year'].astype(int).astype(str)+'-01-01')
        msk = df_sumup_selec.start_date == df_sumup_selec.end_date
        df_sumup_selec.loc[msk, 'end_date'] = pd.to_datetime((df_sumup_selec.loc[msk, 'end_year']+1).astype(int).astype(str)+'-01-01')
        df_sumup_selec['smb_mod'] = np.nan
        
        for i in df_sumup_selec.index:
            df_sumup_selec.loc[i, 'smb_mod'] = df_out.loc[
                df_sumup_selec.loc[i,'start_date']:df_sumup_selec.loc[i,'end_date'],
                'smb_mweq'].sum()
        df_sumup_selec['station'] = c.station
        df_all_smb = pd.concat((df_all_smb,df_sumup_selec.reset_index()), ignore_index=True)
    df_all_smb.to_csv('side analysis/SMB_evaluation.csv')

df_all_smb = pd.read_csv('side analysis/SMB_evaluation.csv')
df_all_smb['start_date'] = pd.to_datetime(df_all_smb['start_date'])
df_all_smb['end_date'] = pd.to_datetime(df_all_smb['end_date'])
fig, ax = plt.subplots(1,2,figsize=(9,9))
plt.subplots_adjust(bottom=0.6,left=0.1,right=0.95,top=0.95)
for stid in ds_aws.stid.data[1:]:
    df_sumup_selec = df_all_smb.loc[df_all_smb.station==stid,:]
    msk = (df_sumup_selec.end_date-df_sumup_selec.start_date) <= pd.to_timedelta('7 day')
    if msk.any():
        df_sumup_daily = df_sumup_selec.loc[msk]        
        ax[0].plot(df_sumup_daily.smb, df_sumup_daily.smb_mod, marker='.', ls='None',label=stid)

    df_sumup_annual = df_sumup_selec.loc[~msk]
    ax[1].plot(df_sumup_annual.smb, df_sumup_annual.smb_mod, marker='.', ls='None',label=stid)
lim1=[-0.2, 0.15]
ax[0].plot(lim1,lim1, color='k')
ax[0].set_xlim(lim1)
ax[0].set_ylim(lim1)
lim2= [-8, 2.5]
ax[1].plot(lim2,lim2, color='k')
ax[1].set_xlim(lim2)
ax[1].set_ylim(lim2)
ax[1].legend(loc='lower center',ncol=5, bbox_to_anchor=(-0.2,-1.6))
ax[0].set_title('Daily/weekly measurements', loc='left')
ax[1].set_title('Seasonal/annual measurements', loc='left')
ax[0].set_xlabel('Observed SMB (m w.e.)')
ax[0].set_ylabel('Modelled SMB (m w.e.)')
ax[1].set_xlabel('Observed SMB (m w.e.)')
ax[1].set_ylabel('Modelled SMB (m w.e.)')
fig.savefig('side analysis/SMB_evaluation_SUMup2024_HH.png', dpi=300)

#%% 10 m temperature at all sites
import lib.plot as lplt
output_path= './output'

compute_T10m_evaluation = True

if compute_T10m_evaluation:
    with xr.open_dataset("./input/weather data/CARRA_at_AWS.nc") as ds:
        ds_aws = ds.load()
        station_list = ds.stid.load().values
    unwanted= ['NUK_K', 'MIT', 'ZAK_A', 'ZAK_L', 'ZAK_Lv3', 'ZAK_U', 'ZAK_Uv3', # local glaciers
                'LYN_L', 'LYN_T', 'FRE',  # local glaciers
                'KAN_B', 'NUK_B','WEG_B', # off-ice AWS
                'DY2','NSE','SDL','NAU','NAE','SDM','TUN','HUM','SWC', 'JAR', 'NEM',  # redundant
                'T2_08','S10','Swiss Camp 10m','SW2','SW4','Crawford Point 1', 'G1','EGP', # redundant
                'CEN1','THU_L2'
                ] 
    station_list = [s for s in station_list if s not in unwanted]
    station_list = [s for s in station_list if 'v3' not in s]
    ds_aws=ds_aws.where(ds_aws.stid.isin(station_list), drop=True)
    
    df_sumup, df_meta = lplt.load_sumup_temperature()
    df_sumup = df_sumup.loc[df_sumup.depth==10, :]
    
    df_all_temps = pd.DataFrame()
    for stid in station_list:
        run_name = stid + '_100_layers_3H'
        if not os.path.isfile(output_path+'/'+ run_name+'/constants.csv'):
            continue
        tmp =pd.read_csv(output_path+'/'+ run_name+'/constants.csv', dtype={'key':str})
        tmp['value_num'] = pd.to_numeric(tmp.value, errors='coerce')
        tmp.loc[tmp.value_num.notnull(),'value'] = tmp.loc[tmp.value_num.notnull(),'value_num']
        tmp = tmp.set_index('key')[['value']]
        c = Struct(**tmp.to_dict()['value'] )
        print(c.station)
        filename = c.output_path + run_name + "/" + c.station + "_T_ice.nc"
        if not os.path.isfile(filename): continue
        df_out = (xr.open_dataset(filename).to_dataframe().unstack('level'))
        df_out.columns = df_out.columns.map('{0[0]}_{0[1]}'.format)
        df_out['t_i_10m'] = lplt.interpolate_temperature_fast(
            df_out.index, df_out[[v for v in df_out.columns if 'depth' in v]].values, 
            df_out[[v for v in df_out.columns if 'T_ice' in v]].values-273.15, 
        )
        
        df_sumup_selec, df_meta_selec = lplt.select_sumup(df_sumup, df_meta, c)
        if len(df_sumup_selec)==0: continue
        df_sumup_selec['duration'] = df_sumup_selec.duration.fillna(3).values
        df_sumup_selec['end_date'] =  df_sumup_selec.timestamp.values
        df_sumup_selec['start_date'] = df_sumup_selec.timestamp.values - pd.to_timedelta(
            df_sumup_selec.duration, unit='h')
        for i in df_sumup_selec.index:
            df_sumup_selec.loc[i, 'temperature_mod'] = df_out.loc[
                df_sumup_selec.loc[i,'start_date']:df_sumup_selec.loc[i,'end_date'],
                't_i_10m'].mean()
        df_sumup_selec = df_sumup_selec.set_index('timestamp')
         
        ind = df_sumup_selec.index.intersection(df_out.index)
        df_sumup_selec = df_sumup_selec.loc[ind]
        
        df_sumup_selec['station'] = c.station
        df_all_temps = pd.concat((df_all_temps, df_sumup_selec.reset_index()), ignore_index=True)
    df_all_temps.to_csv('./side analysis/T10m_evaluation.csv')    
    
df_all_temps = pd.read_csv('./side analysis/T10m_evaluation.csv')    

fig, ax = plt.subplots(1,1,figsize=(7,10))
plt.subplots_adjust(bottom=0.5,left=0.2,right=0.8,top=0.98)
for stid in ds_aws.stid.data[1:]:
    print(stid)
    df_sumup_selec = df_all_temps.loc[df_all_temps.station==stid,:].set_index('index')
    ax.plot(df_sumup_selec.temperature, df_sumup_selec.temperature_mod, marker='.',
            ls='None',label=stid)
lim1=[-35, 0.15]
ax.plot(lim1,lim1, color='k')
ax.set_xlim(lim1)
ax.set_ylim(lim1)
ax.grid(lim1)
ax.legend(loc='lower center',ncol=4, bbox_to_anchor=(0.5,-1.05))
ax.set_xlabel('Observed T10m (m w.e.)')
ax.set_ylabel('Modelled T10m (m w.e.)')
fig.savefig('side analysis/T10m_evaluation_SUMup2024.png', dpi=300)

# %% Grain size evaluation
ds_gs = pd.read_csv(r'C:\Users\bav\GitHub\SUMup\grain-size\output\grain_size_ssa_compilation.csv')
output_path = './output/'
run_name = 'Summit_100_layers_3H'
plt.close('all')
ds = xr.open_dataset(output_path+run_name+'/Summit_dgrain.nc')

for tag in ['Lomonaco','Tedesco']:
    fig, ax = plt.subplots(1,2)
    plt.subplots_adjust(top=0.8)
    ds_selec = ds_gs.loc[ds_gs.reference_short.str.startswith(tag),:]
    tmp = ds.sel(time=ds_selec.timestamp.iloc[0],method='nearest')
    ax[0].plot(tmp.dgrain, -tmp.depth,label='model')
    ax[1].plot(tmp.dgrain, -tmp.depth, label='model')
    ax[0].plot(ds_selec.grain_diameter_mm, -ds_selec.start_depth_m,  marker='o',ls='None',label=ds_selec.reference_short.iloc[0])
    ax[1].plot(ds_selec.grain_diameter_mm, -ds_selec.start_depth_m,  marker='o',ls='None',label=ds_selec.reference_short.iloc[0])
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")
    for i in [0,1]:
        ax[i].set_ylabel('Depth (m)')
        ax[i].set_xlabel('Grain diameter (mm)')
        ax[i].grid()
    ax[1].set_ylim(-4,0)
    ax[1].set_xlim(0,3)
    ax[0].set_xlim(0,5)
    ax[0].legend(title='Summit '+ds_selec.timestamp.iloc[0], loc='lower right',
                 bbox_to_anchor=(1.5,1))
    fig.savefig('side analysis/grain_size_evaluation_'+tag+'.png', dpi=300)

run_name = 'EastGRIP_100_layers_3H'
ds = xr.open_dataset(output_path+run_name+'/EastGRIP_dgrain.nc')

for tag in ['Montagnat']:
    fig, ax = plt.subplots(1,1)
    ax=[ax]
    plt.subplots_adjust(top=0.8)
    ds_selec = ds_gs.loc[ds_gs.reference_short.str.startswith(tag),:]
    tmp = ds.sel(time=ds_selec.timestamp.iloc[0],method='nearest')
    ax[0].plot(tmp.dgrain, -tmp.depth, label='model')
    for name in ds_selec.name.unique():
        ax[0].plot(ds_selec.loc[ds_selec.name==name, 'grain_diameter_mm'],
                   -ds_selec.loc[ds_selec.name==name, 'start_depth_m'], 
                   marker='o',ls='None',
                   label=name)

    ax[0].set_ylabel('Depth (m)')
    ax[0].set_xlabel('Grain diameter (mm)')
    ax[0].grid()

    ax[0].set_xlim(0,2.5)
    ax[0].set_ylim(-3.1,0.1)
    ax[0].legend(title='EastGRIP '+ds_selec.timestamp.iloc[0], loc='lower right',
                 bbox_to_anchor=(0.7,1))
    fig.savefig('side analysis/grain_size_evaluation_'+tag+'.png', dpi=300)

# %% LWC at several sites
station_list = ['H2', 'KAN_U', 'FA-13','Summit']
from lib.plot import load_sumup_temperature, select_sumup, plot_var
df_sumup, df_meta = load_sumup_temperature()
# %% 
for stid in station_list:
    # %%
    stid='FA-13'
    output_path = './output/'
    run_name = stid+'_100_layers_3H'
    
    print(run_name)
    try: 
        tmp =pd.read_csv(output_path+'/'+ run_name+'/constants.csv', dtype={'key':str})
    except:
        # continue
        pass
    tmp['value_num'] = pd.to_numeric(tmp.value, errors='coerce')
    tmp.loc[tmp.value_num.notnull(),'value'] = tmp.loc[tmp.value_num.notnull(),'value_num']
    tmp = tmp.set_index('key')[['value']]
    c = Struct(**tmp.to_dict()['value'] )
    c.RunName=run_name
        
    df_sumup_selec, df_meta_selec = select_sumup(df_sumup, df_meta, c)

    # infiltration evaluation
    plot_var(c.station, c.output_path, c.RunName, 'slwc', zero_surf=True, 
                 df_sumup=df_sumup_selec, ylim=[20], tag='_SUMup2024_slwc')

#%% SMB at some sites

station_list = ['H2', 'KAN_U', 'FA-13','Summit']

plt.figure()
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
    
    df_in, c = io.load_surface_input_data(c, resample=False)
    if output_path != c.output_path:
        print('Warning: Output has been moved from',c.output_path,'to',output_path)
        c.output_path = output_path
        
    #  loading surface variables
    try:
        df_out = xr.open_dataset(c.output_path+run_name+'/'+c.station+'_surface.nc').to_dataframe()
        df_in = df_in.loc[df_out.index[0]:df_out.index[-1],:]
    except Exception as e:
        print(e)
    df_melt_a = df_out.melt_mweq.resample('Y').sum()
    df_melt_a.plot(label=stid+' mean: %0.2f m a-1'%df_melt_a.mean().item())
plt.ylabel('Annual surface melt (m w.e.)')
plt.grid()
plt.legend()