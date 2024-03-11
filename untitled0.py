# -*- coding: utf-8 -*-
"""
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
import matplotlib
import matplotlib.dates as mdates
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 
myFmt = mdates.DateFormatter('%Y-%m-%d')


def plot_model(ds_mod, varname = 'lwc', max_depth = 3,cmap='gist_ncar_r', ax = None):
    depth = np.asarray(ds_mod['depth'][:]).T
    if len(depth) == 200:
        depth = np.tile(depth.T, (len(ds_mod.time), 1)).T
    
    var = np.asarray(ds_mod[varname][:]).T 
    var[depth>max_depth]=np.nan

    depth = np.concatenate((depth[:1,:]*0, depth), axis=0)
    # var = np.concatenate((var, var[-1,:]), axis=0)
    # depth = np.concatenate((depth, depth[:,-1:]), axis=1)
    
    time = ds_mod['time'].values
    # time = np.concatenate((time, time[-1:]))
     
    time_grid = np.expand_dims(time,0)
    time_grid = np.repeat(time_grid,depth.shape[0],axis=0)
    
    if ax is None:
        ax = plt.gca()
    cmap = matplotlib.cm.get_cmap(cmap).copy()
    cmap.set_under('white')
    im = ax.pcolor(time_grid,
              depth,
              var,
              shading ='flat',
              cmap = cmap,
              vmin=1e-3)
    
    plt.colorbar(im, label = 'Liquid water content (mm in layer)', ax=ax)

        
    ax.set_ylim(max_depth, 0) #np.min(depth))
    # ax.xaxis.set_major_formatter(myFmt)
    ax.set_ylabel('Depth (m)')
    # ax[count].invert_yaxis()


def plot_profile_evolution(ds_mod, varname = "lwc", max_depth = 3,cmap='viridis', ax=None):
    max_depth = min(max_depth, ds_mod.depth.max())
    z = np.arange(0, max_depth, 0.01)  # regular grid, every 1 cm
    
    def regrid_profile(ds_mod, varname):
        imax = np.argmax(ds_mod.depth.values > z[-1]) + 1
        var = ds_mod[varname].values[0:imax]
        thickness = np.diff(np.insert(ds_mod.depth.values[0:imax], 0, 0))

        return xr.DataArray(np.interp(z,
                                      ds_mod.depth.values[0:imax], 
                                      var,
                                      right=np.nan), dims='z', coords=[-z])
    snowtemp = ds_mod.groupby('time').map(regrid_profile, args=(varname, ))

    time = matplotlib.dates.date2num(ds_mod.time[0]), matplotlib.dates.date2num(ds_mod.time[-1])
    
    if ax is None:
        ax = plt.gca()
    print(snowtemp.values.min(), snowtemp.values.max())
    im = ax.imshow(snowtemp.values.T[::-1, :], origin='lower', extent=(time[0], time[1], z[-1], 0), aspect='auto', cmap=cmap)
    plt.colorbar(im, label = 'Liquid water content (mm m$^{-3}$)', ax=ax)

    ax.xaxis_date()
    # ax.set_title(varname)


# import time
# current_time = time.time()
# plot_model(lwc)
# # plot_profile_evolution(lwc)
# print(current_time-time.time()) 


# ds_mod = xr.open_dataset(f"BaptisteFirn/RetMIP_GEUS_Dye-2__3hourly_columns.nc")



# %% DYE-2
# Loading heilig
df_heilig = pd.read_csv('BaptisteFirn/observations/Heilig_Dye-2_thermistor.csv')

df_fc = pd.read_csv('BaptisteFirn/observations/FirnCover_DYE-2.csv')
temp_label = ['rtd'+str(i) for i in range(0,24)]
p = np.poly1d([ 1.03093649, -0.49950273])
print('Decreasing FirnCover temperatures by on average: %0.2f degC' % np.mean(df_fc[temp_label] .values-p(df_fc[temp_label].values) ))

df_fc[temp_label]  = p(df_fc[temp_label].values) 

# loading upGPR
df_upgpr = pd.read_csv('BaptisteFirn/observations/PercolationDepth_m_Perco2016.asc',header = None)
df_upgpr.columns = ['matlab_time','perc_depth']
from datetime import datetime, timedelta
df_upgpr['time'] = [datetime.fromordinal(int(t)) + timedelta(days=t%1) - timedelta(days = 366) for t in df_upgpr.matlab_time.values]
df_upgpr=df_upgpr.set_index('time')['perc_depth']
# ds_mod = xr.open_dataset(f"BaptisteFirn/RetMIP_GEUS_Dye-2_16_3hourly_columns.nc")
ds_mod = xr.open_dataset(f"BaptisteFirn/DYE-2_slwc_bin_1.nc")
ds_mod['depth'] = ds_mod['Depth']
ds_mod['lwc'] = ds_mod['slwc']*1000

ds_mod['time'] = ds_mod.time-pd.Timedelta(days=1)

# ds_mod['thickness'] = (('time','level'), np.diff(np.insert(ds_mod['depth'].values,0, ds_mod['depth'].values[:,0]*0,axis= 1)))
# ds_mod['lwc'] = ds_mod['lwc']/ds_mod['thickness']

for yr in range(2015,2018):
    date_start = str(yr)+'-05-01'
    date_end = str(yr)+'-10-01'
    ds_tmp = ds_mod.sel(time = slice(date_start,date_end))

    # fig, ax = plt.subplots(1,1,figsize=(10,7))
    # plot_profile_evolution(ds_tmp, ax=ax, max_depth = 2.5)
    
    fig, ax = plt.subplots(1,1,figsize=(10,5))
    plot_model(ds_tmp, ax = ax,max_depth=2.5, cmap='winter_r')

    
    # FirnCover
    for i in range(0,24):
        tmp = df_fc['depth_'+str(i)].copy()
        tmp.loc[df_fc['rtd'+str(i)]<-0.2] = np.nan
        tmp.loc[df_fc['rtd'+str(i)].isnull()] = np.nan
        ax.plot(pd.to_datetime(df_fc.date), df_fc['depth_'+str(i)], color='yellow')    
        ax.plot(pd.to_datetime(df_fc.date), tmp, linestyle='None', marker='o',color='orange')

    ax.plot(np.nan,np.nan, linestyle='None', marker='o',color='yellow',label='FirnCover sensor depth')
    ax.plot(np.nan,np.nan, linestyle='None', marker='o',color='orange',label='FirnCover temp >-0.2degC')
    
    # Heilig
    if yr>=2016:
        for i in range(1,9):
            tmp = df_heilig['depth_'+str(i)].copy()
            tmp.loc[df_heilig['temp_'+str(i)]<-0.2] = np.nan
            tmp.loc[df_heilig['temp_'+str(i)].isnull()] = np.nan
            ax.plot(pd.to_datetime(df_heilig.time), df_heilig['depth_'+str(i)], color='pink')    
            ax.plot(pd.to_datetime(df_heilig.time), tmp, linestyle='None', marker='o',color='r')
        ax.plot(np.nan,np.nan, linestyle='None', marker='o',color='pink',label='Heilig sensor depth')
        ax.plot(np.nan,np.nan, linestyle='None', marker='o',color='r',label='Heilig temp >-0.2degC')
    df_upgpr.plot(color='k',label='upgpr percolation front')

    import matplotlib.patches as mpatches
    pmark = mpatches.Patch(edgecolor='k',facecolor='white', label='No water')
    leg = fig.legend(handles=[pmark], bbox_to_anchor=(0.416, -0.3, 0.5, 0.5), fontsize=12, frameon=False, markerscale=10)
    plt.title('DYE-2')
    plt.legend(loc ="lower left")
    ax.set_xlim(date_start,date_end) 
    ax.set_ylim(2.5,0)
    fig.savefig('plots/JoG_model_'+str(yr))
    
# %% NASA-SE, CP1, Saddle ----> Not enough melt
site_list = ['NASA-SE', 'CP1', 'Saddle']

site = 'NASA-SE'

df_fc = pd.read_csv('BaptisteFirn/observations/FirnCover_'+site+'.csv')
temp_label = ['rtd' + str(i) for i in range(24)]
levels = np.arange(1,25)
p = np.poly1d([ 1.03093649, -0.49950273])
print('Decreasing FirnCover temperatures by on average: %0.2f degC'%np.mean(df_fc[temp_label] .values-p(df_fc[temp_label].values) ))
df_fc[temp_label]  = p(df_fc[temp_label].values) 

ds_mod = xr.open_dataset('C://Data_save/Data JoG 2020/Corrected/'+site+'_0_SiCL_pr0.001_Ck1.00_darcy_wh0.10/'+site+'_slwc_bin_1.nc')

ds_mod['depth'] = ds_mod['Depth']
ds_mod['lwc'] = ds_mod['slwc']*1000

ds_mod['time'] = ds_mod.time-pd.Timedelta(days=1)

# ds_mod['thickness'] = (('time','level'), np.diff(np.insert(ds_mod['depth'].values,0, ds_mod['depth'].values[:,0]*0,axis= 1)))
# ds_mod['lwc'] = ds_mod['lwc']/ds_mod['thickness']
for yr in range(2015,2018):
    date_start = str(yr)+'-05-01'
    date_end = str(yr)+'-10-01'
    ds_tmp = ds_mod.sel(time = slice(date_start,date_end))

    # fig, ax = plt.subplots(1,1,figsize=(10,7))
    # plot_profile_evolution(ds_tmp, ax=ax, max_depth = 2.5)
    
    fig, ax = plt.subplots(1,1,figsize=(10,5))
    plot_model(ds_tmp, ax = ax, max_depth=2.5, cmap='winter_r')

    # FirnCover
    for i in range(0,24):

        ax.plot(pd.to_datetime(df_fc.date), df_fc['depth_'+str(i)], color='yellow')    
        ax.plot(pd.to_datetime(df_fc.loc[df_fc['rtd'+str(i)]>-0.5].date),
                df_fc.loc[df_fc['rtd'+str(i)]>-0.5,'depth_'+str(i)],
                linestyle='None', marker='o',color='orange')
        ax.plot(pd.to_datetime(df_fc.loc[df_fc['rtd'+str(i)].isnull()].date),
                df_fc.loc[df_fc['rtd'+str(i)].isnull(),'depth_'+str(i)],
                linestyle='None', marker='o',color='lightgray')

    ax.plot(np.nan,np.nan, linestyle='None', marker='o',color='yellow',label='FirnCover sensor depth')
    ax.plot(np.nan,np.nan, linestyle='None', marker='o',color='orange',label='FirnCover temp >-0.2degC')
    # ax.plot(np.nan,np.nan, linestyle='None', marker='o',color='lightgray',label='no recording')

    import matplotlib.patches as mpatches
    pmark = mpatches.Patch(edgecolor='k',facecolor='white', label='No water')
    leg = fig.legend(handles=[pmark], bbox_to_anchor=(0.416, -0.35, 0.5, 0.5), fontsize=12, frameon=False, markerscale=10)
    plt.title(site)

    plt.legend(loc ="lower left")
    ax.set_xlim(pd.to_datetime(date_start), pd.to_datetime(date_end)) 
    ax.set_ylim(2.5,0.1)
    fig.savefig('plots/JoG_model_'+site+'_'+str(yr))
    
# %% KAN_U
site = 'KAN_U'

# Loading firn cover
df_fc = pd.read_csv('BaptisteFirn/observations/FirnCover_'+site+'.csv')
temp_label = ['rtd' + str(i) for i in range(0, 24)]
p = np.poly1d([1.03093649, -0.49950273])
print('Decreasing FirnCover temperatures by on average: %0.2f degC' % np.mean(df_fc[temp_label].values - p(df_fc[temp_label].values)))
df_fc[temp_label] = p(df_fc[temp_label].values)

# Loading PROMICE
df_p = pd.read_csv('BaptisteFirn/observations/'+site+'_PROMICE_thermistor.csv')
for i in range(8):
    df_p['rtd'+str(i)] = df_p['rtd'+str(i)].interpolate(limit= 24*30*2).values
# loading Splaz
df_list = []
for k, note in enumerate(["SPLAZ_main", "SPLAZ_2", "SPLAZ_3"]):
    df = pd.read_csv('BaptisteFirn/observations/'+note+'.csv')
    df_list = df_list + [df]


ds_mod = xr.open_dataset('C://Data_save/PROMICE_model_output/'+site+'_0_IWC_CL_100_layers/'+site+'_slwc.nc')

ds_mod['depth'] = ds_mod['Depth']
ds_mod['lwc'] = ds_mod['slwc']*1000

ds_mod['time'] = ds_mod.time-pd.Timedelta(days=1)

time_mod = ds_mod.time.sel(time=slice('2012-05-01',None)).values
ice_slab_depth = ds_mod.sel(time=slice('2012-05-01',None)).depth.isel(level=-1).values
ice_slab_depth = ice_slab_depth - ice_slab_depth[0] + 2
plt.close('all')
# ds_mod['thickness'] = (('time','level'), np.diff(np.insert(ds_mod['depth'].values,0, ds_mod['depth'].values[:,0]*0,axis= 1)))
# ds_mod['lwc'] = ds_mod['lwc']/ds_mod['thickness']
for yr in range(2009,2018):
    if yr < 2015:
        max_depth = 5
    else:
        max_depth = 2.5
    date_start = str(yr)+'-05-01'
    date_end = str(yr)+'-10-01'
    ds_tmp = ds_mod.sel(time = slice(date_start,date_end))

    # fig, ax = plt.subplots(1,1,figsize=(10,7))
    # plot_profile_evolution(ds_tmp, ax=ax, max_depth = 2.5)
    
    fig, ax = plt.subplots(1,1,figsize=(10,5))
    plot_model(ds_tmp, ax = ax, max_depth=max_depth, cmap='winter_r')

    # FirnCover
    if yr >=2015:
        for i in range(24):
            ax.plot(pd.to_datetime(df_fc.date), df_fc['depth_'+str(i)], color='yellow')    
            ax.plot(pd.to_datetime(df_fc.loc[df_fc['rtd'+str(i)]>-0.2].date),
                    df_fc.loc[df_fc['rtd'+str(i)]>-0.2,'depth_'+str(i)],
                    linestyle='None', marker='o',color='orange')
            ax.plot(pd.to_datetime(df_fc.loc[df_fc['rtd'+str(i)].isnull()].date),
                    df_fc.loc[df_fc['rtd'+str(i)].isnull(),'depth_'+str(i)],
                    linestyle='None', marker='o',color='lightgray')
    
        ax.plot(np.nan,np.nan, linestyle='None', marker='o',color='yellow',label='FirnCover sensor depth')
        ax.plot(np.nan,np.nan, linestyle='None', marker='o',color='orange',label='FirnCover temp >-0.2degC')
        # ax.plot(np.nan,np.nan, linestyle='None', marker='o',color='lightgray',label='no recording')
    
    # PROMICE
    # if yr >=2015:
    for i in range(8):
        ax.plot(pd.to_datetime(df_p.date), df_p['depth_'+str(i)], color='lightgreen')    
        ax.plot(pd.to_datetime(df_p.loc[df_p['rtd'+str(i)]>-0.2].date),
                df_p.loc[df_p['rtd'+str(i)]>-0.2,'depth_'+str(i)],
                linestyle='None', marker='o',color='green')
        ax.plot(pd.to_datetime(df_p.loc[df_p['rtd'+str(i)].isnull()].date),
                df_p.loc[df_p['rtd'+str(i)].isnull(),'depth_'+str(i)],
                linestyle='None', marker='o',color='lightgray')

    ax.plot(np.nan,np.nan, linestyle='None', marker='o',color='lightgreen',label='PROMICE sensor depth')
    ax.plot(np.nan,np.nan, linestyle='None', marker='o',color='green',label='PROMICE temp >-0.2degC')
    
    
    # SPLAZ
    if (yr < 2013)&(yr > 2011):
        for i in range(3):
            df_splaz = df_list[i]
            for i in range(len([col for col in df_splaz.columns if 'depth' in col])):
                
                ax.plot(pd.to_datetime(df_splaz.date), 
                        df_splaz['depth_'+str(i)], color='pink')    
                ax.plot(pd.to_datetime(df_splaz.loc[df_splaz['rtd'+str(i)]>-0.2].date),
                        df_splaz.loc[df_splaz['rtd'+str(i)]>-0.2,'depth_'+str(i)],
                        linestyle='None', marker='o',color='red')
                ax.plot(pd.to_datetime(df_splaz.loc[df_splaz['rtd'+str(i)].isnull()].date),
                        df_splaz.loc[df_splaz['rtd'+str(i)].isnull(),'depth_'+str(i)],
                        linestyle='None', marker='o',color='lightgray')
        
        ax.plot(np.nan,np.nan, linestyle='None', marker='o',color='pink',label='SPLAZ sensor depth')
        ax.plot(np.nan,np.nan, linestyle='None', marker='o',color='red',label='SPLAZ temp >-0.2degC')
        # ax.plot(np.nan,np.nan, linestyle='None', marker='o',color='lightgray',label='no recording')
    plt.plot(time_mod, ice_slab_depth,'k',label='top ice slab')
    import matplotlib.patches as mpatches
    pmark = mpatches.Patch(edgecolor='k',facecolor='white', label='No water')
    leg = fig.legend(handles=[pmark], bbox_to_anchor=(0.416, -0.38, 0.5, 0.5), fontsize=12, frameon=False, markerscale=10)
    plt.title(site)
    plt.legend(loc ="lower left")
    ax.set_xlim(pd.to_datetime(date_start), pd.to_datetime(date_end)) 
    ax.set_ylim(max_depth,0.1)
    fig.savefig('plots/RS_model_'+site+'_'+str(yr))

# %%     
    fig, ax = plt.subplots(1,1,figsize=(10,5))

    # if yr >=2015:
    for i in range(8):
        ax.plot(pd.to_datetime(df_p.date), df_p['rtd'+str(i)].interpolate(limit= 24*30*2))
        # ax.plot(pd.to_datetime(df_p.loc[df_p['rtd'+str(i)]>-0.2].date),
        #         df_p.loc[df_p['rtd'+str(i)]>-0.2,'depth_'+str(i)],
        #         linestyle='None', marker='o',color='green')
        # ax.plot(pd.to_datetime(df_p.loc[df_p['rtd'+str(i)].isnull()].date),
        #         df_p.loc[df_p['rtd'+str(i)].isnull(),'depth_'+str(i)],
        #         linestyle='None', marker='o',color='lightgray')