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
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import matplotlib
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

def plot_var(site, output_path, run_name, var_name, ylim=[], zero_surf=True,
             df_sumup=[], tag='', year=None):
    print('plotting',var_name, 'from',run_name)
    filename = output_path+"/" + run_name + "/" + site + "_" + var_name + ".nc"
    ds = xr.open_dataset(filename).transpose()
    
    if year:
        if len(year) == 2:
            ds = ds.sel(time=slice(str(year[0]), str(year[1])))
            tag='_'+str(year[0])+'_'+str(year[1])
        else:
            ds = ds.sel(time=str(year))
            tag='_'+str(year)
    # if len(df_sumup)>0:
    #     ds = ds.sel(time=slice(df_sumup.timestamp.min(),
    #                       df_sumup.timestamp.max())
    #                 ).where((ds.depth<df_sumup.depth.max()+15).all(dim='time'), drop=True)
    
    if zero_surf:
        ds['surface_height'] = 0 * ds.depth.isel(level=-1)
    else:
        ds['surface_height'] = -(ds.depth.isel(level=-1)
                                -ds.depth.isel(level=-1).isel(time=0)
                                -(ds.depth.isel(level=-1)
                                  .diff(dim='time')
                                  .where(ds.depth.isel(level=-1)
                                         .diff(dim='time')>6,0)
                                  .cumsum()))
        ds['surface_height'].values[0] = 0
        ds['depth'] = ds.depth + ds.surface_height

    ds = ds.resample(time='6H').nearest()

    if var_name == "slwc":
        # change unit to mm / m3
        ds[var_name] = ds[var_name] * 1000 / ds.depth
    if var_name == "T_ice":
        # change unit to mm / m3
        ds[var_name] = ds[var_name] -273.15
    
    # default plot infos
    label = var_name
    cmap = 'magma'
    vmin = np.percentile(ds[var_name], 5)
    vmax = np.percentile(ds[var_name], 95)
            
    # updating for pre-set values
    plot_info = pd.read_csv('lib/plot_info.csv', skipinitialspace=True)
    if var_name in plot_info.variable_name.to_list():
        plot_info = plot_info.set_index('variable_name')
        label = plot_info.loc[var_name, 'label']
        cmap = plot_info.loc[var_name, 'cmap']
        if ~np.isnan(plot_info.loc[var_name].vmin):
            vmin = plot_info.loc[var_name, 'vmin']
            vmax = plot_info.loc[var_name, 'vmax']
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    plt.subplots_adjust(left=0.07, right=0.99, top=0.95, bottom=0.1, hspace=0.2)
    fig.suptitle(site)
    
    im = ax.pcolormesh(
                   ds.time.expand_dims(dim={"level": ds.level.shape[0]+1}).transpose(), 
                   np.hstack([ds.surface_height.values.reshape([-1,1]),
                                  ds.depth.values]),
                   ds[var_name].isel(time=slice(1,None)),
                   shading='flat',
                   cmap = cmap, vmin=vmin, vmax=vmax,
                   zorder=0
                   )
    
    if not zero_surf:
        ax.plot(ds.time, ds.surface_height, linewidth=2, color="k")
        
    plt.colorbar(im, label=label, ax=ax)
    ax.invert_yaxis()
    if ylim:
        if len(ylim)==1: ax.set_ylim(ylim, ax.get_ylim()[1])
        if len(ylim)==2: ax.set_ylim(np.max(ylim), np.min(ylim))
    ax.set_ylabel("Depth (m)")
    
    # adding SUMup observations as colored scatter
    if len(df_sumup)>0:
        if var_name == 'T_ice':
            depth_var = 'depth'
            sumup_var = 'temperature'
        if var_name == 'density_bulk':
            depth_var = 'midpoint'
            sumup_var = 'density'
            
        plt.plot(df_sumup.timestamp,
            df_sumup[depth_var],
            marker='o',ls='None', markersize=6, color='lightgray',zorder=1)
        plt.scatter(df_sumup.timestamp,
                    df_sumup[depth_var],
                    12, df_sumup[sumup_var],
                    vmin=vmin, vmax=vmax,
                    cmap=cmap, zorder=2)
        ax.set_ylim(df_sumup[depth_var].max()+2, 0)
        ax.set_xlim(df_sumup.timestamp.min()-pd.Timedelta('100D'),
                    df_sumup.timestamp.max()+pd.Timedelta('100D'))
    
    fig.savefig(output_path+"/" + run_name + "/" + site + "_" + var_name +tag+ ".png")
    return fig, ax

def plot_var_start_end(c, var_name='T_ice', ylim=[], to_file=False):
    site = c.station
    output_path = c.output_path
    run_name = c.RunName
    print('plotting',var_name, 'from',run_name)
    filename = output_path+"/" + run_name + "/" + site + "_" + var_name + ".nc"
    ds = xr.open_dataset(filename).transpose()
    ds = ds.resample(time='6H').nearest()
    
    if var_name == "slwc":
        # change unit to mm / m3
        ds[var_name] = ds[var_name] * 1000 / ds.depth
    if var_name == "T_ice":
        # change unit to mm / m3
        ds[var_name] = ds[var_name] -273.15
        
        T10m = pd.read_csv('../../Data/Firn temperature/output/10m_temperature_dataset_monthly.csv')
        ds_T10m = xr.open_dataset('../../Data/Firn temperature/output/T10m_prediction.nc')

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    # plt.subplots_adjust(left=0.07, right=0.99, top=0.95, bottom=0.1, hspace=0.2)
    plt.plot( ds[var_name].isel(time=0), 
             ds.depth.isel(time=0),
             marker='o',
             color='tab:blue',
             label=ds[var_name].time.isel(time=0).dt.strftime("%Y-%m-%d").item(),
             )

    plt.plot(ds[var_name].sel(time='2023-09-01 00:00:00'), 
             ds.depth.sel(time='2023-09-01 00:00:00'),
             marker='o',
             color='tab:red',
             label='2023-09-01 00:00:00',
             )
    if var_name == "T_ice":
        T10m.loc[T10m.site==c.station, :].plot(x='temperatureObserved',
                                               y='depthOfTemperatureObservation',
                                               ax=plt.gca(),
                                               marker='o', ls='None')
        plt.axvline(float(c.Tdeep)-273.15, ls='--')
        plt.axvline(ds.isel(level=0).median().T_ice, ls='-.')
        plt.axvline(ds.isel(level=1).median().T_ice, ls=':')
        plt.axvline(ds.isel(level=10).median().T_ice, ls=':')
        plt.axvline(ds.isel(level=20).median().T_ice, ls=':')
        plt.axvline(    ds_T10m.sel(latitude = T10m.loc[T10m.site==c.station, 'latitude'].mean(),
                        longitude = T10m.loc[T10m.site==c.station, 'longitude'].mean(),
                        method='nearest').sel(time=slice('1980','1990')).T10m.mean().item(), c='tab:red')
        plt.axvline(    ds_T10m.sel(latitude = T10m.loc[T10m.site==c.station, 'latitude'].mean(),
                        longitude = T10m.loc[T10m.site==c.station, 'longitude'].mean(),
                        method='nearest').sel(time=slice('2000','2020')).T10m.mean().item(), 
                    ls='--', c='tab:red')

    plt.legend()
    plt.title(site)
    plt.gca().invert_yaxis()
    if ylim:
        if len(ylim)==1: ax.set_ylim(ylim, ax.get_ylim()[1])
        if len(ylim)==2: ax.set_ylim(np.max(ylim), np.min(ylim))
    ax.set_ylabel("Depth (m)")
    if to_file:
        (ds[[var_name]+['depth']]
         .sel(time='2023-09-01 00:00:00')
         .to_dataframe()[[var_name]+['depth']]
         .set_index('depth')
         .to_csv('input/initial state/spin up/'+c.station+'_initial_'+var_name+'.csv'))
        
    fig.savefig(output_path+"/" + run_name + "/" + site + "_" + var_name + "_start_end.png")
    return fig, ax

def plot_movie(site, output_path,  run_name, var_name, ylim=[]):
    print('plotting',var_name, 'from',run_name)
    filename = output_path+"/" + run_name + "/" + site + "_" + var_name + ".nc"
    ds = xr.open_dataset(filename).transpose()
    ds = ds.resample(time='6H').nearest()
    
    if var_name == "slwc":
        # change unit to mm / m3
        ds[var_name] = ds[var_name] * 1000 / ds.depth
    if var_name == "T_ice":
        # change unit to mm / m3
        ds[var_name] = ds[var_name] -273.15
        
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from tqdm import tqdm
    
    # Assuming ds is your dataset and var_name is the variable name
    # Adjust the parameters accordingly
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    
    def update(frame):
        ax.clear()
        ax.plot(ds[var_name].isel(time=frame), 
                ds.depth.isel(time=frame),
                marker='o',
                color='tab:blue',
                label=ds[var_name].time.isel(time=frame),
                )
        ax.set_title(f'Timestamp: {ds[var_name].time.isel(time=frame).values}')
        ax.set_xlim(ds[var_name].min(), ds[var_name].max())
        ax.set_ylim(80,0)
        
    animation = FuncAnimation(fig, update, 
                              frames=range(0,len(ds.time),24*7),
                              interval=100, repeat=False)
    
    # Save the animation as an MP4 file
    animation.save(output_path + '/'+run_name+'/'+var_name+'.gif', fps=30, writer='pillow')
    plt.show()

def track_horizon(time, H_surf, depth_act, compaction, date_start, depth_start, step=1):
    ind_start = (np.abs(time - date_start)).argmin()

    length_out = len(time)
    depth_hor = np.empty(length_out) * np.nan
    depth_hor[ind_start] = depth_start

    for i in range(ind_start + step, len(time), step):
        depth_hor[i] = max(0, depth_hor[i - step] + (H_surf[i] - H_surf[i - step]))
        depth_mod = depth_act[:, i]
        comp_mod = compaction[:, i] * step

        # compaction in all layers below the horizon
        ind_next = int(
            interp1d(
                np.insert(depth_mod, 0, 0), np.arange(len(depth_mod) + 1), kind="next"
            )(depth_hor[i])
        )
        ind_prev = int(
            interp1d(
                np.insert(depth_mod, 0, 0), np.arange(len(depth_mod) + 1), 
                kind="previous", fill_value=0
            )(depth_hor[i])
        )
        comp_tot = np.sum(comp_mod[ind_next:]) 

        # plus compaction within the layer where the horizon is
        comp = (
            (depth_mod[ind_next] - depth_hor[i])
            / (depth_mod[ind_next] - depth_mod[ind_next - 1])
            * comp_mod[ind_next - 1]
        )

        # comp_tot = comp_tot + comp

        depth_hor[i] = depth_hor[i] + comp_tot
        
    # interpolating between the steps
    if np.sum(np.isnan(depth_hor)) > 0:
        depth_hor[np.isnan(depth_hor)] = interp1d(
            np.argwhere(~np.isnan(depth_hor)).transpose()[0],
            depth_hor[~np.isnan(depth_hor)],
            kind="linear",
            fill_value="extrapolate",
        )(np.argwhere(np.isnan(depth_hor)).transpose()[0])
        depth_hor[:ind_start] = np.nan
    return depth_hor

def evaluate_compaction(c):
    site = c.station
    output_path = c.output_path
    run_name = c.RunName
    filename = output_path+"/" + run_name + "/" + site + "_compaction.nc"
    ds = nc.Dataset(filename)
    compaction = ds["compaction"][:]
    time_org = np.asarray(ds["time"][:])
    time = np.datetime64("1900-01-01T00") + np.round((time_org) * 24 * 3600) * np.timedelta64(1, "s")
    depth_act = np.asarray(ds["depth"][:])
    H_surf = depth_act[-1, :] - depth_act[-1, 0]
    df_comp_info = pd.read_csv(
        "side analysis/Firn viscosity/Compaction_Instrument_Metadata.csv"
    ).set_index("sitename")
    
    if site == 'DY2': site='DYE-2'
    if site == 'NSE': site='NASA-SE'
    if site == 'SDL': site='Saddle'
    if site == 'SUM': site='Summit'
    if site == 'EastGRIP': site='EGP'
    
    if site not in df_comp_info.index.unique():
        return None
    
    df_comp_info = (df_comp_info.loc[site, ["instrument_ID",
                                        "installation_daynumber_YYYYMMDD",
                                        "borehole_top_from_surface_m",
                                        "borehole_bottom_from_surface_m"] ]
                                .reset_index(drop=True)
                                .set_index("instrument_ID") )

    df_comp = pd.read_csv("side analysis/Firn viscosity/borehole_shortening_m.csv")
    df_comp.date = pd.to_datetime(df_comp.date)
    df_comp = df_comp.set_index(["instrument_id", "date"])

    fig1, ax = plt.subplots(1, 1)  
    plot_var(c.station, output_path, run_name, 'density_bulk')
    ax.plot(time, -H_surf, label="Surface")

    fig2, ax2 = plt.subplots(len(df_comp_info.index), figsize=(10, 25), sharex=True)
    fig2.suptitle(site)
    fig2.subplots_adjust(left=0.1, right=0.99, top=0.9, hspace=0.3)

    cmap = matplotlib.cm.get_cmap("Spectral")

    for i, ID in enumerate(df_comp_info.index):
        print(site, ID)
        if ID not in df_comp.index.get_level_values(0).unique():
            print("No data")
            ax2[i].set_title("Instrument " + str(ID) + ": no data")
            continue
        date_start = pd.to_datetime(
            str(df_comp_info.loc[ID, "installation_daynumber_YYYYMMDD"])
        ).to_datetime64()
        depth_top = df_comp_info.loc[ID, "borehole_top_from_surface_m"]
        depth_bot = -df_comp_info.loc[ID, "borehole_bottom_from_surface_m"]

        depth_1 = track_horizon(time, H_surf, depth_act, compaction,  date_start, depth_top, step=12)
        depth_2 = track_horizon(time, H_surf, depth_act, compaction,  date_start, depth_bot, step=12)

        ax.plot(time, depth_1 - H_surf, color=cmap(i / len(df_comp_info.index)), label="_no_legend_")
        ax.plot(time, depth_2 - H_surf, color=cmap(i / len(df_comp_info.index)), label="Instrument " + str(ID))

        # ax2[i] = plt.subplot(2,1,1)
        # ax2[i].plot(time, (depth_2-depth_1) - (depth_2[0]-depth_1[0]))
        # df_comp.loc[ID,'borehole_shortening_m'].plot(ax=ax1)
        # ax2[i].set_title(site + 'Instrument '+str(ID))
        # ax1 = plt.subplot(2,1,2)
        df_comp.loc[ID, "borehole_shortening_m"].diff().plot(
            ax=ax2[i], label="Observation"
        )
        ax2[i].plot(time[:-1], np.diff(depth_2 - depth_1) * 24, label="Simulated")
        ax2[i].set_ylim(-0.004, 0.001)
        tmp = df_comp.loc[
            np.isin(df_comp.index.get_level_values(0), df_comp_info.index), :
        ].index.get_level_values(1)
        ax2[i].set_xlim(tmp.min(), tmp.max())
        ax2[i].set_title("Instrument " + str(ID))
    ax.set_title(site)
    ax.legend()
    ax.grid()
    ax.set_ylabel("Depth (m)")
    ax.set_ylim(np.nanmin(depth_2 - H_surf), -np.nanmax(H_surf))
    ax2[i].legend()
    fig2.text(0.03, 0.5, "Compaction rate (m d$^{-1}$)",
        ha="center", va="center", rotation="vertical")
    fig1.savefig(output_path+"/" + run_name + "/" + site + "_compaction_1.png")
    fig2.savefig(output_path+"/" + run_name + "/" + site + "_compaction_2.png")


def find_summer_surface_depths(c):
    site = c.station
    output_path = c.output_path
    run_name = c.RunName
    filename = output_path+"/" + run_name + "/" + site + "_compaction.nc"
    ds = nc.Dataset(filename)
    compaction = ds["compaction"][:]
    time_org = np.asarray(ds["time"][:])
    time = np.datetime64("1900-01-01T00") + np.round((time_org) * 24 * 3600) * np.timedelta64(1, "s")
    depth_act = np.asarray(ds["depth"][:])
    H_surf = depth_act[-1, :] - depth_act[-1, 0]
    df_ssd = pd.DataFrame(index=pd.to_datetime(time))
    df_ssd['H_surf'] = H_surf
    # list_years = [2015, 2016, 2017, 2018, 2019, 2020, 2021] 
    list_years = df_ssd.index.year.unique()
    for i, yr in enumerate(list_years):
        print(yr)
        date_start = np.datetime64(str(yr)+'-09-01')
        df_ssd[str(yr)] = track_horizon(time, H_surf, depth_act, compaction,  date_start, 0, step=48)
        
    cmap = matplotlib.cm.get_cmap("Spectral")
    fig1, ax = plt.subplots(2, 1,figsize=(8,8))  
    
    ax[0].plot(time, -H_surf, label="Surface")
    for i, yr in  enumerate(list_years):
        ax[0].plot(time, df_ssd[str(yr)].values - df_ssd.H_surf, 
                color=cmap(i / len(df_ssd.index.year.unique()[:-1])), 
                label="_no_legend_")

    ax[0].set_title(site)
    ax[0].grid()
    ax[0].set_ylabel("Height (m)")
    ax[0].set_ylim(-np.nanmin(H_surf), -np.nanmax(H_surf))
    
    for i, yr in  enumerate(list_years):
        ax[1].plot(time, -df_ssd[str(yr)].values, 
                color=cmap(i / len(df_ssd.index.year.unique()[:-1])), 
                label="_no_legend_")
    ax[1].grid()
    ax[1].set_ylabel("Depth (m)")

    fig1.savefig(output_path+"/" + run_name + "/" + site + "_summer_surface_depth_step_48.png",dpi=240)
    df_ssd.resample('D').mean().to_csv(output_path+"/" + run_name + "/" + site + "_summer_surface_depth.csv")


def plot_summary(df, c, filetag="summary", var_list=None):
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
                c.output_path + "/" + c.RunName + "/" + c.station + "_summary_" + str(count_fig),
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
        c.output_path + "/" + c.RunName + "/" + c.station + "_summary_" + str(count_fig),
        bbox_inches="tight",
    )
    

from scipy.spatial import distance
from math import sin, cos, sqrt, atan2, radians

def get_distance(point1, point2):
    R = 6370
    lat1 = radians(point1[0])  #insert value
    lon1 = radians(point1[1])
    lat2 = radians(point2[0])
    lon2 = radians(point2[1])

    dlon = lon2 - lon1
    dlat = lat2- lat1

    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1-a))
    distance = R * c
    return distance

def evaluate_temperature_sumup(df_out, c):
    # Evaluating temperature with SUMup 2024
    df_sumup = xr.open_dataset(
        'C:/Users/bav/GitHub/SUMup/SUMup-2024/SUMup 2024 beta/SUMup_2024_temperature_greenland.nc',
        group='DATA').to_dataframe()
    ds_meta = xr.open_dataset(
        'C:/Users/bav/GitHub/SUMup/SUMup-2024/SUMup 2024 beta/SUMup_2024_temperature_greenland.nc',
        group='METADATA')
    
    df_sumup.method_key = df_sumup.method_key.replace(np.nan,-9999)
    # df_sumup['method'] = ds_meta.method.sel(method_key = df_sumup.method_key.values).astype(str)
    df_sumup['name'] = ds_meta.name.sel(name_key = df_sumup.name_key.values).astype(str)
    df_sumup['reference'] = (ds_meta.reference
                             .drop_duplicates(dim='reference_key')
                             .sel(reference_key=df_sumup.reference_key.values)
                             .astype(str))
    df_sumup['reference_short'] = (ds_meta.reference_short
                             .drop_duplicates(dim='reference_key')
                             .sel(reference_key=df_sumup.reference_key.values)
                             .astype(str))
    # df_ref = ds_meta.reference.to_dataframe()
    df_sumup = df_sumup.loc[df_sumup.timestamp>pd.to_datetime('1989')]
    # selecting Greenland metadata measurements
    df_meta = df_sumup.loc[df_sumup.latitude>0, 
                      ['latitude', 'longitude', 'name_key', 'name', 'method_key',
                       'reference_short','reference', 'reference_key']
                      ].drop_duplicates()
    
    query_point = [[c.latitude, c.longitude]]
    all_points = df_meta[['latitude', 'longitude']].values
    df_meta['distance_from_query_point'] = distance.cdist(all_points, query_point, get_distance)
    min_dist = 10 # in km
    df_meta_selec = df_meta.loc[df_meta.distance_from_query_point<min_dist, :]   

    df_sumup = df_sumup.loc[
        df_sumup.latitude.isin(df_meta_selec.latitude)&df_sumup.longitude.isin(df_meta_selec.longitude),:]
    plot_var(c.station, c.output_path, c.RunName, 'T_ice', zero_surf=True, 
                 df_sumup=df_sumup, tag='_SUMup2024')

    fig,ax = plt.subplots(1,1,figsize=(7,7))
    plt.subplots_adjust(bottom=0.4)
    cmap = matplotlib.cm.get_cmap('tab10')

    for count, ref in enumerate(df_meta_selec.reference_short.unique()):
        label = ref
        for n in df_meta_selec.loc[df_meta_selec.reference_short==ref, 'name_key'].drop_duplicates().values:
            df_subset=df_sumup.loc[(df_sumup.name_key == n)&(df_sumup.depth == 10), :]
            if len(df_subset)>0:
                df_subset.plot(ax=ax, x='timestamp', y='temperature',
                              color = cmap(count),
                              marker='o',ls='None',
                              label=label, alpha=0.4, legend=False
                              )

    df_out.t_i_10m.plot(ax=ax,color='tab:red', label='GEUS model')
    ax.set_ylabel('10 m temperature (°C)')
    ax.set_xlabel('')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1))
    plt.title(c.station)
    fig.savefig(c.output_path+c.RunName+'/T10m_evaluation_SUMup2024.png', dpi=120)

def evaluate_density_sumup(c):
    # Evaluating density with SUMup 2024
    df_sumup = xr.open_dataset(
        'C:/Users/bav/GitHub/SUMup/SUMup-2024/SUMup 2024 beta/SUMup_2024_density_greenland.nc',
        group='DATA').to_dataframe()
    ds_meta = xr.open_dataset(
        'C:/Users/bav/GitHub/SUMup/SUMup-2024/SUMup 2024 beta/SUMup_2024_density_greenland.nc',
        group='METADATA')
    
    df_sumup.method_key = df_sumup.method_key.replace(np.nan,-9999)
    # df_sumup['method'] = ds_meta.method.sel(method_key = df_sumup.method_key.values).astype(str)
    df_sumup['profile'] = ds_meta.profile.sel(profile_key = df_sumup.profile_key.values).astype(str)
    df_sumup['reference'] = (ds_meta.reference
                             .drop_duplicates(dim='reference_key')
                             .sel(reference_key=df_sumup.reference_key.values)
                             .astype(str))
    df_sumup['reference_short'] = (ds_meta.reference_short
                             .drop_duplicates(dim='reference_key')
                             .sel(reference_key=df_sumup.reference_key.values)
                             .astype(str))
    # df_ref = ds_meta.reference.to_dataframe()
    df_sumup = df_sumup.loc[df_sumup.timestamp>pd.to_datetime('1989')]
    # selecting Greenland metadata measurements
    df_meta = df_sumup.loc[df_sumup.latitude>0, 
                      ['latitude', 'longitude', 'profile_key', 'profile', 'method_key',
                       'reference_short','reference', 'reference_key']
                      ].drop_duplicates()
    
    query_point = [[c.latitude, c.longitude]]
    all_points = df_meta[['latitude', 'longitude']].values
    df_meta['distance_from_query_point'] = distance.cdist(all_points, query_point, get_distance)
    min_dist = 10 # in km
    df_meta_selec = df_meta.loc[df_meta.distance_from_query_point<min_dist, :]   

    df_sumup = df_sumup.loc[
        df_sumup.latitude.isin(df_meta_selec.latitude)&df_sumup.longitude.isin(df_meta_selec.longitude),:]
    plot_var(c.station, c.output_path, c.RunName, 'density_bulk', zero_surf=True, 
                 df_sumup=df_sumup, tag='_SUMup2024')
    
    # plot each profile
    filename = c.output_path+"/" + c.RunName + "/" + c.station + "_density_bulk.nc"
    ds_mod_dens = xr.open_dataset(filename).transpose()

    profile_list = df_sumup.profile_key.drop_duplicates()
    def new_figure(): 
        fig,ax = plt.subplots(1,6, figsize=(16,7))
        plt.subplots_adjust(left=0.1, right=0.9, top=0.8, wspace=0.2)
        return fig, ax
    fig,ax = new_figure()
    count = 0
    for i, p in enumerate(profile_list):
        df_profile = df_sumup.loc[df_sumup.profile_key == p, :]
        
        if df_meta.loc[df_meta.profile_key == p, 'reference_short'].item() == 'Clerx et al. (2022)':
            df_profile[['start_depth','stop_depth','midpoint']]
        if df_profile[['start_depth','stop_depth','midpoint']].isnull().all().all():
            print('no data in profile', p, 
                  df_meta.loc[df_meta.profile_key == p, 'profile'].item(),
                  df_meta.loc[df_meta.profile_key == p, 'reference_short'].item())
            continue

            
        df_profile.plot(ax=ax[i-count*6], y='start_depth',x='density',
                        drawstyle="steps-pre",
                        label='observation')

        (ds_mod_dens
         .sel(time=df_profile.timestamp.values[0])
         .to_dataframe()
         .plot(ax=ax[i-count*6],y='depth',x='density_bulk',
               drawstyle="steps-pre",
               label='model',
               color='tab:red'))
        if i-count*6 == 0:
            ax[i-count*6].legend(loc='upper left', ncol=2, bbox_to_anchor=(2.5,1.2))
            ax[i-count*6].set_ylabel('Depth (m)')
        else:
            ax[i-count*6].get_legend().remove()
        title =  (pd.to_datetime(df_profile.timestamp.values[0]).strftime('%Y-%m-%d')
                  + '\n' + df_meta.loc[df_meta.profile_key == p, 'profile'].item()
                  + '\n' + df_meta.loc[df_meta.profile_key == p, 'reference_short'].item())
        ax[i-count*6].set_title(title, fontsize=8, fontweight='bold')
        ax[i-count*6].set_xlabel('Density (kg m$^{-3}$)')
        ax[i-count*6].set_ylim(df_profile.start_depth.max()+1, 0)
        ax[i-count*6].set_xlim(100,1000)
        ax[i-count*6].grid()

        if (i-count*6) == 5: 
            fig.savefig(
                c.output_path+c.RunName+'/'+'density_evaluation_SUMup_'+str(count)+'.png', 
                dpi=120)
            count = count +1
            fig,ax = new_figure()
    if (i-count*6) != 5:
        fig.savefig(
            c.output_path+c.RunName+'/'+'density_evaluation_SUMup_'+str(count)+'.png', 
            dpi=120)
        
def evaluate_accumulation_snowfox(df_in, c):
    # SnowFox
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
