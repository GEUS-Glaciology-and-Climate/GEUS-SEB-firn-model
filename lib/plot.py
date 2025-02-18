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
from scipy.interpolate import interp1d
import matplotlib
import warnings
import matplotlib as mpl

from matplotlib import gridspec
from scipy.stats import linregress

warnings.filterwarnings("ignore", category=DeprecationWarning)

def plot_var(site, output_path, run_name, var_name, ylim=[], zero_surf=True,
             df_sumup=[], tag='', year=None, weq_depth=False):
    try:
        print('plotting',var_name, 'from',run_name)
        if var_name != 'density_bulk':
            filename = output_path+"/" + run_name + "/" + site + "_" + var_name + ".nc"
            ds = xr.open_dataset(filename).transpose()
        else:
            filename = output_path+"/" + run_name + "/" + site + "_snowc.nc"
            snowc = xr.open_dataset(filename).transpose()
            filename = output_path+"/" + run_name + "/" + site + "_snic.nc"
            snic = xr.open_dataset(filename).transpose()
            filename = output_path+"/" + run_name + "/" + site + "_rhofirn.nc"
            rhofirn = xr.open_dataset(filename).transpose()
            ds = snowc[['depth']]
            ds[var_name] = (snowc.snowc + snic.snic) / (snowc.snowc / rhofirn.rhofirn + snic.snic / 900)

        if weq_depth:
            filename = output_path+"/" + run_name + "/" + site + "_snowc.nc"
            snowc = xr.open_dataset(filename).transpose()
            filename = output_path+"/" + run_name + "/" + site + "_snic.nc"
            snic = xr.open_dataset(filename).transpose()
            filename = output_path+"/" + run_name + "/" + site + "_slwc.nc"
            slwc = xr.open_dataset(filename).transpose()
            ds['depth'] = (snowc.snowc + snic.snic ).cumsum(axis=1)

        if year:
            if len(year) == 2:
                ds = ds.sel(time=slice(str(year[0]), str(year[1])))
                tag='_'+str(year[0])+'_'+str(year[1])
            else:
                ds = ds.sel(time=str(year[0]))
                tag='_'+str(year)
        # if len(df_sumup)>0:
        #     ds = ds.sel(time=slice(df_sumup.timestamp.min(),
        #                       df_sumup.timestamp.max())
        #                 ).where((ds.depth<df_sumup.depth.max()+15).all(dim='time'), drop=True)

        if zero_surf:
            ds['surface_height'] = 0 * ds.depth.isel(level=-1)
        else:
            ds['surface_height'] = -(ds.depth.isel(level=-1)-ds.depth.isel(level=-1).isel(time=0)
                                    -(ds.depth.isel(level=-1).diff(dim='time').where(
                                        np.abs(ds.depth.isel(level=-1).diff(dim='time'))>3,0)
                                      .cumsum()))
            ds['surface_height'].values[0] = 0
            ds['depth'] = ds.depth + ds.surface_height

        # ds = ds.resample(time='6H').nearest()

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

        if var_name == 'slwc':
            cmap = mpl.cm.get_cmap("winter").copy()
            cmap.set_under('white')

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

        # adding SUMup observations as colored scatter
        if len(df_sumup)>0:
            depth_var = 'depth'
            sumup_var = 'temperature'
            cmap_scatter = cmap
            vmin_scatter = vmin
            vmax_scatter = vmax
            if var_name == 'density_bulk':
                depth_var = 'midpoint'
                sumup_var = 'density'

            if var_name == 'slwc':
                h1 = plt.plot(df_sumup.timestamp,
                    df_sumup[depth_var],
                    marker='o',ls='None', markersize=3, color='lightgray',
                    zorder=1, label='thermistor <-0.2 °C')
                h2 = plt.plot(df_sumup.loc[df_sumup[sumup_var]>-0.2, 'timestamp'],
                            df_sumup.loc[df_sumup[sumup_var]>-0.2, depth_var],
                            marker='o',ls='None', markersize=3, alpha=0.5,
                            color='tab:red', zorder=2, label='thermistor >-0.2 °C')
                plt.legend(loc='lower right')
            else:
                plt.plot(df_sumup.timestamp, df_sumup[depth_var],
                    marker='o',ls='None', markersize=6, color='lightgray',zorder=1)
                plt.scatter(df_sumup.timestamp,
                            df_sumup[depth_var],
                            12, df_sumup[sumup_var],
                            vmin=vmin, vmax=vmax,
                            cmap=cmap, zorder=2)

            ax.set_ylim(df_sumup[depth_var].max()+2, 0)
            ax.set_xlim(df_sumup.timestamp.min()-pd.Timedelta('100D'),
                        df_sumup.timestamp.max()+pd.Timedelta('100D'))

        if len(ylim)==1: ax.set_ylim(ylim[0], ax.get_ylim()[1])
        if len(ylim)==2: ax.set_ylim(np.max(ylim), np.min(ylim))
        ax.set_ylabel("Depth (m)")
        if year:
            if len(year) == 2:
                ax.set_xlim(pd.to_datetime(str(year[0])), pd.to_datetime(str(year[1])))
            else:
                ax.set_xlim(pd.to_datetime(str(year[0])), pd.to_datetime(str(year[0]+1)))
        fig.savefig(output_path+"/" + run_name + "/" + site + "_" + var_name +tag+ ".png")
        plt.close(fig)
    except Exception as e:
        print(c.RunName, e)

def plot_var_start_end(c, var_name='T_ice', ylim=[], to_file=False):
    try:
        site = c.station
        output_path = c.output_path
        run_name = c.RunName
        print('plotting',var_name, 'from',run_name)

        if var_name != 'density_bulk':
            filename = output_path+"/" + run_name + "/" + site + "_" + var_name + ".nc"
            ds = xr.open_dataset(filename).transpose()
        else:
            filename = output_path+"/" + run_name + "/" + site + "_snowc.nc"
            snowc = xr.open_dataset(filename).transpose()
            filename = output_path+"/" + run_name + "/" + site + "_snic.nc"
            snic = xr.open_dataset(filename).transpose()
            filename = output_path+"/" + run_name + "/" + site + "_rhofirn.nc"
            rhofirn = xr.open_dataset(filename).transpose()
            ds = snowc[['depth']]
            ds[var_name] = (snowc.snowc + snic.snic) / (snowc.snowc / rhofirn.rhofirn + snic.snic / 900)

        ds = ds.resample(time='6h').nearest()

        if var_name == "slwc":
            # change unit to mm / m3
            ds[var_name] = ds[var_name] * 1000 / ds.depth
        if var_name == "T_ice":
            # change unit to mm / m3
            ds[var_name] = ds[var_name] -273.15

            ds_T10m = xr.open_dataset('input/T10m_prediction.nc')

        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        # plt.subplots_adjust(left=0.07, right=0.99, top=0.95, bottom=0.1, hspace=0.2)
        plt.plot( ds[var_name].isel(time=0),
                 ds.depth.isel(time=0),
                 marker='o',
                 color='tab:blue',
                 label=ds.time.isel(time=0).dt.strftime("%Y-%m-%d").item(),
                 )

        plt.plot(ds[var_name].isel(time=-1),
                 ds.depth.isel(time=-1),
                 marker='o',
                 color='tab:red',
                 label=ds.time.isel(time=-1).dt.strftime("%Y-%m-%d").item(),
                 )
        if var_name == "T_ice":
            plt.axvline(float(c.Tdeep)-273.15, c='k', label='Tdeep')
            plt.axvline(ds.isel(level=0).median().T_ice, ls='-.',label='median T_ice(level=0)')
            plt.axvline(    ds_T10m.sel(latitude = c.latitude,
                                        longitude = c.longitude,
                                        method='nearest').sel(time=slice('1980','1990')
                                        ).T10m.mean().item(),
                        ls='--',label='T10m reconstructed 1980-1990 avg', c='tab:red')
            plt.axvline(    ds_T10m.sel(latitude = c.latitude,
                                        longitude = c.longitude,
                                        method='nearest').sel(time=slice('2000','2020')
                                        ).T10m.mean().item(),
                        ls=':',label='T10m reconstructed 2000-2020 avg', c='tab:orange')

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
        plt.close(fig)
    except Exception as e:
        print(c.RunName, e)
        
        
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
    depth_hor = np.ones(length_out) * np.nan
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
    try:
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

        ID_list = [i for i in df_comp_info.index if i in df_comp.index.get_level_values(0).unique()]
        df_comp = df_comp.reset_index('instrument_id')
        df_comp = df_comp.loc[np.isin(df_comp.instrument_id,ID_list), :]

        fig1, ax = plt.subplots(1, 1)
        # plot_var(c.station, output_path, run_name, 'density_bulk')
        ax.plot(time, -H_surf, label="Surface")

        cmap = matplotlib.cm.get_cmap("Spectral")

        compaction_mod = pd.DataFrame()
        for i, ID in enumerate(ID_list):
            date_start = pd.to_datetime(
                str(df_comp_info.loc[ID, "installation_daynumber_YYYYMMDD"])
            ).to_datetime64()
            depth_top = -df_comp_info.loc[ID, "borehole_top_from_surface_m"]
            depth_bot = -df_comp_info.loc[ID, "borehole_bottom_from_surface_m"]

            depth_1 = track_horizon(time, H_surf, depth_act, compaction,  date_start, depth_top, step=12)
            depth_2 = track_horizon(time, H_surf, depth_act, compaction,  date_start, depth_bot, step=12)
            inst_mod = pd.DataFrame()
            inst_mod['time'] = time
            inst_mod['depth_top'] = depth_1 - H_surf
            inst_mod['depth_bot'] = depth_2 - H_surf
            inst_mod['compaction'] = np.nan
            inst_mod.iloc[:-1,-1] = -np.diff(depth_2 - depth_1) * 24*3600/c.zdtime *1000
            inst_mod['ID'] = ID
            inst_mod = inst_mod.set_index(['time','ID'])
            inst_mod = inst_mod.loc[inst_mod['compaction'].first_valid_index():inst_mod['compaction'].last_valid_index(),:]
            compaction_mod = pd.concat((compaction_mod, inst_mod))

            ax.plot(time, depth_1 - H_surf, color=cmap(i / len(df_comp_info.index)), label="_no_legend_")
            ax.plot(time, depth_2 - H_surf, color=cmap(i / len(df_comp_info.index)), label="Instrument " + str(ID))
        ax.set_title(site)
        ax.legend()
        ax.grid()
        ax.set_ylabel("Depth (m)")
        ax.set_ylim(np.nanmin(depth_2 - H_surf), -np.nanmax(H_surf))
        fig1.savefig(output_path+"/" + run_name + "/" + site + "_compaction_1.png", dpi=240)

        fig2, ax2 = plt.subplots(len(ID_list), figsize=(7, 10), sharex=True)
        fig2.suptitle(site)
        fig2.subplots_adjust(left=0.1, right=0.99, top=0.94, hspace=0.3)


        for i, ID in enumerate(ID_list):

            date_start = pd.to_datetime(
                str(df_comp_info.loc[ID, "installation_daynumber_YYYYMMDD"])
            ).to_datetime64()
            depth_top = -df_comp_info.loc[ID, "borehole_top_from_surface_m"]
            depth_bot = -df_comp_info.loc[ID, "borehole_bottom_from_surface_m"]

            # ax2[i] = plt.subplot(2,1,1)
            # ax2[i].plot(time, (depth_2-depth_1) - (depth_2[0]-depth_1[0]))
            # df_comp.loc[ID,'borehole_shortening_m'].plot(ax=ax1)
            # ax2[i].set_title(site + 'Instrument '+str(ID))
            # ax1 = plt.subplot(2,1,2)
            (-df_comp.loc[df_comp.instrument_id == ID, "borehole_shortening_m"].rolling('14D', center=True).mean().diff()*1000).plot(
                ax=ax2[i], label="observed"
            )
            tmp = compaction_mod.loc[ (slice(None), ID),:].reset_index('ID')
            ax2[i].plot(tmp.index, tmp['compaction'], label="simulated")

            ax2[i].set_xlim(df_comp.borehole_shortening_m.first_valid_index(),
                            df_comp.borehole_shortening_m.last_valid_index())
            ax2[i].set_ylim(0, 2.3)
            ax2[i].grid()
            ax2[i].set_title("Instrument %s, top: %0.1f m, bottom %0.2f m"%(str(ID), depth_top, depth_bot))

            ax2[i].legend()
        fig2.text(0.03, 0.5, "Compaction rate (mm d$^{-1}$)", ha="center", va="center", rotation="vertical")
        fig2.savefig(output_path+"/" + run_name + "/" + site + "_compaction_2.png", dpi=240)

        fig2, ax2 = plt.subplots(len(ID_list), figsize=(7, 10), sharex=True)
        fig2.suptitle(site)
        fig2.subplots_adjust(left=0.1, right=0.99, top=0.9, hspace=0.3)


        for i, ID in enumerate(ID_list):

            date_start = pd.to_datetime(
                str(df_comp_info.loc[ID, "installation_daynumber_YYYYMMDD"])
            ).to_datetime64()
            depth_top = -df_comp_info.loc[ID, "borehole_top_from_surface_m"]
            depth_bot = -df_comp_info.loc[ID, "borehole_bottom_from_surface_m"]


            (df_comp.loc[df_comp.instrument_id == ID, "borehole_shortening_m"]+depth_bot-depth_top).plot(
                ax=ax2[i], label="Observation"
            )
            tmp = compaction_mod.loc[ (slice(None), ID),:].reset_index('ID')
            ax2[i].plot(tmp.index,
                        compaction_mod.loc[ (slice(None), ID),'depth_bot'] - compaction_mod.loc[ (slice(None), ID),'depth_top'],
                        label="Simulated")

            ax2[i].set_xlim(df_comp.borehole_shortening_m.first_valid_index(),
                            df_comp.borehole_shortening_m.last_valid_index())
            ax2[i].grid()
            ax2[i].set_title("Instrument %s, top: %0.1f m, bottom %0.2f m"%(str(ID), depth_top, depth_bot))

            ax2[i].legend()
        fig2.text(0.03, 0.5, "Borehole length (m)", ha="center", va="center", rotation="vertical")
        fig2.savefig(output_path+"/" + run_name + "/" + site + "_compaction_3.png", dpi=240)
        print(c.RunName, 'plotted compaction')

    except Exception as e:
        print(c.RunName, e)
        
from scipy.optimize import curve_fit
import matplotlib
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
def find_summer_surface_depths(c):
    site = c.station
    output_path = c.output_path
    run_name = c.RunName
    filename = output_path+"/" + run_name + "/" + site + "_compaction.nc"
    ds = xr.open_dataset(filename)
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
    # list_years = [2015, 2016, 2017, 2018, 2019, 2020, 2021]
    list_years = df_ssd.index.year.unique()
    # list_years = np.arange(2010,2024)
    for i, yr in enumerate(list_years):
        print(yr)
        date_start = np.datetime64(str(yr)+'-06-01')
        df_ssd[str(yr)] = track_horizon(time, H_surf, depth_act, compaction,  date_start, 0, step=48)
    df_ssd = df_ssd.resample('M').first()
    for m in range(1,7):
        df_ssd.loc['2024-'+str(m).zfill(2)+'-01',:]=np.nan

    di = df_ssd.index
    df_ssd = df_ssd.reset_index().iloc[:,1:]

    def func(x,  c, d):
        return  c * x + d

    guess = (-0.5, 0.5)
    fit_df_ssd = df_ssd.copy()
    col_params = {}
    # Curve fit each column
    for col in fit_df_ssd.columns:
        # Get x & y
        if col == '2023':
            x = fit_df_ssd['2020'].dropna().iloc[:-7].index.astype(float).values
            y = fit_df_ssd['2020'].dropna().iloc[:-7].values
        else:
            x = fit_df_ssd[col].dropna()[:-7].index.astype(float).values
            y = fit_df_ssd[col].dropna()[:-7].values
        # Curve fit column and get curve parameters
        params = curve_fit(func, x, y, guess)
        # Store optimized parameters
        col_params[col] = params[0]

    # Extrapolate each column
    for col in df_ssd.columns:
        # Get the index values for NaNs in the column
        if col != 'H_surf':
            x = df_ssd.loc[
                pd.isnull(df_ssd[col]) & (di>=pd.to_datetime(col+'-06-01')),
                col].index.astype(float).values
        else:
            x = df_ssd[pd.isnull(df_ssd[col])].index.astype(float).values
        # Extrapolate those points with the fitted function
        df_ssd.loc[x, col] = func(x, *col_params[col]) - func(x, *col_params[col])[0] \
             + df_ssd.loc[df_ssd[col].last_valid_index(), col]
    df_ssd.index = di


    fig1, ax = plt.subplots(1, 1,figsize=(12,8))

    # ax[0].plot(time, -H_surf, label="Surface")
    # for i, yr in  enumerate(list_years):
    #     ax[0].plot(time, df_ssd[str(yr)].values - df_ssd.H_surf,
    #             color=cmap(i / len(df_ssd.index.year.unique()[:-1])),
    #             label="_no_legend_")

    # ax[0].set_title(site)
    # ax[0].grid()
    # ax[0].set_ylabel("Height (m)")
    # ax[0].set_ylim(-np.nanmin(H_surf), -np.nanmax(H_surf))

    for i, yr in  enumerate(list_years):
        ax.plot(df_ssd.index, -df_ssd[str(yr)].values,
                color='k',
                label="_no_legend_")
    ax.grid()
    ax.set_ylim(-12,0.1)
    ax.set_xlim(pd.to_datetime('2010-01-01'),pd.to_datetime('2024-06-01'))
    # ax[0].set_ylim(-H_surf[-1]-0.1, -H_surf[-1]+12)
    ax.set_ylabel("Depth (m)")
    ax.xaxis.set_major_locator(matplotlib.dates.YearLocator(base=1))
    ax.xaxis.set_minor_locator(matplotlib.dates.MonthLocator(bymonthday=1))
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.tick_params(axis='both', which='both', right=True, top=True, labelright=True, labeltop=True, labelrotation=0)
    ax.set_title(c.station)
    fig1.savefig(output_path+"/" + run_name + "/" + site + "_summer_surface_depth_step_48.png",dpi=300)
    df_ssd.to_csv(
        output_path+"/" + run_name + "/" + site + "_summer_surface_depth.csv",
        float_format='%.2f')


def plot_summary(df, c, filetag="summary", var_list=None):
    try:
        def new_fig():
            fig, ax = plt.subplots(8, 1, sharex=True, figsize=(15, 10))
            plt.subplots_adjust(
                left=0.1, right=0.9, top=0.97, bottom=0.1, wspace=0.2, hspace=0.05
            )
            return fig, ax

        if not var_list:
            var_list = df.columns
            df.columns = df.columns.astype(str)
            
        # resampling for faster plotting
        df = df[var_list].resample('D').mean()
        
        fig, ax = new_fig()
        count = 0
        count_fig = 0

        for i, var in enumerate(var_list):
            var = str(var)
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
    except Exception as e:
        print(c.RunName, e)

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
    # try:
        df_sumup, df_meta = load_sumup(var='temperature', name_var='name', c=c)

        if len(df_sumup)==0:
            print('no temperature in SUMup for ',c.station)
            return
        # T_ice evaluation
        plot_var(c.station, c.output_path, c.RunName, 'T_ice', zero_surf=True,
                     df_sumup=df_sumup, tag='_SUMup2024')
        # infiltration evaluation
        plot_var(c.station, c.output_path, c.RunName, 'slwc', zero_surf=True,
                     df_sumup=df_sumup, ylim=[10], tag='_SUMup2024_slwc')

        # T10m evaluation
        fig,ax = plt.subplots(1,1,figsize=(7,7))
        plt.subplots_adjust(bottom=0.4)
        cmap = matplotlib.cm.get_cmap('tab10')

        for count, ref in enumerate(df_meta.reference_short.unique()):
            label = ref
            for n in df_meta.loc[df_meta.reference_short==ref, 'name_key'].drop_duplicates().values:
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
        plt.close(fig)
    # except Exception as e:
        # print(c.RunName, e)

# from scipy.interpolate import interp1d
# from tqdm import tqdm  # Import tqdm for the progress bar
# def interpolate_temperature(dates, depth_cor, temp,  depth=10,
#     min_diff_to_depth=2,  kind="linear" ):
#     depth_cor = depth_cor.astype(float)
#     df_interp = pd.DataFrame()
#     df_interp["date"] = dates
#     df_interp["temperatureObserved"] = np.nan

#     # preprocessing temperatures for small gaps
#     tmp = pd.DataFrame(temp)
#     tmp["time"] = dates
#     tmp = tmp.set_index("time")
#     tmp = tmp.resample("H").mean()

#     temp = tmp.loc[dates].values
#     for i in tqdm(range(len(dates)), desc="Interpolating temperatures", unit="date"):
#         x = depth_cor[i, :].astype(float)
#         y = temp[i, :].astype(float)
#         ind_no_nan = ~np.isnan(x + y)
#         x = x[ind_no_nan]
#         y = y[ind_no_nan]
#         x, indices = np.unique(x, return_index=True)
#         y = y[indices]
#         if len(x) < 2 or np.min(np.abs(x - depth)) > min_diff_to_depth:
#             continue
#         f = interp1d(x, y, kind, fill_value="extrapolate")
#         df_interp.iloc[i, 1] = np.min(f(depth), 0)

#     if df_interp.iloc[:5, 1].std() > 0.1:
#         df_interp.iloc[:5, 1] = np.nan
#     return df_interp
def interpolate_temperature_fast(dates, depth_matrix, temp_matrix,  depth=10,
    min_diff_to_depth=2,  kind="linear" ):
    # Choose the depth you want to interpolate to (e.g., 10 meters)
    target_depth = 10
    N = depth_matrix.shape[0]
    M = depth_matrix.shape[1]
    closest_depth_indices = np.abs(depth_matrix - target_depth).argmin(axis=1)
    closest_depths_idx_1 = np.maximum(0, closest_depth_indices - 1)
    closest_depths_idx_2 = np.minimum(M - 1, closest_depth_indices + 1)
    closest_depths = depth_matrix[np.arange(N), closest_depths_idx_1]
    next_closest_depths = depth_matrix[np.arange(N), closest_depths_idx_2]

    temp_at_closest_depths = temp_matrix[np.arange(N), closest_depths_idx_1]
    temp_at_next_closest_depths = temp_matrix[np.arange(N), closest_depths_idx_2]

    weights = (next_closest_depths - target_depth) / (next_closest_depths - closest_depths)
    temp_at_10m = temp_at_closest_depths + weights * (temp_at_next_closest_depths - temp_at_closest_depths)
    return temp_at_10m

        
def evaluate_temperature_scatter(df_out, c, year = None):
    try:
        df_sumup, df_meta = load_sumup(var='temperature', name_var='name',c=c)

        filename = c.output_path+"/" + c.RunName + "/" + c.station + "_T_ice.nc"
        ds = xr.open_dataset(filename).transpose()
        ds['T_ice'] = ds.T_ice -273.15

        if year:
            if len(year) == 2:
                ds = ds.sel(time=slice(str(year[0]), str(year[1])))
                tag='_'+str(year[0])+'_'+str(year[1])
            else:
                ds = ds.sel(time=str(year))
                tag='_'+str(year)
        df_sumup['T_ice_mod'] = np.nan

        time_mat = ds.time.expand_dims(dim={"level": ds.level.shape[0]}).transpose().values
        depth_mat = ds.depth.values

        for i, time_obs in enumerate(df_sumup.timestamp.unique()):
            if i % 1000 == 0: print(np.trunc(i/len(df_sumup.timestamp.unique())*100),'%')
            depth_obs = df_sumup.loc[df_sumup.timestamp == time_obs,'depth'].values
            ds_at_ts = ds.sel(time=time_obs, method='nearest')
            interp_T_ice = (ds_at_ts.set_coords('depth').swap_dims(level='depth')
             .interp(depth=depth_obs, method='linear')
             .T_ice.values)
            df_sumup.loc[df_sumup.timestamp == time_obs,'T_ice_mod'] = interp_T_ice

        # plotting
        fig,ax = plt.subplots(1,1,figsize=(7,7))
        plt.subplots_adjust(bottom=0.4)
        sc = ax.scatter(df_sumup.temperature, df_sumup.T_ice_mod,
                   5, df_sumup.depth, marker='.',
                alpha=0.8,ls='None', cmap='spring')
        min_T = df_sumup[['temperature','T_ice_mod']].min().min()
        max_T = df_sumup[['temperature','T_ice_mod']].max().max()
        ax.plot([min_T, max_T], [min_T, max_T], c='k')
        cb1 = plt.colorbar(sc, label='depth of observation (m)')
        cb1.ax.invert_yaxis()
        RMSE = np.sqrt(np.mean((df_sumup.T_ice_mod - df_sumup.temperature)**2))
        ME = (df_sumup.T_ice_mod - df_sumup.temperature).mean().item()
        N = (df_sumup.T_ice_mod - df_sumup.temperature).count().item()
        ax.grid()

        # Annotate with RMSE and ME
        ax.annotate(f'RMSE: {RMSE:.2f}\nME: {ME:.2f}\nN: {N:.0f}',
                      xy=(0.05, 0.95),
                      xycoords='axes fraction',
                      horizontalalignment='left', verticalalignment='top',
                      fontsize=10, bbox=dict(boxstyle="round,pad=0.3",
                                 alpha=0.5, edgecolor='black', facecolor='white'))

        ax.set_ylabel('Modelled subsurface temperature (°C)')
        ax.set_xlabel('Observed subsurface temperature (°C)')
        plt.title(c.station)
        fig.savefig(c.output_path+c.RunName+'/T10m_evaluation_SUMup2024_scatter.png', dpi=120)
        plt.close(fig)
    except Exception as e:
        print(c.RunName,e)


def evaluate_density_sumup(c):
    # try:
        # Evaluating density with SUMup 2024
        df_sumup, df_meta = load_sumup(var='density',name_var='profile', c=c)

        profile_list = df_sumup.profile_key.drop_duplicates()
        if len(profile_list) == 0:
            print('no density profile in SUMup for',c.station)
            return None


        plot_var(c.station, c.output_path, c.RunName, 'density_bulk', zero_surf=True,
                     df_sumup=df_sumup, tag='_SUMup2024')

        # plot each profile

        filename = c.output_path+"/" + c.RunName + "/" + c.station + "_snowc.nc"
        snowc = xr.open_dataset(filename).transpose()
        filename = c.output_path+"/" + c.RunName + "/" + c.station + "_snic.nc"
        snic = xr.open_dataset(filename).transpose()
        filename = c.output_path+"/" + c.RunName + "/" + c.station + "_rhofirn.nc"
        rhofirn = xr.open_dataset(filename).transpose()
        ds_mod_dens = snowc[['depth']]
        ds_mod_dens['density_bulk'] = (snowc.snowc + snic.snic) / (snowc.snowc / rhofirn.rhofirn + snic.snic / 900)

        def new_figure():
            fig,ax = plt.subplots(1,6, figsize=(16,7))
            plt.subplots_adjust(left=0.1, right=0.9, top=0.8, wspace=0.2)
            return fig, ax
        fig,ax = new_figure()
        count = 0
        for i, p in enumerate(profile_list):
            df_profile = df_sumup.loc[df_sumup.profile_key == p, :].sort_values(by='start_depth')

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
             .sel(time=df_profile.timestamp.values[0], method='nearest')
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
                plt.close(fig)
                count = count +1
                fig,ax = new_figure()
        if (i-count*6) != 5:
            fig.savefig(
                c.output_path+c.RunName+'/'+'density_evaluation_SUMup_'+str(count)+'.png',
                dpi=120)
            plt.close(fig)
    # except Exception as e:
        # print(c.RunName,e)

def load_sumup(var='SMB', name_var='name', c=None):
    with xr.open_dataset(f'../SUMup-2024/SUMup_2024_{var}_greenland.nc',
                         group='DATA') as ds:
        df_sumup = ds.to_dataframe()
        if 'timestamp' in df_sumup.columns:
            df_sumup = df_sumup.loc[
            pd.to_numeric(df_sumup.timestamp.astype(str).str[:4], errors='coerce') > 1989, :]
            # df_sumup['timestamp'] = pd.to_datetime(df_sumup['timestamp'], errors='coerce')
        else:
            # df_sumup['start_date'] = pd.to_datetime(df_sumup['start_date'], errors='coerce')
            # df_sumup['end_date'] = pd.to_datetime(df_sumup['end_date'], errors='coerce')
            df_sumup = df_sumup.loc[df_sumup.end_year > 1989]
            
        df_sumup = df_sumup.loc[df_sumup.latitude>0]
        
        # selecting Greenland metadata measurements
        df_meta = df_sumup.loc[:,
                          ['latitude', 'longitude', name_var+'_key', 'method_key',
                            'reference_key']
                          ].drop_duplicates()
        
        query_point = [[c.latitude, c.longitude]]
        all_points = df_meta[['latitude', 'longitude']].values
        df_meta['distance_from_query_point'] = distance.cdist(all_points, query_point, get_distance)
        min_dist = 5 # in km
        df_meta = df_meta.loc[df_meta.distance_from_query_point<min_dist, :]

        df_sumup = df_sumup.loc[
                    df_sumup.latitude.isin(df_meta.latitude) \
                        & df_sumup.longitude.isin(df_meta.longitude),:]

    print(c.RunName, 'found', len(df_sumup),var, 'measurements in SUMup')
    ds_meta = xr.open_dataset(
        f'../SUMup-2024/SUMup_2024_{var}_greenland.nc',
        group='METADATA')

    # decoding strings as utf-8
    for v in [name_var,'reference','reference_short','method']:
        ds_meta[v] = ds_meta[v].str.decode('utf-8')

    df_sumup.method_key = df_sumup.method_key.replace(np.nan,-9999)
    # df_sumup['method'] = ds_meta.method.sel(method_key = df_sumup.method_key.values).astype(str)
    df_sumup[name_var] = ds_meta[name_var].sel({name_var+'_key': df_sumup[name_var+'_key'].values}).astype(str)
    df_sumup['reference'] = (ds_meta.reference
                             .drop_duplicates(dim='reference_key')
                             .sel(reference_key=df_sumup.reference_key.values)
                             .astype(str))
    df_sumup['reference_short'] = (ds_meta.reference_short
                             .drop_duplicates(dim='reference_key')
                             .sel(reference_key=df_sumup.reference_key.values)
                             .astype(str))   

    df_meta[name_var] = ds_meta[name_var].sel({name_var+'_key': df_meta[name_var+'_key'].values}).astype(str)
    df_meta['reference'] = (ds_meta.reference
                             .drop_duplicates(dim='reference_key')
                             .sel(reference_key=df_meta.reference_key.values)
                             .astype(str))
    df_meta['reference_short'] = (ds_meta.reference_short
                             .drop_duplicates(dim='reference_key')
                             .sel(reference_key=df_meta.reference_key.values)
                             .astype(str))                                
    return df_sumup, df_meta

def plt_smb(ax, df_sumup):
    colors = plt.cm.tab10(np.linspace(0, 1, df_sumup.reference_short.nunique()))
    ref_colors = dict(zip(df_sumup.reference_short.unique(), colors))

    for i in df_sumup.index:
        ref = df_sumup.loc[i, 'reference_short']
        color = ref_colors[ref]
        
        ax[0].plot([df_sumup.loc[i, 'start_date'], df_sumup.loc[i, 'end_date']],
                   df_sumup.loc[i, 'smb'] * np.array([1., 1.]),
                   color=color, marker='x', 
                   label=ref if i == df_sumup[df_sumup.reference_short == ref].index[0] else "")
        ax[0].plot([df_sumup.loc[i, 'start_date'], df_sumup.loc[i, 'end_date']],
                   df_sumup.loc[i, 'smb_mod'] * np.array([1., 1.]),
                   color='k', marker='o',label='__nolegend__')
    ax[0].plot([np.nan, np.nan],
               df_sumup.loc[i, 'smb_mod'] * np.array([1., 1.]),
               color='k', marker='o',label='model')

    ax[0].set_xlabel('Year')
    ax[0].set_ylabel('SMB (m w.e.)')
    ax[0].legend()

    for ref, color in ref_colors.items():
        subset = df_sumup[df_sumup.reference_short == ref]
        ax[1].plot(subset.smb, subset.smb_mod, marker='.', ls='None', color=color, label=ref)

    ax[1].plot([df_sumup.smb_mod.min(), df_sumup.smb_mod.max()],
               [df_sumup.smb_mod.min(), df_sumup.smb_mod.max()],
               color='k')
    ax[1].set_xlabel('Observed SMB (m w.e.)')
    ax[1].set_ylabel('Modelled SMB (m w.e.)')

    return ax



def evaluate_smb_sumup(df_out, c):
    # try:
        df_sumup, df_meta = load_sumup(var='SMB',name_var='name', c=c)

        msk = df_sumup.start_date.isnull()
        df_sumup.loc[msk, 'start_date'] = pd.to_datetime(df_sumup.loc[msk, 'start_year'].astype(int).astype(str)+'-01-01')
        msk = df_sumup.end_date.isnull()
        df_sumup.loc[msk, 'end_date'] = pd.to_datetime(df_sumup.loc[msk, 'end_year'].astype(int).astype(str)+'-01-01')
        msk = df_sumup.start_date == df_sumup.end_date
        df_sumup.loc[msk, 'end_date'] = pd.to_datetime((df_sumup.loc[msk, 'end_year']+1).astype(int).astype(str)+'-01-01')

        df_sumup['smb_mod'] = np.nan
        for i in df_sumup.index:
            df_sumup.loc[i, 'smb_mod'] = df_out.loc[
                df_sumup.loc[i,'start_date']:df_sumup.loc[i,'end_date'],
                'smb_mweq'].sum()

        msk = (df_sumup.end_date-df_sumup.start_date) <= pd.to_timedelta('7 day')
        if msk.any():
            fig = plt.figure(figsize=(12,10))
            gs  = gridspec.GridSpec(2,2, width_ratios=[0.7, 0.3])
            ax=[plt.subplot(g) for g in gs]

            df_sumup_selec = df_sumup.loc[msk]
            ax[0], ax[1] =plt_smb([ax[0], ax[1]], df_sumup_selec)

            df_sumup_selec = df_sumup.loc[~msk]
            ax[2], ax[3] =plt_smb([ax[2], ax[3]], df_sumup_selec)

            ax[0].set_title('Daily/weekly measurements', loc='left')
            ax[2].set_title('Seasonal/annual measurements', loc='left')
        else:
            fig = plt.figure(figsize=(12,5))
            gs  = gridspec.GridSpec(1,2, width_ratios=[0.7, 0.3])
            ax=[plt.subplot(g) for g in gs]
            ax=plt_smb(ax, df_sumup)
        fig.suptitle(c.station)
        fig.savefig(c.output_path+c.RunName+'/'+c.station+'_SMB_evaluation_SUMup2024.png', dpi=120)
        plt.close(fig)
    # except Exception as e:
        # print(c.RunName,e)

def evaluate_accumulation_snowfox(df_in, c):
    # SnowFox
    if c.station in ['KAN_M', 'QAS_M', 'QAS_U','TAS_A','THU_U2']:
        try:
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
            plt.close(fig)
        except Exception as e:
            print(c.RunName,e)
           

def plot_observed_vars(df_obs, df_out, c, var_list = ['t_surf','LRout','LHF','SHF','t_i_10m']):
    try:
        fig = plt.figure(figsize=(12, 17))
        gs = gridspec.GridSpec(len(var_list), 2, width_ratios=[3, 1])

        df_obs = df_obs[~df_obs.index.duplicated(keep='first')]
        df_out = df_out[~df_out.index.duplicated(keep='first')]
        common_idx = df_obs.index.intersection(df_out.index)
        for i, var in enumerate(var_list):
            if var not in df_obs.columns:
                df_obs[var] = np.nan
            ax1 = plt.subplot(gs[i, 0])
            ax2 = plt.subplot(gs[i, 1])
            # first plot
            ME = np.mean(df_out.loc[common_idx, var] - df_obs.loc[common_idx, var])
            RMSE = np.sqrt(np.mean((df_out.loc[common_idx, var] - df_obs.loc[common_idx, var])**2))
            df_obs[var].plot(ax=ax1, label='AWS',marker='.',markersize=2)
            df_out[var].plot(ax=ax1, alpha=0.7, label='SEB model')
            ax1.set_ylabel(var)
            ax1.set_xlim(common_idx.min(), common_idx.max())
            ax1.grid()
            if i == 0:  ax1.set_title(c.station+'\n\n')
            ax1.legend()

            # second plot
            ax2.plot(df_obs.loc[common_idx,var], df_out.loc[common_idx,var],
                     color='k',alpha=0.1,marker='.',ls='None')
            ax2.set_xlabel('AWS')
            ax2.set_ylabel('SEB model')
            common_idx = df_obs.loc[df_obs[var].notnull()].index.intersection(df_out.loc[df_out[var.replace('_uncor','')].notnull()].index)

            try:
                slope, intercept, r_value, p_value, std_err = linregress(
                    df_obs.loc[common_idx, var], df_out.loc[common_idx, var])
                max_val = max(df_obs.loc[common_idx,var].max(), df_out.loc[common_idx,var].max())
                min_val = min(df_obs.loc[common_idx,var].min(), df_out.loc[common_idx,var].min())
                ax2.plot([min_val, max_val], [min_val, max_val], 'k-', label='1:1 Line')
                regression_line = slope * df_obs[var] + intercept
                ax2.plot(df_obs[var], regression_line, 'r-', label='Linear Regression')
            except:
                pass
            if i == 0: ax2.legend(loc='lower right')
            ax2.grid()

            # Annotate with RMSE and ME
            ax2.annotate(f'{var}\nRMSE: {RMSE:.2f}\nME: {ME:.2f}',
                         xy=(1.05, 0.95), xycoords='axes fraction',
                         horizontalalignment='left', verticalalignment='top',
                         fontsize=10, bbox=dict(boxstyle="round,pad=0.3",
                                                edgecolor='black', facecolor='white'))

        fig.savefig(c.output_path+c.RunName+'/SEB_evaluation_vs_AWS.png', dpi=120)
        plt.close(fig)
    except Exception as e:
        print(c.RunName, e)
        
        
def plot_smb_components(df_out, c):
    try:
        fig, ax = plt.subplots(2,1, sharex=True)
        ax[0].plot(df_out.index, df_out.smb_mweq.cumsum(), color='k',lw=2,label='SMB')
        ax[0].fill_between(df_out.index, df_out.sublimation_mweq.cumsum()*0,
                           -df_out.sublimation_mweq.cumsum(), label='sublimation_mweq')
        ax[0].fill_between(df_out.index, -df_out.sublimation_mweq.cumsum(),
                           (df_out.snowfall_mweq - df_out.sublimation_mweq).cumsum(), label='snowfall_mweq')
        ax[0].fill_between(df_out.index, (df_out.snowfall_mweq - df_out.sublimation_mweq).cumsum(),
               (df_out.snowfall_mweq + df_out.rainfall_mweq - df_out.sublimation_mweq).cumsum(),
               label='rainfall_mweq')
        ax[0].fill_between(df_out.index, df_out.snowfall_mweq.cumsum()*0,
                           -df_out.runoff_mweq.cumsum(), label='Runoff')
        ax[0].legend()

        ax[1].plot(df_out.index, df_out.melt_mweq.cumsum(), label='melt')
        ax[1].plot(df_out.index, df_out.runoff_mweq.cumsum(), color='tab:red', label='runoff')
        ax[1].plot(df_out.index, df_out.refreezing_mweq.cumsum(), label='refreezing')
        ax[1].plot(df_out.index, df_out.snowthick, label='snowthickness')
        ax[1].legend()

        fig.savefig(c.output_path+c.RunName+'/'+c.station+'_SMB.png', dpi=120)
        plt.close(fig)
    except Exception as e:
        print(c.RunName, e)