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

def plot_var(site, run_name, var_name, ylim=[], zero_surf=True):
    print('plotting',var_name, 'from',run_name)
    filename = "Output/" + run_name + "/" + site + "_" + var_name + ".nc"
    ds = xr.open_dataset(filename).transpose()
    ds = ds.resample(time='6H').nearest()
    
    if not zero_surf:
        ds['surface_height'] = -ds.depth.isel(level=-1) + ds.depth.isel(level=-1).isel(time=0)
        ds['depth'] = ds.depth + ds.surface_height
    
    if var_name == "slwc":
        # change unit to mm / m3
        ds[var_name] = ds[var_name] * 1000 / ds.depth
    
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
                   cmap = cmap, vmin=vmin, vmax=vmax
                   )
    
    if not zero_surf:
        ax.plot(ds.time, ds.surface_height, linewidth=2, color="k")
        
    plt.colorbar(im, label=label, ax=ax)
    ax.invert_yaxis()
    if ylim:
        if len(ylim)==1: ax.set_ylim(ylim, ax.get_ylim()[1])
        if len(ylim)==2: ax.set_ylim(np.max(ylim), np.min(ylim))
    ax.set_ylabel("Depth (m)")
    
    fig.savefig("output/" + run_name + "/" + site + "_" + var_name + ".png")
    return fig, ax


def evaluate_compaction(site, run_name):
    filename = "Output/" + run_name + "/" + site + "_compaction.nc"
    ds = nc.Dataset(filename)
    compaction = ds["compaction"][:]
    time_org = np.asarray(ds["time"][:])
    time = np.datetime64("1900-01-01T00") + np.round(
        (time_org - 1) * 24 * 3600
    ) * np.timedelta64(1, "s")

    depth_act = np.asarray(ds["depth"][:])

    H_surf = depth_act[-1, :] - depth_act[-1, 0]

    df_comp_info = pd.read_csv(
        "Firn viscosity/Compaction_Instrument_Metadata.csv"
    ).set_index("sitename")

    df_comp_info = (
        df_comp_info.loc[
            site,
            [
                "instrument_ID",
                "installation_daynumber_YYYYMMDD",
                "borehole_top_from_surface_m",
                "borehole_bottom_from_surface_m",
            ],
        ]
        .reset_index(drop=True)
        .set_index("instrument_ID")
    )

    df_comp = pd.read_csv("Firn viscosity/borehole_shortening_m.csv")
    df_comp.date = pd.to_datetime(df_comp.date)
    df_comp = df_comp.set_index(["instrument_id", "date"])

    fig1, ax = plt.subplots(1, 1)  # plot_var(site, run_name, 'density_bulk')
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

        depth_1 = track_horizon(
            time, H_surf, depth_act, compaction, date_start, depth_top, step=12
        )
        depth_2 = track_horizon(
            time, H_surf, depth_act, compaction, date_start, depth_bot, step=12
        )

        ax.plot(
            time,
            depth_1 - H_surf,
            color=cmap(i / len(df_comp_info.index)),
            label="_no_legend_",
        )
        ax.plot(
            time,
            depth_2 - H_surf,
            color=cmap(i / len(df_comp_info.index)),
            label="Instrument " + str(ID),
        )

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
    fig2.text(
        0.03,
        0.5,
        "Compaction rate (m d$^{-1}$)",
        ha="center",
        va="center",
        rotation="vertical",
    )

    fig1.savefig("Output/" + run_name + "/" + site + "_compaction_1.png")
    fig2.savefig("Output/" + run_name + "/" + site + "_compaction_2.png")


def plot_summary(df, c, filetag="summary", var_list=None):
    def new_fig():
        fig, ax = plt.subplots(7, 1, sharex=True, figsize=(15, 10))
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
        
        # if var == "L":    #Changing the y-axis for L
        #     ax[count].set_ylim((-30000, 30000))


        count = count + 1

        if count == len(ax):
            ax[0].set_title(c.station)
            plt.savefig(
                c.output_path + "/" + c.RunName + "/" + "summary_" + str(count_fig),
                bbox_inches="tight",
            )
            count_fig = count_fig + 1
            fig, ax = new_fig()
            count = 0
    if count < 6:
        count = count - 1
        ax[count].xaxis.set_tick_params(which="both", labelbottom=True)

        for k in range(count + 1, len(ax)):
            ax[k].set_axis_off()
        ax[0].set_title(c.station)
        
        plt.savefig(
            c.output_path + "/" + c.RunName + "/" + filetag + "_" + str(count_fig),
            bbox_inches="tight",
        )
