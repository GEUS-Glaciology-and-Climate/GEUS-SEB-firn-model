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
    filename = "Output/" + run_name + "/" + site + "_" + var_name + ".nc"

    ds = nc.Dataset(filename)

    time_org = np.asarray(ds["time"][:])
    time = np.datetime64("1900-01-01T00") + np.round(
        (time_org - 1) * 24 * 3600
    ) * np.timedelta64(1, "s")

    depth = np.asarray(ds["depth"][:])

    if ~zero_surf:
        surface_height = -depth[-1, :] + depth[-1, 0]
        depth = depth + surface_height

    var_array = np.asarray(ds[var_name][:])
    if var_name == "slwc":
        var_array = var_array * 1000 / depth
    # now in mm / m3

    depth = np.concatenate((depth[:1, :] * 0, depth), axis=0)
    depth = np.concatenate((depth, depth[:, -1:]), axis=1)

    time = np.concatenate((time, time[-1:]))

    time_grid = np.expand_dims(time, 0)
    time_grid = np.repeat(time_grid, depth.shape[0], axis=0)

    myFmt = mdates.DateFormatter("%Y-%m-%d")

    label_list = dict(
        slwc="Liquid water content (mm m$^{-3}$)",
        density_bulk="Bulk density (kg m$^{-3}$)",
        rhofirn="Firn density (kg m$^{-3}$)",
        T_ice="Subsurface temperature ($^{o}C$)",
    )

    cmap_list = dict(
        slwc="gist_ncar_r", density_bulk="Blues", rhofirn="Blues", T_ice="magma"
    )

    fig, ax = plt.subplots(1, 1, figsize=(15, 30))
    plt.subplots_adjust(left=0.07, right=0.99, top=0.95, bottom=0.1, hspace=0.2)
    count = 0
    fig.suptitle(site)

    # plotting firn model
    im = ax.pcolormesh(
        time_grid[:, 0::6],
        depth[:, 0::6],
        var_array[:, 0::6][:, :-1],
        cmap=cmap_list[var_name],
        vmin=np.percentile(var_array, 5),
        vmax=np.percentile(var_array, 95),
    )
    if not zero_surf:
        ax.plot(time[0::6], surface_height[0::6], linewidth=2, color="k")
    plt.colorbar(im, label=label_list[var_name], ax=ax)
    if len(ylim) > 0:
        ax.set_ylim(ylim)
    else:
        ax.invert_yaxis()

    ax.xaxis.set_major_formatter(myFmt)
    ax.set_ylabel("Depth (m)")

    fig.savefig("Output/" + run_name + "/" + site + "_" + var_name + ".png")
    return fig, ax


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
        comp_tot = np.sum(comp_mod[ind_next:])

        # plus compaction within the layer where the horizon is
        comp = (
            (depth_mod[ind_next] - depth_hor[i])
            / (depth_mod[ind_next] - depth_mod[ind_next - 1])
            * comp_mod[ind_next - 1]
        )

        comp_tot = comp_tot + comp

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
        fig, ax = plt.subplots(6, 1, sharex=True, figsize=(15, 10))
        plt.subplots_adjust(
            left=0.1, right=0.9, top=0.97, bottom=0.1, wspace=0.2, hspace=0.05
        )
        return fig, ax

    if not var_list:
        var_list = df.columns
    fig, ax = new_fig()
    count = 0
    count_fig = 0
    for i, var in enumerate(var_list):
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
        count = count + 1

        if count == 6:
            ax[0].set_title(c.station)
            plt.savefig(
                c.OutputFolder + "/" + c.RunName + "/" + "summary_" + str(count_fig),
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
            c.OutputFolder + "/" + c.RunName + "/" + filetag + "_" + str(count_fig),
            bbox_inches="tight",
        )
