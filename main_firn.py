# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import pandas as pd
import numpy as np
import lib_subsurface as sub
import lib_initialization as ini
from progressbar import progressbar
import xarray as xr
import lib_io as io
import os


def run_GEUS_model(site, filename):
    c = ini.ImportConst()
    c.station = site
    c.rh2oice = c.rho_water / c.rho_ice
    c.zdtime = 3600
    NumLayer = 100
    c.num_lay = NumLayer
    c.z_max = 50
    c.dz_ice = c.z_max / NumLayer
    c.verbose = 1
    c.OutputFolder = "C://Data_save/Output firn model"

    df_info = pd.read_csv("Input/Constants/StationInfo.csv", sep=";")
    c.Tdeep = (
        df_info.loc[
            df_info["station name"] == c.station, "deep firn temperature (degC)"
        ].values[0]
        + 273.15
    )

    # c.lim_new_lay = c.accum_AWS/c.new_lay_frac;
    c.lim_new_lay = 0.02
    c.rho_fresh_snow = 315

    # df_aws = pd.read_csv('./Input/'+c.station+'_high-res_meteo.csv', sep=';')
    df_aws = pd.read_csv(filename, sep=";")
    df_aws.time = pd.to_datetime(df_aws.time)
    df_aws = df_aws.set_index("time").resample("H").nearest()

    print("Number of NaN: ")
    print(df_aws.isnull().sum())

    time = df_aws.index.values

    (
        rhofirn,
        snowc,
        snic,
        slwc,
        dgrain,
        tsoil,
        grndc,
        grndd,
        compaction,
        zrfrz,
        zsupimp,
        ts,
        zrogl,
        pgrndcapc,
        pgrndhflx,
        dH_comp,
        snowbkt,
    ) = ini.IniVar(time, c)

    ts = df_aws.Tsurf_K.values
    ts_out = ts.copy() * np.nan
    net_accum = df_aws.acc_subl_mmweq.values / 1000
    melt = df_aws.melt_mmweq.values / 1000

    for i in progressbar(range(0, len(time))):
        # print(i)
        (
            snowc[:, i],
            snic[:, i],
            slwc[:, i],
            tsoil[:, i],
            zrfrz[:, i],
            rhofirn[:, i],
            zsupimp[:, i],
            dgrain[:, i],
            zrogl[i],
            ts_out[i],
            grndc[:, i],
            grndd[:, i],
            pgrndcapc[i],
            pgrndhflx[i],
            dH_comp[i],
            snowbkt[i],
            compaction[:, i],
        ) = sub.subsurface(
            ts[i].copy(),
            grndc[:, i - 1].copy(),
            grndd[:, i - 1].copy(),
            slwc[:, i - 1].copy(),
            snic[:, i - 1].copy(),
            snowc[:, i - 1].copy(),
            rhofirn[:, i - 1].copy(),
            tsoil[:, i - 1].copy(),
            dgrain[:, i - 1].copy(),
            net_accum[i].copy(),
            0,
            melt[i].copy(),
            c.Tdeep,
            snowbkt[i - 1].copy(),
            c,
        )

    # writing output
    thickness_act = snowc * (c.rho_water / rhofirn) + snic * (c.rho_water / c.rho_ice)
    depth_act = np.cumsum(thickness_act, 0)
    density_bulk = (snowc + snic) / (snowc / rhofirn + snic / c.rho_ice)

    c.RunName = c.station + "_" + str(c.num_lay) + "_layers"
    i = 0
    succeeded = 0

    

    while succeeded == 0:
        try:
            os.mkdir(c.OutputFolder + "/" + c.RunName)
            succeeded = 1
        except:
            if i == 0:
                c.RunName = c.RunName + "_" + str(i)
            else:
                c.RunName = c.RunName[: -len(str(i - 1))] + str(i)
            i = i + 1

    io.write_2d_netcdf(snowc, "snowc", depth_act, time, c)
    io.write_2d_netcdf(snic, "snic", depth_act, time, c)
    io.write_2d_netcdf(slwc, "slwc", depth_act, time, c)
    io.write_2d_netcdf(density_bulk, "density_bulk", depth_act, time, c)
    io.write_2d_netcdf(rhofirn, "rhofirn", depth_act, time, c)
    io.write_2d_netcdf(tsoil, "T_ice", depth_act, time, c)
    # io.write_2d_netcdf(rhofirn, 'rho_firn_only', depth_act, time, RunName)
    # io.write_2d_netcdf(rfrz, 'rfrz', depth_act, time, RunName)
    # io.write_2d_netcdf(dgrain, 'dgrain', depth_act, time, RunName)
    io.write_2d_netcdf(compaction, "compaction", depth_act, time, c)
    print(c.RunName)
    return c.RunName


if __name__ == "__main__":
    site_list = [
        "CEN",
        "FA",
    ]  # ['EGP', 'DYE-2','CP1','Saddle', 'Summit', 'KAN_U', 'NASA-SE']
    # other sites to add:   'FA'

   
    for site in site_list:
        print(site)
        filename = "Firn viscosity/Input files/" + site + ".csv"
        run_name = run_GEUS_model(site, filename)
        print(run_name)
    
    run_GEUS_model(site, filename)
# %% debbuging
# i = 0
# (pts, pgrndc, pgrndd, pslwc, psnic, psnowc, prhofirn, ptsoil, pdgrain, zsn, zraind, zsnmel, pTdeep, psnowbkt, c) = (ts[i].copy(), grndc[:,i-1].copy(), grndd[:, i-1].copy(),
#                             slwc[:, i-1].copy(), snic[:, i-1].copy(), snowc[:, i-1].copy(),
#                             rhofirn[:, i-1].copy(), tsoil[:, i-1].copy(), dgrain[:, i-1].copy(),
#                             net_accum[i].copy(), 0, melt[i].copy(), c.Tdeep, snowbkt[i-1].copy(), c)
