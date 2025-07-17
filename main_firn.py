# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""

import lib_subsurface as sub
import lib_initialization as ini
import lib_io as io
import lib_plot as lpl
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import os
import pandas as pd

from joblib import Parallel, delayed
from lib_initialization import ImportConst, load_json
from progressbar import progressbar

def run_GEUS_model_opt(site):
    # Read paths for weather input and output file
    parameters = load_json()
    output_path_firn = str(parameters['output_path_firn'])
    weather_station = str(parameters['weather_station'])
    station_info = str(parameters['constants']['station_info']['path'])
    weather_data_input_path_unformatted = str(parameters['weather_data']['weather_input_path'])
    weather_data_input_path = weather_data_input_path_unformatted.format(site) 

    # Define constant c
    c = ImportConst()
    c.station = site
    c.rh2oice = c.rho_water / c.rho_ice
    c.zdtime = 3600
    NumLayer = 200
    c.num_lay = NumLayer
    c.z_max = 50
    c.dz_ice = c.z_max / NumLayer
    c.verbose = 1
    c.OutputFolder = output_path_firn
    # Load station info
    df_info = pd.read_csv(station_info, sep=";")
    c.Tdeep = (
        df_info.loc[
            df_info["station name"] == c.station, "deep firn temperature (degC)"
        ].values[0]
        + 273.15
    )
    # c.lim_new_lay = c.accum_AWS/c.new_lay_frac;
    c.lim_new_lay = 0.02
    c.rho_fresh_snow = 315
   
    # New code:
    # Load weather input data
    df_aws = io.load_promice(weather_data_input_path)[:6000]
    df_aws = df_aws.set_index("time").resample("H").mean()
    time = df_aws.index.values

    # Initialize variables
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

    if 'acc_subl_mmweq' not in df_aws.columns:
        print("You do not have preprocessed input data in the right format. Values from the SEB model will be used.")
        # Read parameter values from main SEB firn, UPDATE
        path = r'C:\Users\brink\Documents\Exjobb\GEUS-SEB-firn-model\Input\Weather data\param_from_old_SEB_' + site + '_200_lay.csv'
        df_SEB_output = pd.read_csv(path, index_col=False)
        snowfall = df_SEB_output['snowfall']
        Tsurf = df_SEB_output['Tsurf']
        sublimation_mweq = df_SEB_output['sublimation_mweq']
        melt_mweq = df_SEB_output['melt_mweq']
        Tsurf_out = Tsurf.copy() * np.nan
        net_accum = [snowfall[i] - sublimation_mweq[i] for i in range(min(len(snowfall), len(sublimation_mweq)))]
    else:
        print("The data contains preprocessed information, the SEB model will not be called.")
        # Original code
        # ts = df_aws.Tsurf_K.values
        # ts_out = ts.copy() * np.nan
        # net_accum = df_aws.acc_subl_mmweq.values / 1000
        # melt = df_aws.melt_mmweq.values / 1000

    for i in progressbar(range(0, len(time))):
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
            Tsurf_out[i],
            grndc[:, i],
            grndd[:, i],
            pgrndcapc[i],
            pgrndhflx[i],
            dH_comp[i],
            snowbkt[i],
            compaction[:, i],
        ) = sub.subsurface_opt(
            Tsurf[i].copy(),
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
            melt_mweq[i].copy(),
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

    # #Plotting from main_seb_firn.py
    plt.close("all")
    lpl.plot_summary(df_aws, c, 'input_summary', var_list = ['RelativeHumidity1','RelativeHumidity2'])
    lpl.plot_summary(df_SEB_output, c, 'SEB_output')

    print("File created: ", c.RunName)
    return c.RunName

    #Return for testing
    # return (snowc, snic, slwc, tsoil, zrfrz, 
    #         rhofirn, zsupimp, dgrain, grndc, 
    #         grndd, compaction, zrogl, Tsurf_out,
    #         pgrndcapc, pgrndhflx, dH_comp, snowbkt)


# Run new main_firn, parallel
def run_main_firn_parallel(site_list):
    num_cores = multiprocessing.cpu_count()
    run_name_list = Parallel(n_jobs=num_cores, verbose=10)(delayed(run_GEUS_model_opt)(site) for site in site_list)
   
    return(run_name_list)


if __name__ == "__main__":
    site_list = [
         "KAN_M",
         "KAN_U"       
    ] 
    # Run new main_firn, parallelised
    run_main_firn_parallel(site_list)

    # Calls to original code
    # Run old main_firn
    # run_main_firn_old(site_list)



# %% debbuging
# i = 0
# (pts, pgrndc, pgrndd, pslwc, psnic, psnowc, prhofirn, ptsoil, pdgrain, zsn, zraind, zsnmel, pTdeep, psnowbkt, c) = (ts[i].copy(), grndc[:,i-1].copy(), grndd[:, i-1].copy(),
#                             slwc[:, i-1].copy(), snic[:, i-1].copy(), snowc[:, i-1].copy(),
#                             rhofirn[:, i-1].copy(), tsoil[:, i-1].copy(), dgrain[:, i-1].copy(),
#                             net_accum[i].copy(), 0, melt[i].copy(), c.Tdeep, snowbkt[i-1].copy(), c)



# Old script (modified to get it running), used for comparison and testing:
def run_GEUS_model_old(site, filename):
    # New code 
    # Read paths for weather input and output file
    parameters = load_json()
    output_path_firn = str(parameters['output_path_firn'])
    # weather_station = str(parameters['weather_station'])
    station_info = str(parameters['constants']['station_info']['path'])
    weather_data_input_path_unformatted = str(parameters['weather_data']['weather_input_path'])
    weather_data_input_path = weather_data_input_path_unformatted.format(site)

    c = ini.ImportConst()
    c.station = site
    c.rh2oice = c.rho_water / c.rho_ice
    c.zdtime = 3600
    NumLayer = 200
    c.num_lay = NumLayer
    c.z_max = 50
    c.dz_ice = c.z_max / NumLayer
    c.verbose = 1
    #c.OutputFolder = "C://Data_save/Output firn model"
    c.OutputFolder = output_path_firn + '/old_script'

    df_info = pd.read_csv(station_info, sep=";")
    #df_info = pd.read_csv("Input/Constants/StationInfo.csv", sep=";")
    c.Tdeep = (
        df_info.loc[
            df_info["station name"] == c.station, "deep firn temperature (degC)"
        ].values[0]
        + 273.15
    )

    # c.lim_new_lay = c.accum_AWS/c.new_lay_frac;
    c.lim_new_lay = 0.02
    c.rho_fresh_snow = 315

    #Old code for loading data
    # df_aws = pd.read_csv('./Input/'+c.station+'_high-res_meteo.csv', sep=';')
    # df_aws = pd.read_csv(filename, sep=";")
    # df_aws.time = pd.to_datetime(df_aws.time)
    # df_aws = df_aws.set_index("time").resample("H").nearest()

    # New code for loading data:
    # Load weather input data
    df_aws = io.load_promice(weather_data_input_path)[:6000]
    df_aws = df_aws.set_index("time").resample("H").mean()
    time = df_aws.index.values

    #print("Number of NaN: ")
    #print(df_aws.isnull().sum())
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

    # New code added in order to load data from SEB
    if 'acc_subl_mmweq' not in df_aws.columns:
        #print("You do not have preprocessed input data in the right format. Values from the SEB model will be used.")

        # Read parameter values from main SEB firn, UPDATE
        path = r'C:\Users\brink\Documents\Exjobb\GEUS-SEB-firn-model\Input\Weather data\param_from_old_SEB_' + site + '_200_lay.csv'
        df_SEB_output = pd.read_csv(path, index_col=False)

        snowfall = df_SEB_output['snowfall']
        Tsurf = df_SEB_output['Tsurf']
        sublimation_mweq = df_SEB_output['sublimation_mweq']
        melt_mweq = df_SEB_output['melt_mweq']

        Tsurf_out = Tsurf.copy() * np.nan
        net_accum = [snowfall[i] - sublimation_mweq[i] for i in range(min(len(snowfall), len(sublimation_mweq)))]
    
    else:
        print("The data contains preprocessed information, the SEB model will not be called.")
        # Original code
        # ts = df_aws.Tsurf_K.values
        # ts_out = ts.copy() * np.nan
        # net_accum = df_aws.acc_subl_mmweq.values / 1000
        # melt = df_aws.melt_mmweq.values / 1000

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
            #ts_out[i],
            Tsurf_out[i],
            grndc[:, i],
            grndd[:, i],
            pgrndcapc[i],
            pgrndhflx[i],
            dH_comp[i],
            snowbkt[i],
            compaction[:, i],
        ) = sub.subsurface_opt(
            Tsurf[i].copy(),
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
            melt_mweq[i].copy(),
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

    # # Plotting from main_seb_firn.py
    # plt.close("all")
    # lpl.plot_summary(df_aws, c, 'input_summary', var_list = ['RelativeHumidity1','RelativeHumidity2'])
    # lpl.plot_summary(df_SEB_output, c, 'SEB_output')

    print("File created: ",c.RunName)
    return c.RunName

    #For test file:
    # return (snowc, snic, slwc, tsoil, zrfrz, 
    #         rhofirn, zsupimp, dgrain, grndc, 
    #         grndd, compaction, zrogl, Tsurf_out,
    #         pgrndcapc, pgrndhflx, dH_comp, snowbkt)

# Run original version of code
def run_main_firn_old(site_list):
    run_name_list = []
    for site in site_list:
        print(site)
        filename = "Firn viscosity/Input files/" + site + ".csv"
        run_name = run_GEUS_model_old(site, filename)
        run_name_list.append(run_name)
    return run_name_list