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
import os
import json

import lib_initialization as ini
import lib_seb_smb_model as seb
import lib_io as io

import __init__

def run_SEB_firn():
    # Create struct c with all constant values
    c = setConstants()

    with open("parameters.json") as parameter_file:
        parameters = json.load(parameter_file)
    output_path = str(parameters['output_path'])
    weather_data_input_path = str(parameters['weather_data']['weather_input_path'])

    # DataFrame with the weather data is created
    df_aws = io.load_promice(weather_data_input_path)[:6000]
    df_aws = df_aws.set_index("time").resample("H").mean()

    # DataFrame for the surface is created, indexed with time from df_aws
    df_surface = pd.DataFrame()
    df_surface["time"] = df_aws.index
    df_surface = df_surface.set_index("time")

    # The surface values are received
    (
        df_surface["L"],
        df_surface["LHF"],
        df_surface["SHF"],
        df_surface["theta_2m"],
        df_surface["q_2m"],
        df_surface["ws_10m"],
        df_surface["Re"],
        df_surface["melt_mweq"],
        df_surface["sublimation_mweq"],
        snowc,
        snic,
        slwc,
        T_ice,
        zrfrz,
        rhofirn,
        zsupimp,
        dgrain,
        zrogl,
        Tsurf,
        grndc,
        grndd,
        pgrndcapc,
        pgrndhflx,
        dH_comp,
        snowbkt,
        compaction,
    ) = seb.HHsubsurf(df_aws, c)

    thickness_act = snowc * (c.rho_water / rhofirn) + snic * (c.rho_water / c.rho_ice)
    depth_act = np.cumsum(thickness_act, 0)
    density_bulk = (snowc + snic) / (snowc / rhofirn + snic / c.rho_ice)

    # Writing output
    c.RunName = c.station + "_" + str(c.num_lay) + "_layers"
    i = 0
    succeeded = 0
    while succeeded == 0:
        try:
            os.mkdir(output_path + c.RunName)
            succeeded = 1
        except:
            if i == 0:
                c.RunName = c.RunName + "_" + str(i)
            else:
                c.RunName = c.RunName[: -len(str(i - 1))] + str(i)
            i = i + 1
    c.OutputFolder = output_path

    # io.write_2d_netcdf(snowc, 'snowc', depth_act, df_aws.index, c)
    # io.write_2d_netcdf(snic, 'snic', depth_act, df_aws.index, c)
    # io.write_2d_netcdf(slwc, 'slwc', depth_act, df_aws.index, c)
    io.write_2d_netcdf(density_bulk, 'density_bulk', depth_act, df_aws.index, c)
    # io.write_2d_netcdf(rhofirn, 'rhofirn', depth_act, df_aws.index, c)
    io.write_1d_netcdf(df_surface, c)
    # io.write_2d_netcdf(tsoil, 'T_ice', depth_act, df_aws.index, c)
    # io.write_2d_netcdf(rhofirn, 'rho_firn_only', depth_act, df_aws.index, RunName)
    # io.write_2d_netcdf(rfrz, 'rfrz', depth_act, df_aws.index, RunName)
    # io.write_2d_netcdf(dgrain, 'dgrain', depth_act, df_aws.index, RunName)
    # io.write_2d_netcdf(compaction, 'compaction', depth_act, df_aws.index, RunName)

   
# Constant definition
# All constant values are defined in a set of csv file in the Input folder.
# They can be modiefied there or by giving new values in the "param" variable. 
# The values of the constant given in param will overright the ones extracted 
# from the csv files. The fieldnames in param should be the same as is c.
def setConstants():
    c = ini.ImportConst("parameters.json")
    c.station = "KAN_M"
    c.elev = 2000
    c.rh2oice = c.rho_water / c.rho_ice
    c.zdtime = 3600
    c.ElevGrad = 0.1
    c.z_max = 50
    c.dz_ice = 1
    NumLayer = int(c.z_max / c.dz_ice)
    c.num_lay = NumLayer
    c.verbose = 1
    c.Tdeep = 250.15
    # c.lim_new_lay = c.accum_AWS/c.new_lay_frac;
    c.lim_new_lay = 0.02
    c.rho_fresh_snow = 315
    c.snowthick_ini = 1
    c.dz_ice = 1
    c.z_ice_max = 50
    c.dt_obs = 3600   
    c.kappa = np.float64(c.kappa)
    return c 


#run_SEB_firn()


    #%% Plot output
    #import lib_plot as lpl
    #import matplotlib.pyplot as plt

    #plt.close("all")
    # lpl.plot_summary(df_aws, c, 'input_summary', var_list = ['RelativeHumidity1','RelativeHumidity2'])
    # lpl.plot_summary(df_surface, c, 'SEB_output')
    #lpl.plot_var(c.station, c.RunName, "density_bulk", ylim=(10, -5), zero_surf=False)

