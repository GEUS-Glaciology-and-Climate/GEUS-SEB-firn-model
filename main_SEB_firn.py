# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import lib.io as io
from lib.initialization import ImportConst
from lib.seb_smb_model import GEUS_model
import os
import time
import plot_output as po
import xarray as xr
import matplotlib.pyplot as plt
import traceback

# %%
def run_SEB_firn(station='FA-13', silent=False):
    start_time = time.time()
    # importing standard values for constants
    c = ImportConst()
    c.station = station
    c.spin_up = False
    c.verbose = 1
    if silent or c.spin_up:
        c.verbose = 0
    c.use_spin_up_init = True

    # default setting
    # c.surface_input_path = f"/data/CARRA/extracted/list_pixels_minimal_grid/{station}.nc"
    c.surface_input_path = f"./input/weather data/CARRA_at_AWS/{station}.nc"
    c.surface_input_driver = "CARRA"
    # c.output_path = '/data/CARRA-SMB/list_pixels_minimal_grid/'
    c.output_path = './output/new/'
    # c.spin_up_path = '/data/CARRA-SMB/spin up 3H/'
    c.spin_up_path = './output/spin up 3H/'
    # c.initial_state_folder_path = '/data/CARRA-SMB/spin up 3H/'
    c.initial_state_folder_path = './output/spin up 3H/'

    c.freq = '3h'
    if c.surface_input_driver=='CARRA' and c.freq == 'h':
        resample=True
    else:
        resample=False

    c.num_lay = 100
    c.lim_new_lay = 0.05

    # defining run name
    if c.spin_up:
        print('######### spin-up run ##########')
        c.output_path = '/data/CARRA-SMB/spin up 3H/'

    c.RunName = c.station + "_" + str(c.num_lay) + "_layers_"+c.freq

    # if it is a spin up, then checking whether the pckl file is already available
    if c.spin_up and os.path.isfile(c.output_path+c.RunName+'/'+c.station+'_final.pkl'):
        print(c.RunName, 'already exists, skipping')
        return

    # for all files, we can check if the run folder already exist and skip if necessary
    try:
        os.mkdir(c.output_path + c.RunName)
    except Exception as e:
        if os.path.isfile(c.output_path+c.RunName+'/'+c.station+'_rhofirn.nc'):
            print(c.RunName, 'already exists, redoing')

    # if it is not a spin up and there is no spinup run available then skipping
    if not c.spin_up and not os.path.isfile(c.spin_up_path+c.RunName+'/'+c.station+'_final.pkl'):
        print(c.RunName, 'does not have a spin up run yet, skipping')
        return

    # loading input data
    try:
        df_in, c = io.load_surface_input_data(c, resample=resample)
    except Exception as e:
        print(c.station, e); traceback.print_exc()

    freq = pd.infer_freq(df_in.index)
    if freq=='h': freq = '1h'
    c.delta_time = pd.to_timedelta(freq).total_seconds()
    c.zdtime = c.delta_time

    if c.altitude<1500:
        c.new_bottom_lay=30
    # assigning constants specific to this simulation
    c.snowthick_ini = 0.1
    c.z_max = 70

    # assgning deep temperature
    # ds_T10m = xr.open_dataset('../../Data/Firn temperature/output/T10m_prediction.nc')
    with xr.open_dataset('input/T10m_prediction.nc') as ds_T10m:
        c.Tdeep = ds_T10m.sel(latitude = c.latitude,
                        longitude = c.longitude,
                        method='nearest').sel(time=slice('1980','1990')
                      ).T10m.mean().item() + 273.15

    if np.isnan(c.Tdeep): c.Tdeep = 273.15

    # c.lim_new_lay = c.accum_AWS/c.new_lay_frac;

    print(station, c.Tdeep, 'start/end', df_in.index[0], df_in.index[-1])
    # DataFrame for the surface is created, indexed with time from df_aws
    df_out = pd.DataFrame()
    df_out["time"] = df_in.index
    df_out = df_out.set_index("time")

    # %% Running model
    if c.verbose>0: print(c.RunName,'reading inputs took %0.03f sec'%(time.time() -start_time))
    start_time = time.time()
    # The surface values are received
    try:
        (
            df_out["L"], df_out["LHF"], df_out["SHF"], df_out["theta_2m"],
            df_out["q_2m"], df_out["ws_10m"], df_out["Re"], df_out["melt_mweq"],
            df_out["sublimation_mweq"], df_out["SRin"], df_out["SRout"],
            df_out["LRin"], df_out["LRout_mdl"],
            snowc, snic, slwc, T_ice,
            zrfrz, rhofirn, zsupimp, dgrain,
            df_out['zrogl'], Tsurf,
            grndc, grndd, pgrndcapc, pgrndhflx,
            dH_comp, df_out["snowbkt"], compaction, df_out['snowthick']
        ) = GEUS_model(df_in, c)
    except Exception as e:
        print(c.station, e); traceback.print_exc()
        with open(f"{c.output_path}/{c.RunName}/error.txt", "a+") as text_file:
            text_file.write(f"{c.RunName}\n")  # Write RunName
            text_file.write(str(e) + "\n")  # Write error message
            text_file.write(traceback.format_exc() + "\n")  # Write full traceback
        return

    if c.verbose>0: print('\n',c.RunName,'SEB firn model took %0.03f sec'%(time.time() -start_time))
    start_time = time.time()

    # %% Writing output
    start_time = time.time()

    if not c.spin_up:
        c.write(c.output_path + "/" + c.RunName + "/constants.csv")
        thickness_act = snowc * c.rho_water / rhofirn + snic * c.rho_water / c.rho_ice
        depth_act = np.cumsum(thickness_act, 0)
        # density_bulk = (snowc + snic) / (snowc / rhofirn + snic / c.rho_ice)

        depth_lowest_layer = np.insert(np.diff(depth_act[-1,:]),0,0)
        depth_lowest_layer[np.abs(depth_lowest_layer)<=3]=0
        surface_height = (depth_act[-1,:] -depth_act[-1,0])-depth_lowest_layer.cumsum()
        df_out['surface_height'] = surface_height
        df_out['snowfall_mweq'] = df_in.Snowfallmweq
        df_out['rainfall_mweq'] = df_in.Rainfallmweq
        df_out['refreezing_mweq'] = zrfrz.sum(axis=0)
        df_out = df_out.rename(columns={'zrogl': 'runoff_mweq'})
        df_out['smb_mweq'] = df_out.snowfall_mweq + df_out.rainfall_mweq - df_out.runoff_mweq - df_out.sublimation_mweq
        # ds['surface_height'] = (ds.depth.isel(level=-1)
        #         -ds.depth.isel(level=-1).isel(time=0)
        #         -(ds.depth.isel(level=-1).diff(dim='time')
        #           .where(ds.depth.isel(level=-1)
        #                  .diff(dim='time')>6,0).cumsum()))
        # ds['surface_height'].values[0] = 0

        io.write_1d_netcdf(df_out, c)
        io.write_2d_netcdf(snic, 'snic', depth_act, df_in.index, c)
        io.write_2d_netcdf(snowc, 'snowc', depth_act, df_in.index, c)
        io.write_2d_netcdf(slwc, 'slwc', depth_act, df_in.index, c)
        io.write_2d_netcdf(rhofirn, 'rhofirn', depth_act, df_in.index, c)
        # io.write_2d_netcdf(density_bulk, 'density_bulk', depth_act, df_in.index, c)
        io.write_2d_netcdf(T_ice, 'T_ice', depth_act, df_in.index, c)
        # io.write_2d_netcdf(rfrz, 'rfrz', depth_act, df_in.index, RunName)
        # io.write_2d_netcdf(dgrain, 'dgrain', depth_act, df_in.index, c)
        # io.write_2d_netcdf(compaction, 'compaction', depth_act, df_in.index, c)

        if c.verbose>0: print(c.RunName,'writing output files took %0.03f sec'%(time.time() -start_time))
        start_time = time.time()

        # Plot output
        # try:
        po.main(c.output_path, c.RunName)
        # except Exception as e:
        #     print(c.RunName,e); traceback.print_exc()
        # print(c.RunName,'plotting took %0.03f sec'%(time.time() -start_time))
    plt.close('all')
    start_time = time.time()
    return

import os
import multiprocessing

def run_on_core(station, core):
    """Set CPU affinity and run run_SEB_firn for the given station."""
    print(f"Core {core} processing station: {station}")
    os.system(
        f"taskset -c {core} python3 -c 'import main_SEB_firn; "
        f"main_SEB_firn.run_SEB_firn(\"{station}\")'"
    )
    print(f"Core {core} finished station: {station}")

def worker(task_queue):
    """Worker process that executes tasks sequentially on its assigned core."""
    core = task_queue.get()  # Get the assigned core
    while True:
        station = task_queue.get()  # Get the next station to process
        if station is None:
            break  # Stop worker when None is received
        run_on_core(station, core)

def standard_run_parallel(station_list):
    # max_core_usage =  multiprocessing.cpu_count()-1 # all cores except one
    max_core_usage = 21  # Limit to 6 cores
    num_cores = min(len(station_list), max_core_usage)
    task_queues = [multiprocessing.Queue() for _ in range(num_cores)]
    processes = []

    # Start workers with dedicated cores
    for core, task_queue in enumerate(task_queues):
        task_queue.put(core)  # Assign core number to worker
        p = multiprocessing.Process(target=worker, args=(task_queue,))
        p.start()
        processes.append((p, task_queue))

    # Distribute tasks in a round-robin fashion
    for i, station in enumerate(station_list):
        task_queues[i % num_cores].put(station)

    # Send termination signal (None) to workers
    for _, task_queue in processes:
        task_queue.put(None)

    # Wait for all workers to finish
    for p, _ in processes:
        p.join()

if __name__ == "__main__":
    # station_list = [s.replace('.nc', '') for s in os.listdir("/data/CARRA/extracted/list_pixels_minimal_grid/")]
    station_list = [s.replace('.nc', '') for s in os.listdir("./input/weather data/CARRA_at_AWS")]
    station_list.sort()

    standard_run_parallel(station_list)

    # for single station or debugging runs
    # for station in ['TAS_U']:
    # for station in station_list:
        # try:
        #     run_SEB_firn(station)
        # except Exception as e:
        #     print(station,e); traceback.print_exc()
