# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import numpy as np
import pandas as pd
import lib.io as io
from lib.initialization import ImportConst
from lib.seb_smb_model import HHsubsurf
import os
import time
import plot_output as po
import xarray as xr
import multiprocessing

def run_SEB_firn(station='DY2'):
    start_time = time.time()   
    SPIN_UP = False
    ds_carra = xr.open_dataset("./input/weather data/CARRA_at_AWS.nc")
   
    # if  os.path.isfile('./input/initial state/spin up/'+station+'_initial_T_ice.csv'):
    #     print('spin up already done for',station)
    #     return None

    print(station)
    # try:
    # importing standard values for constants
    c = ImportConst()
    c.output_path = 'C:/Users/bav/data_save/output firn model/'
    c.station = station        
    c.surface_input_path = "./input/weather data/CARRA_at_AWS.nc"
    c.surface_input_driver = "CARRA" 
    c.altitude= ds_carra.where(ds_carra.stid==c.station, drop=True).altitude.item()
    c.latitude= ds_carra.where(ds_carra.stid==c.station, drop=True).latitude.item()
    c.longitude= ds_carra.where(ds_carra.stid==c.station, drop=True).longitude.item()
    if c.altitude<1500:
        c.new_bottom_lay=30
    # assigning constants specific to this simulation
    c.snowthick_ini = 0.1
    c.z_max = 70
    c.num_lay = 100
    c.lim_new_lay = 0.05
    c.initial_state_folder_path = './input/initial state/spin up/'
    

        
    # assgning deep temperature
    ds_T10m = xr.open_dataset('../../Data/Firn temperature/output/T10m_prediction.nc')
    c.Tdeep = ds_T10m.sel(latitude = c.latitude,
                    longitude = c.longitude,
                    method='nearest').sel(time=slice('1980','1990')
                  ).T10m.mean().item() + 273.15
    # c.lim_new_lay = c.accum_AWS/c.new_lay_frac;
    
    # loading input data
    df_in = io.load_surface_input_data(c, resample=False)
    freq = pd.infer_freq(df_in.index)
    if freq=='H': freq = '1H'
    c.delta_time = pd.to_timedelta(freq).total_seconds()
    c.zdtime = c.delta_time 

    # Spin up option
    if SPIN_UP:
        df_in = pd.concat((df_in.loc['1991':'2001',:],
                           df_in.loc['1991':'2001',:],
                           df_in.loc['1991':'2001',:]), ignore_index=True)
        df_in.index=pd.to_datetime('1991-01-01T00:00:00') + pd.to_timedelta(df_in.index.astype(str).to_series() + 'H')
    
    for var in ['AirTemperature2C', 'Albedo', 'ShortwaveRadiationDownWm2', 'LongwaveRadiationDownWm2',
           'AirPressurehPa', 'WindSpeed2ms', 'RelativeHumidity2', 'LongwaveRadiationUpWm2', 
           'ShortwaveRadiationUpWm2', 'Snowfallmweq', 'Rainfallmweq', 'HeightTemperature2m', 
           'HeightHumidity2m', 'HeightWindSpeed2m']:
        assert ~df_in[var].isnull().any(), var+' has NaNs'
    print(station, c.Tdeep, 'start/end', df_in.index[0], df_in.index[-1])
    # DataFrame for the surface is created, indexed with time from df_aws
    df_out = pd.DataFrame()
    df_out["time"] = df_in.index
    df_out = df_out.set_index("time")
    
    # %% Running model
    print('reading inputs took %0.03f sec'%(time.time() -start_time))
    start_time = time.time()
    # The surface values are received 
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
    ) = HHsubsurf(df_in, c)
    
    print('\nHHsubsurf took %0.03f sec'%(time.time() -start_time))
    start_time = time.time()
    
    thickness_act = snowc * (c.rho_water / rhofirn) + snic * (c.rho_water / c.rho_ice)
    depth_act = np.cumsum(thickness_act, 0)
    density_bulk = (snowc + snic) / (snowc / rhofirn + snic / c.rho_ice)
    
    # %% Writing output
    if SPIN_UP:
        c.RunName = c.station + "_" + str(c.num_lay) + "_layers_SU_"+freq
    else:
        c.RunName = c.station + "_" + str(c.num_lay) + "_layers_"+freq
    i = 0
    succeeded = 0
    while succeeded == 0:
        try:
            os.mkdir(c.output_path + c.RunName)
            succeeded = 1
        except:
            if i == 0:
                c.RunName = c.RunName + "_" + str(i)
            else:
                c.RunName = c.RunName[: -len(str(i - 1))] + str(i)
            i = i + 1
    
    io.write_2d_netcdf(snic, 'snic', depth_act, df_in.index, c)
    io.write_2d_netcdf(snowc, 'snowc', depth_act, df_in.index, c)
    io.write_2d_netcdf(slwc, 'slwc', depth_act, df_in.index, c)
    io.write_2d_netcdf(rhofirn, 'rhofirn', depth_act, df_in.index, c)
    io.write_1d_netcdf(df_out, c)
    io.write_2d_netcdf(density_bulk, 'density_bulk', depth_act, df_in.index, c)
    io.write_2d_netcdf(T_ice, 'T_ice', depth_act, df_in.index, c)
    # io.write_2d_netcdf(rfrz, 'rfrz', depth_act, df_in.index, RunName)
    io.write_2d_netcdf(dgrain, 'dgrain', depth_act, df_in.index, c)
    io.write_2d_netcdf(compaction, 'compaction', depth_act, df_in.index, c)
    df_out['zrfrz_sum'] = zrfrz.sum(axis=0)
    
    c.write(c.output_path + "/" + c.RunName + "/constants.csv")
    print('writing output files took %0.03f sec'%(time.time() -start_time))
    start_time = time.time()
    
    # Plot output
    if not SPIN_UP:
        po.main(c.output_path, c.RunName)
        print('plotting took %0.03f sec'%(time.time() -start_time))
    start_time = time.time()
    return c
    # except Exception as e:
    #     print(station,'\n',e,'\n')
    #     return None
    
    
if __name__ == "__main__":
    # c = run_SEB_firn()
    run_SEB_firn('KAN_U')
    # ds_carra = xr.open_dataset("./input/weather data/CARRA_at_AWS.nc")
    # pool = multiprocessing.Pool(5)
    # station_list = ds_carra.stid.values[1:]
    # station_list = 'DYE2,NAE,NSE,NEM,HUM,TUN,CP1,SDM,SDL,KAN_U,NAU,Summit'.split(',')
    # out1, out2, out3 = zip(*pool.map(run_SEB_firn,station_list))
