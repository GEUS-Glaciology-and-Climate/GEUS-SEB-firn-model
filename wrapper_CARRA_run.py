# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""

import xarray as xr
import cfgrib
import matplotlib.pyplot as plt
import os 
import geopandas as gpd
import pandas as pd
import rioxarray
from rasterio.crs import CRS
import warnings
warnings.filterwarnings("ignore", category=FutureWarning) 

ds = xr.open_dataset('model_mask.tif').squeeze('band').band_data.rename('model_mask').drop(['band', 'spatial_ref'])
df_list_pixels = ds.to_dataframe().reset_index()
df_list_pixels = df_list_pixels.loc[df_list_pixels.model_mask == 1, :]
plt.close('all')

# set directory
# os.chdir(r'Y:\CARRA')
print(os.getcwd())

# variable that should be read
folder_list = [
           '2m_relative_humidity', 
           '2m_specific_humidity', 
           '2m_temperature', 
           '10m_wind_speed', 
           'albedo',
           'skin_temperature', 
           'surface_latent_heat_flux',
           'surface_pressure', 
           'surface_sensible_heat_flux',
           'surface_solar_radiation_downwards'
           'surface_thermal_radiation_downwards',
           'time_integral_of_rain_flux', 
           'total_precipitation', 
          ]

for x,y in zip(df_list_pixels.x, df_list_pixels.y):
    break
    #%%
    for folder in folder_list:
        file_list = os.listdir('Y:/CARRA/'+folder)
        path_list = ['Y:/CARRA/'+folder+'/'+f for f in file_list]
        break
    # %% 
        for p in path_list:
            ds_CARRA = xr.open_dataset(p)
            break
        # %% 
    
        CARRA_nearest_point = []
        for f in file_list:
            if not f.endswith('.grib'):
                continue
            print(f)
            tmp = xr.open_dataset(folder+'/'+f, engine='cfgrib')
            tmp['longitude'] = xr.where(tmp['longitude']>180, 
                                        tmp['longitude']-360,
                                        tmp['longitude'])
            # tmp = tmp.isel(step=1)- tmp.isel(step=0)
            tmp['x'] = frac.x
            tmp['y'] = frac.y
            tmp = tmp.rio.write_crs(crs)
     
            # extracting values at the nearest point
            # print(type(tmp.indexes),tmp.indexes)
            tmp_nearest_point = tmp.interp(x = xr.DataArray(df.geometry.x.to_list(), dims='station'),
                                         y = xr.DataArray(df.geometry.y.to_list(), dims='station'), 
                                        method='nearest')
            
        
            if not CARRA_nearest_point:
                CARRA_nearest_point = tmp_nearest_point
            else:
                CARRA_nearest_point = xr.concat((CARRA_nearest_point, tmp_nearest_point), 
                                                dim='time')
                
        CARRA_nearest_point = CARRA_nearest_point.sortby('time')
        CARRA_nearest_point['latitude'] = ('station', df.lat)
        CARRA_nearest_point['longitude'] = ('station', df.lon)
        CARRA_nearest_point['altitude'] = ('station', df.alt)
        
        CARRA_nearest_point['altitude_mod'] = alt_carra.interp(
            x = xr.DataArray(df.geometry.x.to_list(), dims='station'),
            y = xr.DataArray(df.geometry.x.to_list(), dims='station'), 
            method='nearest').orog
        CARRA_nearest_point['name'] = ('station', df.stid)
        CARRA_nearest_point.to_netcdf('bav/PROMICE/CARRA_'+list(tmp.keys())[0]+'_nearest_PROMICE.nc')
        
        # CARRA_nearest_point.to_netcdf('bav/PROMICE/fbri_test/KAN_U_M/CARRA_'+list(tmp.keys())[0]+'_nearest_PROMICE.nc')
