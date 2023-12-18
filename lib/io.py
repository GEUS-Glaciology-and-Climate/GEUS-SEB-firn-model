# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import numpy as np
import xarray as xr
import pandas as pd
import datetime
import pytz

units = {
    "snowc": "m water equivalent",
    "snic": "m water equivalent",
    "slwc": "m water equivalent",
    "rhofirn": "kg m^-3",
    "density_bulk": "kg m^-3",
    "T_ice": "K",
    "compaction": "m per time step",
    "rfrz": "m water equivalent per time step",
    "dgrain": "mm",
}
long_name = {
    "snowc": "Layer snow content",
    "snic": "Layer ice content",
    "slwc": "Layer liquid water content",
    "rhofirn": "Density of snow only",
    "density_bulk": "Bulk density",
    "T_ice": "Subsurface temperature",
    "compaction": "Layer compaction",
    "rfrz": "Amount of water refrozen",
    "dgrain": "Snow grain diameter",
}


def load_promice_old(path_promice):
    """
    Loading PROMICE data for a given path into a DataFrame.
    + adding time index
    + calculating albedo
    + (optional) calculate RH with regard to water

    INTPUTS:
        path_promice: Path to the desired file containing PROMICE data [string]

    OUTPUTS:
        df: Dataframe containing PROMICE data for the desired settings [DataFrame]
    """

    df = pd.read_csv(path_promice, delim_whitespace=True)
    df["time"] = df.Year * np.nan

    try:
        df["time"] = [
            datetime.datetime(y, m, d, h).replace(tzinfo=pytz.UTC)
            for y, m, d, h in zip(
                df["Year"].values,
                df["MonthOfYear"].values,
                df["DayOfMonth"].values,
                df["HourOfDay(UTC)"].values,
            )
        ]
    except:
        df["time"] = pd.to_datetime(
            df["Year"] * 100000 + df["DayOfYear"] * 100 + df["HourOfDayUTC"],
            format="%Y%j%H",
        )

    df.set_index("time", inplace=True, drop=False)

    # set invalid values (-999) to nan
    df[df == -999.0] = np.nan
    # df['Albedo'] = df['ShortwaveRadiationUp(W/m2)'] / df['ShortwaveRadiationDown(W/m2)']
    df.loc[df["Albedo"] > 1, "Albedo"] = np.nan
    df.loc[df["Albedo"] < 0, "Albedo"] = np.nan
    # df['SnowHeight(m)'] = 2.6 - df['HeightSensorBoom(m)']
    # df['SurfaceHeight(m)'] = 1 - df['HeightStakes(m)']

    # df['RelativeHumidity_w'] = RH_ice2water(df['RelativeHumidity(%)'] ,
    #                                                    df['AirTemperature(C)'])
    
    df = df.loc[df.AirPressurehPa.first_valid_index():,:]
    df = df.loc[df.AirTemperature1C.first_valid_index():,:]

    df = df.set_index("time").resample("H").mean()
    df = df.interpolate()
    return df


def load_CARRA_data(*args):
    if len(args) == 1:
        c = args[0]
        surface_input_path = c.surface_input_path
        station = c.station
    elif len(args) == 2:
        surface_input_path, station = args
        
    print("- Reading data from CARRA reanalysis set -")
    aws_ds = xr.open_dataset(surface_input_path)
    
    df_carra = aws_ds.where(aws_ds.stid==station, drop=True).squeeze().to_dataframe()
    df_carra['HeightTemperature2m'] = 2
    df_carra['HeightHumidity2m'] = 2
    df_carra['HeightWindSpeed2m'] = 10

    # converting to a pandas dataframe and renaming some of the columns
    df_carra = df_carra.rename(columns={
                            't2m': 'AirTemperature2C', 
                            'r2': 'RelativeHumidity2', 
                            'si10': 'WindSpeed2ms', 
                            'sp': 'AirPressurehPa', 
                            'ssrd': 'ShortwaveRadiationDownWm2',
                            'ssru': 'ShortwaveRadiationUpWm2',
                            'strd': 'LongwaveRadiationDownWm2',
                            'stru': 'LongwaveRadiationUpWm2',
                            'sf': 'Snowfallmweq',
                            'rf': 'Rainfallmweq',
                            'al': 'Albedo',
                        })

    # Fill null values with 0
    df_carra['ShortwaveRadiationDownWm2'] = df_carra['ShortwaveRadiationDownWm2'].fillna(0)
    df_carra['ShortwaveRadiationUpWm2'] = df_carra['ShortwaveRadiationUpWm2'].fillna(0)
    df_carra = df_carra.resample('H').interpolate()
    df_carra['LongwaveRadiationDownWm2'] = df_carra['LongwaveRadiationDownWm2']+18  # bias adjustment
    # df_carra['AirTemperature2C'] = df_carra['AirTemperature2C']+1.69  # bias adjustment
    # df_carra['ShortwaveRadiationDownWm2'] = df_carra['ShortwaveRadiationDownWm2']*1.3  # bias adjustment
    df_carra['ShortwaveRadiationUpWm2'] = df_carra.ShortwaveRadiationDownWm2*df_carra.Albedo
    
    # calcualting snowfall and rainfall
    df_carra['Snowfallmweq'] = 0  # df_carra['Snowfallmweq'] /3
    df_carra['Rainfallmweq'] = 0  # df_carra['Rainfallmweq'] /3
    cut_off_temp = 0
    df_carra.loc[df_carra.AirTemperature2C < cut_off_temp,
                 'Snowfallmweq'] = df_carra.loc[
                     df_carra.AirTemperature2C < cut_off_temp,'tp'] /3/ 1000/1.5
    df_carra.loc[df_carra.AirTemperature2C >= cut_off_temp,
                 'Rainfallmweq'] = df_carra.loc[
                     df_carra.AirTemperature2C >= cut_off_temp,'tp'] /3/ 1000/1.5
    return df_carra


def load_surface_input_data(c):
    if c.surface_input_driver  == 'AWS_old':
        return load_promice_old(c.surface_input_path)
    # if c.surface_input_driver  == 'AWS':
    #     return load_promice(c.surface_input_path)
    if c.surface_input_driver  == 'CARRA':
        return load_CARRA_data(c)
    
    print('Driver', c.surface_input_driver , 'not recognized')
    return None
    
    
def write_2d_netcdf(data, name_var, depth_act, time, c):
    levels = np.arange(data.shape[0])
    time_days_since = (
        pd.to_timedelta(time - np.datetime64("1900-01-01", "ns")).total_seconds().values
        / 3600
        / 24
    )

    foo = xr.DataArray(
        data, coords=[levels, time_days_since], dims=["level", "time"], name=name_var
    )
    foo.attrs["units"] = units[name_var]
    foo.attrs["long_name"] = long_name[name_var]
    foo.time.attrs["units"] = "days since 1900-01-01"
    foo.level.attrs["units"] = "index of layer (0=top)"

    depth = xr.DataArray(
        depth_act,
        coords=[levels, time_days_since],
        dims=["level", "time"],
        name="depth",
    )
    depth.attrs["units"] = "m below the surface"
    depth.attrs["long_name"] = "Depth of layer bottom"
    depth.time.attrs["units"] = "days since 1900-01-01"
    depth.level.attrs["units"] = "index of layer (0=top)"
    
    float_encoding = {"dtype": "float32", "zlib": True,"complevel": 9}
    int_encoding = {"dtype": "int32", "zlib": True,"complevel": 9}

    ds = xr.merge([foo, depth])
    ds.to_netcdf(
        c.output_path + "/" + c.RunName + "/" + c.station + "_" + name_var + ".nc",
        encoding = {'depth': float_encoding,
                    name_var: float_encoding,
                    'level':int_encoding}
    )


def write_1d_netcdf(data, c, var_list=None, time=None, name_file="surface"):
    var_info = pd.read_csv("Input/Constants/output_variables_info.csv").set_index(
        "short_name"
    )
    if not time:
        time = data.index
    if not var_list:
        var_list = data.columns
    time_days_since = (
        pd.to_timedelta(time - np.datetime64("1900-01-01", "ns")).total_seconds().values
        / 3600
        / 24
    )
    float_encoding = {"dtype": "float32", "zlib": True,"complevel": 9}

    for name_var in var_list:
        foo = xr.DataArray(
            data[name_var].values,
            coords=[time_days_since],
            dims=["time"],
            name=name_var,
        )
        foo.attrs["units"] = var_info.loc[name_var, "units"]
        foo.attrs["long_name"] = var_info.loc[name_var, "long_name"]
        foo.encoding = float_encoding
        foo.time.attrs["units"] = "days since 1900-01-01"
        if name_var == var_list[0]:
            ds = foo
        else:
            ds = xr.merge([ds, foo])
    ds.to_netcdf(
        c.output_path + "/" + c.RunName + "/" + c.station + "_" + name_file + ".nc"
    )
