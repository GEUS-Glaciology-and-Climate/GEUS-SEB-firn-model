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
encoding = {
    "compaction": {"dtype": "int32", "scale_factor": 1e-8, "zlib": True, "complevel": 9},
    "dgrain": {"dtype": "int16", "scale_factor": 0.001, "zlib": True, "complevel": 9},
    "rhofirn": {"dtype": "int16", "scale_factor": 0.001, "zlib": True, "complevel": 9},
    "T_ice": {"dtype": "int16", "scale_factor": 0.1, "zlib": True, "complevel": 9},
    "level": {"dtype": "int32", "zlib": True,"complevel": 9},
    "time": {"zlib": True,"complevel": 9},
    "**": {"dtype": "float32", "zlib": True, "complevel": 9},  # Default for all other variables
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
        df["time"] = [datetime.datetime(y, m, d, h).replace(tzinfo=pytz.UTC)
            for y, m, d, h in zip(df["Year"].values,
                df["MonthOfYear"].values, df["DayOfMonth"].values,
                df["HourOfDay(UTC)"].values)]
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


def load_CARRA_data(*args, resample=True):
    if len(args) == 1:
        c = args[0]
        surface_input_path = c.surface_input_path
        station = c.station
    elif len(args) == 2:
        surface_input_path, station = args

    # print("- Reading data from CARRA reanalysis set -", surface_input_path)
    with xr.open_dataset(surface_input_path) as ds:
        aws_ds = ds.where(ds.stid==c.station, drop=True).load()

    if 'altitude' in aws_ds.data_vars:
        c.altitude= aws_ds.altitude.item()
    else:
        c.altitude= aws_ds.altitude_mod.item()
    c.latitude= aws_ds.latitude.item()
    c.longitude= aws_ds.longitude.item()
    if c.longitude>180:
        c.longitude = c.longitude-360

    df_carra = aws_ds.squeeze().to_dataframe()
    df_carra['HeightTemperaturem'] = 2
    df_carra['HeightHumiditym'] = 2
    df_carra['HeightWindSpeedm'] = 10

    # converting to a pandas dataframe and renaming some of the columns
    df_carra = df_carra.rename(columns={
                            't2m': 'AirTemperatureC',
                            'r2': 'RelativeHumidity',
                            'sh2': 'SpecificHumiditykgkg',
                            'si10': 'WindSpeedms',
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

    if resample:
        df_carra = df_carra.resample('H').interpolate()
    else:
        df_carra = df_carra.apply(pd.to_numeric, errors='coerce').resample(pd.infer_freq(df_carra.index)).interpolate()
    df_carra['LongwaveRadiationDownWm2'] = df_carra['LongwaveRadiationDownWm2']+18  # bias adjustment
    # df_carra['AirTemperature2C'] = df_carra['AirTemperature2C']+1.69  # bias adjustment
    # df_carra['ShortwaveRadiationDownWm2'] = df_carra['ShortwaveRadiationDownWm2']*1.3  # bias adjustment
    df_carra['ShortwaveRadiationUpWm2'] = df_carra.ShortwaveRadiationDownWm2*df_carra.Albedo

    # calcualting snowfall and rainfall
    df_carra['Snowfallmweq'] = 0.
    df_carra['Rainfallmweq'] = 0.

    T_min = 0.5  # Temperature below which all precipitation is snow
    T_max = 2.5  # Temperature above which all precipitation is rain

    # Compute snowfall fraction after Hock and Holmgren 2005
    snow_fraction = np.clip((T_max - df_carra.AirTemperatureC) / (T_max - T_min), 0, 1)

    # Assign snowfall and rainfall
    df_carra['Snowfallmweq'] = (df_carra['tp'] / 1000.) * snow_fraction
    df_carra['Rainfallmweq'] = (df_carra['tp'] / 1000.) * (1 - snow_fraction)

    if (df_carra.index[1] - df_carra.index[0]) == pd.Timedelta('1 hours'):
        df_carra['Snowfallmweq'] = df_carra['Snowfallmweq'] / 3
        df_carra['Rainfallmweq'] = df_carra['Rainfallmweq'] / 3

    df_carra['Rainfallmweq'] = df_carra['Rainfallmweq'].fillna(0)
    df_carra['Snowfallmweq'] = df_carra['Snowfallmweq'].fillna(0)
    return df_carra

def load_surface_input_data(c, resample=True):
    df_surf = None
    if c.surface_input_driver  == 'AWS_old':
        df_surf = load_promice_old(c.surface_input_path)
    # if c.surface_input_driver  == 'AWS':
    #     return load_promice(c.surface_input_path)
    if c.surface_input_driver  == 'CARRA':
        df_surf = load_CARRA_data(c, resample=resample)
    if c.surface_input_driver  == 'CARRA_grid':
        df_surf = load_CARRA_grid(c)

    # Spin up option
    if c.spin_up:
        df_surf = pd.concat((df_surf.loc['1991':'2001',:],
                           df_surf.loc['1991':'2001',:],
                           df_surf.loc['1991':'2001',:],
                           df_surf.loc['1991':'2001',:],
                           df_surf.loc['1991':'2001',:]), ignore_index=True)
        if c.freq == 'h':
            df_surf.index=pd.to_datetime('1991-01-01T00:00:00') + pd.to_timedelta(df_surf.index.astype(str).to_series() + c.freq)
        elif c.freq == '3h':
            df_surf.index=pd.to_datetime('1991-01-01T00:00:00') + pd.to_timedelta((df_surf.index.to_series()*3).astype(str) + 'h')

    for var in ['AirTemperatureC', 'ShortwaveRadiationDownWm2', 'LongwaveRadiationDownWm2',
           'AirPressurehPa', 'WindSpeedms', 'SpecificHumiditykgkg',
           'ShortwaveRadiationUpWm2', 'Snowfallmweq', 'Rainfallmweq', 'HeightTemperaturem',
           'HeightHumiditym', 'HeightWindSpeedm']:
         if df_surf[var].isnull().any():
             print('!!!!!!!!!!!!!!!!!!')
             print(var, 'at',c.station, 'has NaNs')
             print('!!!!!!!!!!!!!!!!!!')
             return

    if df_surf is None:
        print('Driver', c.surface_input_driver , 'not recognized')
    return df_surf, c


def write_2d_netcdf(data, name_var, depth_act, time, c):
    levels = np.arange(data.shape[0])

    # Convert time to CF-compliant format (datetime64[ns] for NetCDF)
    time_cf = xr.DataArray(pd.to_datetime(time), dims="time", name="time")

    foo = xr.DataArray(
        data,
        coords={"level": levels, "time": time_cf},
        dims=["level", "time"],
        name=name_var,
    )
    foo.attrs["units"] = units[name_var]
    foo.attrs["long_name"] = long_name[name_var]

    depth = xr.DataArray(
        depth_act,
        coords={"level": levels, "time": time_cf},
        dims=["level", "time"],
        name="depth",
    )
    depth.attrs["units"] = "m below the surface"
    depth.attrs["long_name"] = "Depth of layer bottom"

    ds = xr.Dataset({name_var: foo, "depth": depth})

    # ds["time"].attrs["units"] = "days since 1900-01-01"
    # ds["time"].attrs["calendar"] = "standard"  # Change to "proleptic_gregorian" if needed
    ds["level"].attrs["units"] = "index of layer (0=top)"

    ds.attrs.update({
        "title": f"Simulated {long_name[name_var].lower()} from the GEUS SEB-firn model",
        "contact": "bav@geus.dk",
        "production_date": datetime.date.today().isoformat(),
        "run_name": c.RunName,
        "latitude": c.latitude,
        "longitude": c.longitude,
        "altitude": c.altitude,
        "station": c.station,
    })

    for attr in ["pixel", "year", "month"]:
        if hasattr(c, attr):
            ds.attrs[attr] = getattr(c, attr)

    output_path = f"{c.output_path}/{c.RunName}/{c.station}_{name_var}.nc"

    filtered_encoding = {
        var: encoding.get(var, encoding.get("**", {})) for var in ds.variables
    }

    ds.to_netcdf(output_path, encoding=filtered_encoding)
    ds.close()


def write_1d_netcdf(data, c, var_list=None, time=None, name_file="surface"):
    var_info = pd.read_csv("input/constants/output_variables_info.csv").set_index("short_name")

    if time is None:
        time = data.index
    if var_list is None:
        var_list = data.columns

    # Convert time to CF-compliant format (datetime64[ns] for NetCDF)
    time_cf = xr.DataArray(pd.to_datetime(time), dims="time", name="time")

    ds = xr.Dataset(coords={"time": time_cf})

    for name_var in var_list:
        foo = xr.DataArray(
            data[name_var].values,
            coords={"time": ds["time"]},
            dims=["time"],
            name=name_var,
        )
        foo.attrs["units"] = var_info.loc[name_var, "units"]
        foo.attrs["long_name"] = var_info.loc[name_var, "long_name"]
        ds[name_var] = foo

    # ds["time"].attrs["units"] = "days since 1900-01-01"
    # ds["time"].attrs["calendar"] = "standard"  # Use "proleptic_gregorian" if required

    ds.attrs.update({
        "title": "Simulated surface variables from the GEUS SEB-firn model",
        "contact": "bav@geus.dk",
        "production_date": datetime.date.today().isoformat(),
        "run_name": c.RunName,
        "latitude": c.latitude,
        "longitude": c.longitude,
        "altitude": c.altitude,
        "station": c.station,
    })

    for attr in ["pixel", "year", "month"]:
        if hasattr(c, attr):
            ds.attrs[attr] = getattr(c, attr)

    output_path = f"{c.output_path}/{c.RunName}/{c.station}_{name_file}.nc"

    filtered_encoding = {
        var: encoding.get(var,
                          encoding.get("**", {})) for var in ds.variables
    }

    ds.to_netcdf(output_path, encoding=filtered_encoding)
    ds.close()
