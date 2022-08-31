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


def load_promice(path_promice):
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

    if "hour" in path_promice:
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

    elif "day" in path_promice:
        df["time"] = [
            datetime.datetime(y, m, d).replace(tzinfo=pytz.UTC)
            for y, m, d in zip(
                df["Year"].values, df["MonthOfYear"].values, df["DayOfMonth"].values
            )
        ]
    elif "month" in path_promice:
        df["time"] = [
            datetime.datetime(y, m).replace(tzinfo=pytz.UTC)
            for y, m in zip(df["Year"].values, df["MonthOfYear"].values)
        ]

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

    return df


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

    ds = xr.merge([foo, depth])
    ds.to_netcdf(
        c.OutputFolder + "/" + c.RunName + "/" + c.station + "_" + name_var + ".nc"
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

    for name_var in var_list:
        foo = xr.DataArray(
            data[name_var].values,
            coords=[time_days_since],
            dims=["time"],
            name=name_var,
        )
        foo.attrs["units"] = var_info.loc[name_var, "units"]
        foo.attrs["long_name"] = var_info.loc[name_var, "long_name"]
        foo.time.attrs["units"] = "days since 1900-01-01"
        if name_var == var_list[0]:
            ds = foo
        else:
            ds = xr.merge([ds, foo])
    ds.to_netcdf(
        c.OutputFolder + "/" + c.RunName + "/" + c.station + "_" + name_file + ".nc"
    )
