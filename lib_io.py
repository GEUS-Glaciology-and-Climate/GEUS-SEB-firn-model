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

units = {'snowc': 'm water equivalent',
         'snic': 'm water equivalent',
         'slwc': 'm water equivalent',
         'rhofirn': 'kg m^-3',
         'density_bulk': 'kg m^-3',
         'T_ice': 'K',
         'compaction': 'm per time step',
         'rfrz': 'm water equivalent per time step',
         'dgrain': 'mm'}
long_name = {'snowc': 'Layer snow content',
             'snic': 'Layer ice content',
             'slwc': 'Layer liquid water content',
             'rhofirn': 'Density of snow only',
             'density_bulk': 'Bulk density',
             'T_ice': 'Subsurface temperature',
             'compaction': 'Layer compaction',
             'rfrz': 'Amount of water refrozen',
             'dgrain': 'Snow grain diameter'}


def write_2d_netcdf(data, name_var, depth_act, time, c):
    levels = np.arange(data.shape[0])
    time_days_since = pd.to_timedelta(time - np.datetime64('1900-01-01', 'ns')).total_seconds().values / 3600 / 24

    foo = xr.DataArray(data, coords=[levels, time_days_since], dims=["level",  "time"], name=name_var)
    foo.attrs["units"] = units[name_var]
    foo.attrs["long_name"] = long_name[name_var]
    foo.time.attrs["units"] = "days since 1900-01-01"
    foo.level.attrs["units"] = "index of layer (0=top)"

    depth = xr.DataArray(depth_act, coords=[levels, time_days_since],
                         dims=["level", "time"], name="depth")
    depth.attrs["units"] = "m below the surface"
    depth.attrs["long_name"] = "Depth of layer bottom"
    depth.time.attrs["units"] = "days since 1900-01-01"
    depth.level.attrs["units"] = "index of layer (0=top)"

    ds = xr.merge([foo, depth])
    ds.to_netcdf('Output/'+c.RunName+'/'+c.station+'_'+name_var+'.nc')
