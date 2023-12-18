# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import xarray as xr
ds = xr.open_zarr("https://snow.univ-grenoble-alpes.fr/opendata/SMOS-Greenland-25km.zarr")

compression = {'zlib': True, 'complevel': 4}

ds.to_netcdf("SMOS-Greenland-25km.nc", encoding={var: compression for var in ds.variables}, format="NETCDF4")

