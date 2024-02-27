# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""

import xarray as xr
import rioxarray 
from scipy.ndimage import label
import numpy as np
import matplotlib.pyplot as plt

ds = xr.open_dataset('input/weather data/fractions.west.nc')

ds_fraction = ds.fraction_permanent_snow #.rio.to_raster('CARRA_permanent_snow.tif')
ds_elevation = ds.orography 

ds_fraction = (ds_fraction>0)*1

ds_label = ds_fraction.copy()
ds_label.data = label(ds_fraction.data)[0]

GrIS_mask = (ds_label==1)
ds_elevation['i'] = ds_elevation*0
ds_elevation['j'] = ds_elevation*0
ds_elevation.i.data, ds_elevation.j.data = np.meshgrid(np.arange(len(ds_elevation.x)), np.arange(len(ds_elevation.y)))
step = 8
msk = (ds_elevation<1200) | (
    (ds_elevation<2000) & (
    ((ds_elevation.i % 2 ==0) & (ds_elevation.j % 2 ==0)) | \
    ((ds_elevation.i % 2 ==1) & (ds_elevation.j % 2 ==1))  ) ) | (
        (ds_elevation>2000) & (
        ((ds_elevation.i % step ==0) & (ds_elevation.j % step ==0)) | \
        ((ds_elevation.i % step ==step/2) & (ds_elevation.j % step ==step/2))  ) )
proc_mask = GrIS_mask.where(msk, other=0)
plt.figure()
proc_mask.plot()

print(proc_mask.sum().item())
print(proc_mask.sum().item() * 200 /60 /60 /24/30)

proc_mask.rename('mask').rio.to_raster('model_mask.tif')
