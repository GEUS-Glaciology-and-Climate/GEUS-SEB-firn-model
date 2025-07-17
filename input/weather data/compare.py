# -*- coding: utf-8 -*-
"""
Created on %(date)s
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""

import xarray as xr

# ds1 = xr.open_dataset('FA-13.nc')
# ds2 = xr.open_dataset('CARRA_at_AWS.nc').isel(station=-3)
# %%
ds1 = xr.open_dataset('../../output/new/FA-13_100_layers_3H/FA-13_surface.nc')
ds2 = xr.open_dataset('../../output/FA-13_100_layers_3H/FA-13_surface.nc')

# %%
import matplotlib.pyplot as plt

for var in ds1.data_vars:
    plt.figure()
    # try:
    ds1[var].cumsum().plot(label='new', lw=2)
    ds2[var].cumsum().plot(label='old')
    plt.legend()
    # except:
    #     pass
