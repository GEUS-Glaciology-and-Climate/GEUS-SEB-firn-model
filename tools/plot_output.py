# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import netCDF4 as nc
import matplotlib.dates as mdates
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

import lib_plot as lp

site = "NASA-SE"
run_name = "NASA-SE_100_layers"

# %% Plotting 2D output
plt.close("all")
for var in ["density_bulk", "rhofirn", "T_ice", "slwc"]:
    print("Plotting " + var)
    lp.plot_var(site, run_name, var, zero_surf=False)


# %% Plotting layer distribution
var_name = "rhofirn"
filename = "Output/" + run_name + "/" + site + "_" + var_name + ".nc"
var = xr.open_dataset(filename)

fig, ax = plt.subplots(1, 1, figsize=(15, 30))
for i in range(len(var.level)):
    var.depth.sel(level=i).plot(ax=ax)
ax.invert_yaxis()
fig.savefig("Output/" + run_name + "/" + site + "_layers.png")

# %% Plotting compaction
