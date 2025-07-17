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
import xarray as xr
import matplotlib.pyplot as plt


# var1 = 'rhofirn'
# var2 = 'rho_firn_only'
# var1 = 'T_ice'
# var2 = 'T_ice'
# var1 = 'snowc'
# var2 = 'snowc'
# var1 = 'snic'
# var2 = 'snic'
# var1 = 'slwc'
# var2 = 'slwc'
var1 = "density_bulk"
var2 = "rho"
filename_python = "Output/IMAU_aws4_51_layers_7/IMAU_aws4_" + var1 + ".nc"
filename_matlab = (
    "../GEUS model/Output/IMAU_aws4_0_IWC_CL_50_layers_4/IMAU_aws4_" + var2 + ".nc"
)

var_p = xr.open_dataset(filename_python)
var_m = xr.open_dataset(filename_matlab).transpose()

fig, ax = plt.subplots(1, 1)
cb = ax.pcolor((var_p[var1].values - var_m[var2].values), cmap="gist_ncar")
ax.set_title("diff")
fig.colorbar(cb, ax=ax)
print(np.mean((var_p[var1].values - var_m[var2].values)))
# i=-1
# fig, ax = plt.subplots(1,1)
# for i in range(20):
#     # i=i+1
#     var_m.rho_firn_only.sel(level = i+1).plot(ax=ax)
#     var_p.rhofirn.sel(level = i).shift(time=0).plot(ax=ax, linestyle='--')
# k=0
# fig, ax = plt.subplots(1,1)
# for i in range(20):
#     # i=i+1
#     var_m.Depth.sel(level = i+1).plot(ax=ax)
#     var_p.depth.sel(level = i).shift(time=0).plot(ax=ax, linestyle='--')
# k=0
# %%

print("python")
print(var_p.isel(time=k).to_dataframe()[[var1, "depth"]].head(5))
print("matlab")
print(var_m.isel(time=k).to_dataframe()[[var2, "Depth"]].head(5))
k = k + 1
print("diff")
print((var_p.rhofirn.isel(time=k).values - var_m.rho_firn_only.isel(time=k).values)[:5])
# %%
fig, ax = plt.subplots(1, 3)
var_p.depth.plot(ax=ax[0])
ax[0].set_title("python")

var_m.Depth.plot(ax=ax[1])
ax[0].set_title("matlab")

(var_p.depth - var_m.Depth).plot(ax=ax[2], cmap="seismic")
ax[2].set_title("diff")

fig, ax = plt.subplots(1, 1)
for i in range(20):
    var_m.Depth.sel(level=i + 1).plot(ax=ax)
    var_p.depth.sel(level=i).plot(ax=ax, linestyle="--")
