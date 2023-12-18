# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""

import lib.io as io
from lib.initialization import Struct
import numpy as np
import matplotlib.pyplot as plt
c= {
    'station': 'KAN_M',
    'surface_input_path': "./input/weather data/CARRA_at_AWS.nc",
    }
c = Struct(**c)

# loading CARRA
c.surface_input_driver = "CARRA" 
df_carra = io.load_surface_input_data(c)

# loading old file
c.surface_input_path = "./input/weather data/data_"+c.station+"_combined_hour.txt"
c.surface_input_driver = "AWS_old" 
df_old = io.load_surface_input_data(c)

# some filtering
df_carra.loc[df_carra.Albedo<0.1, 'Albedo'] = np.nan
df_old.loc[df_old.Albedo<0.1, 'Albedo'] = np.nan

#%%
def new_fig(var_list):
    fig, ax = plt.subplots(len(var_list), 1, sharex=True, figsize=(15, 10))
    if len(var_list) == 1: ax = [ax]
    plt.subplots_adjust(
        left=0.1, right=0.9, top=0.97, bottom=0.1, wspace=0.2, hspace=0.05
    )
    return fig, ax

# plt.close('all')

var_list = ['Snowfallmweq']

fig, ax = new_fig(var_list)
count = 0
count_fig = 0
for i, var in enumerate(var_list):
    var = str(var)           
    if "_origin" in var.lower(): continue
    if var not in df_old.columns: continue
    if 'Height' in var: continue
    print(var)
    if 'fall' in var:
        (df_carra[var].cumsum() - df_carra[var].cumsum().loc[df_old.index[0]]).plot(ax=ax[count], label='CARRA')
        df_old[var].cumsum().plot(ax=ax[count], label='AWS_old')
    else:
        df_carra[var].plot(ax=ax[count], label='CARRA', marker='.')
        df_old[var].plot(ax=ax[count], label='AWS_old', marker='.')
    ax[count].set_ylabel(var)
    ax[count].grid()
    # ax[count].set_xlim((df_old.index[0], df_old.index[-1]))
    
    # if var == "L":    #Changing the y-axis for L
    #     ax[count].set_ylim((-30000, 30000))


    count = count + 1

    if count == len(ax) &(var !=var_list[-1]):
        ax[0].set_title(c.station)
        count_fig = count_fig + 1
        fig, ax = new_fig(var_list)
        count = 0

    
