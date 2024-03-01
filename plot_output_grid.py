# -*- coding: utf-8 -*-
"""
Created on %(date)s
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
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import matplotlib

output_folder = 'C:/Users/bav/data_save/output firn model/'
list_folder = [folder for folder in os.listdir(output_folder) if 'grid' in folder]

var = 'surface'

list_path = []
pixel = 120096
for folder in list_folder:
    year = folder.split('_')[1]
    month = folder.split('_')[2]
    list_path.append(
        output_folder+folder+'/pixel_'+str(pixel) +'_' + year + '_' + month +'_100_layers_3H' \
            +'/pixel_'+str(pixel)+'_'+year+'_'+month+'_'+var+'.nc')
    
ds_out = xr.open_mfdataset(list_path)

plt.figure()
ds_out.LRout_mdl.plot()
