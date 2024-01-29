# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import lib.plot as lpl
from lib.initialization import Struct
import pandas as pd
import lib.io as io
output_path= 'C:/Users/bav/data_save/output firn model/spin up/'
run_name = 'CEN2_100_layers_SU'
#%%
def main(output_path, run_name):
    # %% Loading data
    print(run_name)
    tmp =pd.read_csv(output_path+'/'+ run_name+'/constants.csv', dtype={'key':str})
    tmp['value_num'] = pd.to_numeric(tmp.value, errors='coerce')
    tmp.loc[tmp.value_num.notnull(),'value'] = tmp.loc[tmp.value_num.notnull(),'value_num']
    tmp = tmp.set_index('key')[['value']]
    c = Struct(**tmp.to_dict()['value'] )
    c.RunName=run_name
    # df_in = io.load_surface_input_data(c)
    if output_path != c.output_path:
        print('Warning: Output has been moved from',c.output_path,'to',output_path)
        c.output_path = output_path

    # df_in = pd.concat((df_in.loc['1991':'2001',:],
    #                    df_in.loc['1991':'2001',:],
    #                    df_in.loc['1991':'2001',:]), ignore_index=True)
    # df_in.index=pd.to_datetime('1991-01-01T00:00:00') + pd.to_timedelta(df_in.index.astype(str).to_series() + 'H')
    # #  loading and plotting surface variables
    # try:
    #     df_out = xr.open_dataset(c.output_path+run_name+'/'+c.station+'_surface.nc').to_dataframe()
    #     df_in = df_in.loc[df_out.index[0]:df_out.index[-1],:]
    #     df_in.AirTemperature2C.resample('3M').mean().plot()
    # except Exception as e:
    #     print(e)

    # %% Start/end plots
    lpl.plot_var_start_end(c, 'T_ice',to_file=True)
    lpl.plot_var_start_end(c, 'density_bulk',to_file=True)   
    lpl.plot_var_start_end(c, 'rhofirn',to_file=True)   
    lpl.plot_var_start_end(c, 'dgrain',to_file=True)   
    lpl.plot_var_start_end(c, 'snic',to_file=True)
    try:
        lpl.plot_var_start_end(c, 'snowc',to_file=True)   
    except:
        pass

# %%
import os    
if __name__ == "__main__":
    for run_name in os.listdir('C:/Users/bav/data_save/output firn model/spin up/'):
        main(output_path=output_path, run_name=run_name)
