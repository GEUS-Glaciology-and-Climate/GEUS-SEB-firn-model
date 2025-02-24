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
import os
import lib.io as io
from multiprocessing import Pool

def name_alias(stid):
    rename = {'South Dome':'SDM', 'Saddle':'SDL', 'NASA-U': 'NAU',
                'NASA-E': 'NAE', 'NEEM': 'NEM', 'EastGRIP': 'EGP',
                'DYE-2': 'DY2', 'Tunu-N':'TUN', 'CEN1':'CEN', 'CEN2':'CEN',
                'JAR1':'JAR', 'NASA-SE':'NSE','GITS':'CEN','Humboldt':'HUM',
                # ['Summit', 'DMI'],
                # ['Summit', 'NOAA']
}
    if stid in rename.keys():
        return rename[stid]
    else:
        return stid
# output_path= 'C:/Users/bav/data_save/output firn model/spin up 3H/'
output_path = './output/new/'
run_name = '10450_100_layers_3h'
#%%
def main(output_path, run_name):
    # %% Loading data
    print(run_name)
    tmp =pd.read_csv(output_path+'/'+ run_name+'/constants.csv', dtype={'key':str})
    # converting all numerical fields to numeric, except station
    tmp['value_num'] = pd.to_numeric(tmp.value, errors='coerce')
    msk = (tmp.value_num.notnull() & (tmp.key!='station'))
    tmp.loc[msk,'value'] = tmp.loc[msk,'value_num']
    tmp = tmp.set_index('key')[['value']]
    # making it a structure
    c = Struct(**tmp.to_dict()['value'] )
    c.RunName=run_name

    if c.surface_input_driver=='CARRA' and c.zdtime == 3600:
        print('resample')
        resample=True
    else:
        resample=False

    df_in, c = io.load_surface_input_data(c, resample=resample)
    if output_path != c.output_path:
        print('Warning: Output has been moved from',c.output_path,'to',output_path)
        c.output_path = output_path

    #  loading surface variables
    try:
        df_out = xr.open_dataset(c.output_path+run_name+'/'+c.station+'_surface.nc').to_dataframe()
        df_in = df_in.loc[df_out.index[0]:df_out.index[-1],:]
    except Exception as e:
        print(c.RunName, e)

    # %% plotting surface variables
    lpl.plot_summary(df_out, c, 'SEB_output')

    # %% plotting subsurface variables
    for var in ['compaction','T_ice','density_bulk','slwc','dgrain']:
        if len(df_in) <300:
            ylim =   [10]
        else:
            ylim = []
        lpl.plot_var(c.station, c.output_path, c.RunName, var, ylim=ylim, zero_surf=False)

    # if c.station in ['DY2', 'KAN_U','CP1']:
        # lpl.plot_var(c.station, c.output_path, c.RunName, 'slwc',
                     # zero_surf=True, ylim=(8,0), year = (2012, 2024))


    # %% Surface height evaluation
    # extracting surface height

    path_aws_l4 = '../thredds-data/level_3_sites/csv/hour/'
    # path_aws_l4 = 'C:/Users/bav/GitHub/PROMICE data/thredds/level_3_sites/hour/'
    if os.path.isfile(path_aws_l4+name_alias(c.station)+'_hour.csv'):
        df_obs = pd.read_csv(path_aws_l4+name_alias(c.station)+'_hour.csv')
        obs_avail = True
    else:
        path_aws_l4 = '../GC-Net-Level-1-data-processing/L1/hour/'
        if os.path.isfile(path_aws_l4+c.station.replace(' ','')+'.csv'):
            import nead
            df_obs = nead.read(path_aws_l4+c.station.replace(' ','')+'_daily.csv').to_dataframe()
            df_obs = df_obs.rename(columns={'timestamp':'time',
                                            'HS_combined':'z_surf_combined',
                                            'T10m': 't_i_10m',
                                            'LHF':'dlhf_u',
                                            'SHF': 'dshf_u',
                                            'OLWR':'LRout',
                                            'Tsurf':'t_surf',
                                            })
            obs_avail = True
        else:
            print(c.RunName, ': no weather observation was found')
            obs_avail = False

            # return []
        # else:
        #     tmp = pd.DataFrame()
    # if len(df_obs)>0:
    #     df_obs.time= pd.to_datetime(df_obs.time).dt.tz_convert(None)
    #     if len(tmp)>0:
    if obs_avail:
        df_obs.time= pd.to_datetime(df_obs.time)
        df_obs = df_obs.set_index('time')
        df_obs = df_obs.resample(pd.infer_freq(df_out.index)).mean()

        fig = plt.figure()
        tmp = (df_obs.z_surf_combined -df_out.surface_height).mean()
        plt.plot(df_obs.index, df_obs.z_surf_combined-tmp,
                  marker='.',ls='None', label='AWS')
        plt.plot(df_out.index, df_out.surface_height,color='tab:red',
                 label='model')
        plt.legend()
        plt.ylabel('Surface height (m)')
        plt.title(c.station)
        fig.savefig(c.output_path+c.RunName+'/'+c.station+'_surface_height.png', dpi=120)

    # %% calculating modelled t_i_10m

    filename = c.output_path + run_name + "/" + c.station + "_T_ice.nc"
    df = (xr.open_dataset(filename).to_dataframe().unstack('level'))
    df.columns = df.columns.map('{0[0]}_{0[1]}'.format)
    # df_10m = interpolate_temperature(
    #     df.index, df[[v for v in df.columns if 'depth' in v]].values,
    #     df[[v for v in df.columns if 'T_ice' in v]].values-273.15,
    # )
    # df_out['t_i_10m'] = df_10m.temperatureObserved.values
    from lib.plot import interpolate_temperature_fast
    df_out['t_i_10m'] = interpolate_temperature_fast(
        df.index, df[[v for v in df.columns if 'depth' in v]].values,
        df[[v for v in df.columns if 'T_ice' in v]].values-273.15,
    )


    # %% plotting ['t_surf','LRout','LHF','SHF','t_i_10m']
    df_out['t_surf']  =  ((df_out.LRout_mdl - (1 -  c.em) * df_out.LRin) / c.em / 5.67e-8)**0.25 - 273.15
    df_out['LRout'] = df_out.LRout_mdl

    if obs_avail:
        for var1, var2 in zip(['LHF','SHF','LRout'],
                                ['dlhf_u', 'dshf_u', 'ulr']):
            if var2 in df_obs.columns:
                df_obs[var1] = df_obs[var2]
            else:
                df_obs[var1] = np.nan
        lpl.plot_observed_vars(df_obs, df_out, c, var_list = ['t_surf','LRout','LHF','SHF','t_i_10m'])

    lpl.plot_smb_components(df_out, c)
    lpl.evaluate_temperature_sumup(df_out, c)
    # lpl.evaluate_temperature_scatter(df_out, c, year = None)
    lpl.evaluate_density_sumup(c)
    lpl.evaluate_smb_sumup(df_out, c)
    lpl.evaluate_accumulation_snowfox(df_in, c)
    lpl.plot_var_start_end(c, 'T_ice')
    lpl.plot_var_start_end(c, 'density_bulk')
    # lpl.plot_movie(c.station, c.output_path, c.RunName, 'T_ice')
    # lpl.plot_movie(c.station, c.output_path, c.RunName, 'density_bulk')
    lpl.evaluate_compaction(c)

    plt.close('all')
    # try:
    #     lpl.find_summer_surface_depths(c)
    # except Exception as e:
    #     print(e)
    #     pass

# %%
if __name__ == "__main__":
    # for run_name in os.listdir('output/new/'):
    #     main(output_path=output_path, run_name=run_name)
    main(output_path=output_path, run_name=run_name)

    # run_name_list = os.listdir('output/new/')

    # def main_wrapper(run_name):
    #     main(output_path=output_path, run_name=run_name)

    # with Pool(7, maxtasksperchild=1) as pool:
    #     pool.map(main_wrapper, run_name_list, chunksize=1)
