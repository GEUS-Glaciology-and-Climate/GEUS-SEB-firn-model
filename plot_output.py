# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import xarray as xr
import lib.plot as lpl
from lib.initialization import Struct
import pandas as pd
import os
import lib.io as io

output_path= 'C:/Users/bav/data_save/output firn model/'
run_name = 'pixel_121162_100_layers_3H'
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
    df_in, c = io.load_surface_input_data(c)
    if output_path != c.output_path:
        print('Warning: Output has been moved from',c.output_path,'to',output_path)
        c.output_path = output_path
        
    #  loading surface variables
    try:
        df_out = xr.open_dataset(c.output_path+run_name+'/'+c.station+'_surface.nc').to_dataframe()
        df_in = df_in.loc[df_out.index[0]:df_out.index[-1],:]
    except Exception as e:
        print(e)
    # %% plotting surface variables
    try:
        lpl.plot_summary(df_out, c, 'SEB_output')
    except Exception as e:
        print(e)
        
    # %% plotting subsurface variables
    for var in ['compaction','T_ice','density_bulk','slwc','dgrain']:
        try:
            if len(df_in) <300: 
                ylim = [10]
            else:
                ylim = []
            lpl.plot_var(c.station, c.output_path, c.RunName, var, ylim=ylim, zero_surf=False)
        except Exception as e:
            print('lpl.plot_var(c.station, c.output_path, c.RunName, var, zero_surf=False)')
            print(var)
            print(e)
            pass

    if c.station in ['DY2']:
            lpl.plot_var(c.station, c.output_path, c.RunName, 'slwc', 
                         zero_surf=True, ylim=(8,0), year = (2012, 2024))

    # %% Start/end plots
    lpl.plot_var_start_end(c, 'T_ice')
    lpl.plot_var_start_end(c, 'density_bulk')

    # %% Mass balance components
    fig = plt.figure()
    ax=plt.gca()
    df_out.melt_mweq.cumsum().plot(ax=ax, label='melt')
    df_out.Snowfallmweq.cumsum().plot(ax=ax, label='Snowfall')
    df_out.zrogl.cumsum().plot(ax=ax, label='runoff')
    df_out.zrfrz_sum.cumsum().plot(ax=ax, label='refreezing')
    df_out.snowthick.plot(ax=ax, label='snowthickness')
    plt.legend()
    fig.savefig(c.output_path+c.RunName+'/'+c.station+'_SMB.png', dpi=120)
    
    
    # %% Surface height evaluation
    # extracting surface height
    if 'surface_height' not in df_out.columns:
        filename = c.output_path + run_name + "/" + c.station + "_T_ice.nc"
        ds = xr.open_dataset(filename).transpose()
        ds['surface_height'] = (ds.depth.isel(level=-1)
                -ds.depth.isel(level=-1).isel(time=0)
                -(ds.depth.isel(level=-1).diff(dim='time')
                  .where(ds.depth.isel(level=-1)
                         .diff(dim='time')>6,0).cumsum()))
        ds['surface_height'].values[0] = 0
        df_out['surface_height'] = ds['surface_height'].values
        del ds
    # try:
    #     path_aws_l3 = 'C:/Users/bav/GitHub/PROMICE data/aws-l3-dev/level_3/'
    #     df_obs = pd.read_csv(path_aws_l3+c.station+'/'+c.station+'_hour.csv')
    # except:
    #     print('cannot find L3 AWS file')
    #     df_obs = pd.DataFrame()
    #     pass
    
    path_aws_l4 = '../PROMICE/PROMICE-AWS-toolbox/out/L4/'
    if os.path.isfile(path_aws_l4+c.station+'_L4_ext.csv'):
        df_obs = pd.read_csv(path_aws_l4+c.station+'_L4_ext.csv')
    elif os.path.isfile(path_aws_l4+c.station+'_L4.csv'):
        df_obs = pd.read_csv(path_aws_l4+c.station+'_L4.csv')
    else:
        path_aws_l4 = '../PROMICE/GC-Net-Level-1-data-processing/L1/daily/'
        if os.path.isfile(path_aws_l4+c.station+'_daily.csv'):
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
        else:
            print('No observation for', c.station)
            return []
        # else:
        #     tmp = pd.DataFrame()
    # if len(df_obs)>0:
    #     df_obs.time= pd.to_datetime(df_obs.time).dt.tz_convert(None)
    #     if len(tmp)>0:
            
    df_obs.time= pd.to_datetime(df_obs.time).dt.tz_convert(None)
    df_obs = df_obs.set_index('time')
    # df_obs = df_obs.resample(pd.infer_freq(df_out.index)).mean()
    
    fig = plt.figure()
    tmp = (df_obs.z_surf_combined -df_out.surface_height).mean()
    plt.plot(df_obs.index, df_obs.z_surf_combined-tmp, 
              marker='.',ls='None', label='AWS')
    plt.plot(df_out.index, df_out.surface_height,color='tab:red',
             label='model')
    plt.legend()
    plt.ylabel('Surface height (m)')
    plt.title(c.station)
    fig.savefig(c.output_path+c.RunName+'/'+c.station+'surface_height.png', dpi=120)
    
    # %% calculating modelled t_i_10m
    # from scipy.interpolate import interp1d
    # from tqdm import tqdm  # Import tqdm for the progress bar  
    # def interpolate_temperature(dates, depth_cor, temp,  depth=10,
    #     min_diff_to_depth=2,  kind="linear" ):
    #     depth_cor = depth_cor.astype(float)
    #     df_interp = pd.DataFrame()
    #     df_interp["date"] = dates
    #     df_interp["temperatureObserved"] = np.nan
    
    #     # preprocessing temperatures for small gaps
    #     tmp = pd.DataFrame(temp)
    #     tmp["time"] = dates
    #     tmp = tmp.set_index("time")
    #     tmp = tmp.resample("H").mean()

    #     temp = tmp.loc[dates].values
    #     for i in tqdm(range(len(dates)), desc="Interpolating temperatures", unit="date"):
    #         x = depth_cor[i, :].astype(float)
    #         y = temp[i, :].astype(float)
    #         ind_no_nan = ~np.isnan(x + y)
    #         x = x[ind_no_nan]
    #         y = y[ind_no_nan]
    #         x, indices = np.unique(x, return_index=True)
    #         y = y[indices]
    #         if len(x) < 2 or np.min(np.abs(x - depth)) > min_diff_to_depth:
    #             continue
    #         f = interp1d(x, y, kind, fill_value="extrapolate")
    #         df_interp.iloc[i, 1] = np.min(f(depth), 0)
    
    #     if df_interp.iloc[:5, 1].std() > 0.1:
    #         df_interp.iloc[:5, 1] = np.nan
    #     return df_interp
    def interpolate_temperature_fast(dates, depth_matrix, temp_matrix,  depth=10,
        min_diff_to_depth=2,  kind="linear" ):
        # Choose the depth you want to interpolate to (e.g., 10 meters)
        target_depth = 10
        N = depth_matrix.shape[0]
        M = depth_matrix.shape[1]
        closest_depth_indices = np.abs(depth_matrix - target_depth).argmin(axis=1)
        closest_depths_idx_1 = np.maximum(0, closest_depth_indices - 1)
        closest_depths_idx_2 = np.minimum(M - 1, closest_depth_indices + 1)
        closest_depths = depth_matrix[np.arange(N), closest_depths_idx_1]
        next_closest_depths = depth_matrix[np.arange(N), closest_depths_idx_2]
        
        temp_at_closest_depths = temp_matrix[np.arange(N), closest_depths_idx_1]
        temp_at_next_closest_depths = temp_matrix[np.arange(N), closest_depths_idx_2]
        
        weights = (next_closest_depths - target_depth) / (next_closest_depths - closest_depths)
        temp_at_10m = temp_at_closest_depths + weights * (temp_at_next_closest_depths - temp_at_closest_depths)
        return temp_at_10m
    
    filename = c.output_path + run_name + "/" + c.station + "_T_ice.nc"
    df = (xr.open_dataset(filename).to_dataframe().unstack('level'))
    df.columns = df.columns.map('{0[0]}_{0[1]}'.format)
    # df_10m = interpolate_temperature(
    #     df.index, df[[v for v in df.columns if 'depth' in v]].values, 
    #     df[[v for v in df.columns if 'T_ice' in v]].values-273.15, 
    # )
    # df_out['t_i_10m'] = df_10m.temperatureObserved.values
    df_out['t_i_10m'] = interpolate_temperature_fast(
        df.index, df[[v for v in df.columns if 'depth' in v]].values, 
        df[[v for v in df.columns if 'T_ice' in v]].values-273.15, 
    )
    

    # %% plotting ['t_surf','LRout','LHF','SHF','t_i_10m']
    df_out['t_surf']  =  ((df_out.LRout_mdl - (1 -  c.em) * df_out.LRin) / c.em / 5.67e-8)**0.25 - 273.15
    df_out['LRout'] = df_out.LRout_mdl
    df_obs['LHF'] = df_obs.dlhf_u
    df_obs['SHF'] = df_obs.dshf_u
    if 'ulr' in df_obs.columns:
        df_obs['LRout'] = df_obs.ulr
    else:
        df_obs['LRout'] = np.nan
    var_list = ['t_surf','LRout','LHF','SHF','t_i_10m']
    
    from matplotlib import gridspec
    from scipy.stats import linregress
    fig = plt.figure(figsize=(12, 17))
    gs = gridspec.GridSpec(len(var_list), 2, width_ratios=[3, 1]) 
    
    df_obs = df_obs[~df_obs.index.duplicated(keep='first')]
    df_out = df_out[~df_out.index.duplicated(keep='first')]
    common_idx = df_obs.index.intersection(df_out.index)
    for i, var in enumerate(var_list): 
        if var not in df_obs.columns:
            df_obs[var] = np.nan
        ax1 = plt.subplot(gs[i, 0])
        ax2 = plt.subplot(gs[i, 1])
        # first plot
        ME = np.mean(df_out.loc[common_idx, var] - df_obs.loc[common_idx, var])
        RMSE = np.sqrt(np.mean((df_out.loc[common_idx, var] - df_obs.loc[common_idx, var])**2))                           
        df_obs[var].plot(ax=ax1, label='AWS',marker='.',markersize=2)
        df_out[var].plot(ax=ax1, alpha=0.7, label='SEB model')
        ax1.set_ylabel(var)
        ax1.set_xlim(common_idx.min(), common_idx.max())
        ax1.grid()
        if i == 0:  ax1.set_title(c.station+'\n\n')
        ax1.legend()
    
        # second plot
        ax2.plot(df_obs.loc[common_idx,var], df_out.loc[common_idx,var], 
                 color='k',alpha=0.1,marker='.',ls='None')
        ax2.set_xlabel('AWS')
        ax2.set_ylabel('SEB model')        
        common_idx = df_obs.loc[df_obs[var].notnull()].index.intersection(df_out.loc[df_out[var.replace('_uncor','')].notnull()].index)
    
        try:
            slope, intercept, r_value, p_value, std_err = linregress(
                df_obs.loc[common_idx, var], df_out.loc[common_idx, var])
            max_val = max(df_obs.loc[common_idx,var].max(), df_out.loc[common_idx,var].max())
            min_val = min(df_obs.loc[common_idx,var].min(), df_out.loc[common_idx,var].min())
            ax2.plot([min_val, max_val], [min_val, max_val], 'k-', label='1:1 Line')
            regression_line = slope * df_obs[var] + intercept
            ax2.plot(df_obs[var], regression_line, 'r-', label='Linear Regression')
        except:
            pass
        if i == 0: ax2.legend(loc='lower right')
        ax2.grid()
        
        # Annotate with RMSE and ME
        ax2.annotate(f'{var}\nRMSE: {RMSE:.2f}\nME: {ME:.2f}', 
                     xy=(1.05, 0.95), xycoords='axes fraction', 
                     horizontalalignment='left', verticalalignment='top',
                     fontsize=10, bbox=dict(boxstyle="round,pad=0.3",
                                            edgecolor='black', facecolor='white'))

    fig.savefig(c.output_path+c.RunName+'/SEB_evaluation_vs_AWS.png', dpi=120)
    
   

    
    # %% Movies
    # lpl.plot_movie(c.station, c.output_path, c.RunName, 'T_ice')
    # lpl.plot_movie(c.station, c.output_path, c.RunName, 'density_bulk')
    lpl.evaluate_temperature_sumup(df_out, c)
    lpl.evaluate_density_sumup(c)
    lpl.evaluate_accumulation_snowfox(df_in, c)
    try:
        lpl.evaluate_compaction(c)
    except Exception as e:
        print(e)
        pass
    try:
        lpl.find_summer_surface_depths(c)
    except Exception as e:
        print(e)
        pass

# %%
import os    
if __name__ == "__main__":
    # for run_name in os.listdir('C:/Users/bav/data_save/output firn model/'):
    main(output_path=output_path, run_name=run_name)
