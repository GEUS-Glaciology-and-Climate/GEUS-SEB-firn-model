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
import lib.plot as lpl

# output_path= 'C:/Users/bav/data_save/output firn model/spin up 3H/'
output_path = './output/new/'
run_name = 'CP1_100_layers_3h'
#%%
def main(output_path, run_name):
    # %% Loading data
    print(run_name)
    df_out, df_in, c = lpl.load_model_input_output(output_path, run_name)

    # plotting surface variables
    lpl.plot_summary(df_out, c, 'SEB_output')

    # %% plotting subsurface variables
    for var in ['density_bulk','slwc','rfrz']:
        lpl.plot_var(c.station, c.output_path, c.RunName, var,
                     ylim=[10] if len(df_out) <300 else [],
                     zero_surf=False)

    # if c.station in ['DY2', 'KAN_U','CP1']:
        # lpl.plot_var(c.station, c.output_path, c.RunName, 'slwc',
                     # zero_surf=True, ylim=(8,0), year = (2012, 2024))

    # %%
    df_obs, obs_avail = lpl.plot_surface_height_evaluation(df_out, c)

    if obs_avail:
        lpl.plot_observed_vars(df_obs, df_out, c, var_list = ['t_surf','LRout','LHF','SHF','t_i_10m'])

    lpl.plot_smb_components(df_out, c)
    lpl.evaluate_temperature_sumup(df_out, c)
    # lpl.evaluate_temperature_scatter(df_out, c, year = None)
    lpl.evaluate_density_sumup(c,path_to_SUMup='../../Data/SUMup/2025')
    lpl.evaluate_smb_sumup(df_out, c, path_to_SUMup='../../Data/SUMup/2025')
    lpl.evaluate_accumulation_snowfox(df_in, c)
    lpl.plot_var_start_end(c, 'T_ice')
    lpl.plot_var_start_end(c, 'density_bulk')
    # lpl.plot_movie(c.station, c.output_path, c.RunName, 'T_ice')
    # lpl.plot_movie(c.station, c.output_path, c.RunName, 'density_bulk')
    lpl.evaluate_compaction(c)

    plt.close('all')
    # try:
        # lpl.find_summer_surface_depths(c)
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
