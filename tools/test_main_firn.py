import __init__
import netCDF4 as nc
import numpy as np
import pandas as pd
from main_firn import run_main_firn_old, run_GEUS_model_opt, run_main_firn_parallel, run_GEUS_model_old
from joblib import Parallel, delayed

def compare_old_opt_firn():
    print("--- Test firn old running ---")
    site_list = ["KAN_U"]

    for site in site_list:
        print(site)
        filename = "Firn viscosity/Input files/" + site + ".csv"
        (snowc_old, snic_old, slwc_old, tsoil_old, zrfrz_old, 
        rhofirn_old, zsupimp_old, dgrain_old, grndc_old, 
        grndd_old, compaction_old, zrogl_old, Tsurf_out_old,
        pgrndcapc_old, pgrndhflx_old, dH_comp_old, snowbkt_old) = run_GEUS_model_old(site, filename)

    output_old = [snowc_old, snic_old, slwc_old, tsoil_old, zrfrz_old, 
        rhofirn_old, zsupimp_old, dgrain_old, grndc_old, 
        grndd_old, compaction_old, zrogl_old, Tsurf_out_old,
        pgrndcapc_old, pgrndhflx_old, dH_comp_old, snowbkt_old]

    # Run new main_firn, parallel
    res = Parallel(n_jobs=8, verbose=10)(delayed(run_GEUS_model_opt)(site) for site in site_list)
    output_opt = res[0]

    diff_exist = False
    # Assert that the values from original and optimized code are the same
    for i in range(len(output_opt)):
        if (output_old[i] != output_opt[i]).all():
            print("Difference in objects for i: ", i)
            diff_exist = True

    if diff_exist == False:
        print("Test done. No differences exist.")
   


def compare_ncfiles():
    print("--- Test firn running by comparing output nc files ---")
    site_list = ["KAN_U", "KAN_M"]
    variables = ['_density_bulk.nc','_compaction.nc','_rhofirn.nc', '_slwc.nc', '_snic.nc', '_snowc.nc', '_T_ice.nc']
    empty_list = [None] * len(site_list)
    (density_opt, compaction_opt, rhofirn_opt, slwc_opt, snic_opt, snowc_opt, T_ice_opt) = (empty_list, empty_list, empty_list, empty_list, empty_list, empty_list, empty_list)
    (density_old, compaction_old, rhofirn_old, slwc_old, snic_old, snowc_old, T_ice_old) = (empty_list, empty_list, empty_list, empty_list, empty_list, empty_list, empty_list)

    # Running in parallel, nc files are created
    runname_opt = Parallel(n_jobs=8, verbose=10)(delayed(run_GEUS_model_opt)(site) for site in site_list)
    
    # Reading the created nc files
    for i in range(0, len(site_list)):
        folder_path = 'C:/Users/brink/Documents/Exjobb/GEUS-SEB-firn-model/Output main_firn/' + runname_opt[i] + '/'
        for variable in variables:
            var_opt_path =  folder_path + site_list[i] + variable
            var_opt_ds = nc.Dataset(var_opt_path)
            if variable == '_density_bulk.nc':
                density_opt[i] = var_opt_ds['density_bulk'][:]
            if variable == '_compaction.nc':
                compaction_opt[i] = var_opt_ds['compaction'][:]
            if variable == '_rhofirn.nc':
                rhofirn_opt[i] = var_opt_ds['rhofirn'][:]
            if variable == '_slwc.nc':
                slwc_opt[i] = var_opt_ds['slwc'][:]
            if variable == '_snic.nc':
                snic_opt[i] = var_opt_ds['snic'][:]
            if variable == '_snowc.nc':
                snowc_opt[i] = var_opt_ds['snowc'][:]
            if variable == '_T_ice.nc':
                T_ice_opt[i] = var_opt_ds['T_ice'][:]

    # Running old script, in serie, nc files are created
    runname_old = [None] * len(site_list)
    for i in range(0, len(site_list)):
        filename = "Firn viscosity/Input files/" + site_list[i] + ".csv"
        print("Old code, for site:", site_list[i], ' and run name before assigned: ', runname_old[i])
        runname_old[i] = run_GEUS_model_old(site_list[i], filename) 
    
    # Reading the created nc files
    for i in range(0, len(site_list)):
        folder_path = 'C:/Users/brink/Documents/Exjobb/GEUS-SEB-firn-model/Output main_firn/old_script/' + runname_old[i] + '/'
        for variable in variables:
            var_old_path =  folder_path + site_list[i] + variable
            var_old_ds = nc.Dataset(var_old_path)
            if variable == '_density_bulk.nc':
                density_old[i] = var_old_ds['density_bulk'][:]
            if variable == '_compaction.nc':
                compaction_old[i] = var_old_ds['compaction'][:]
            if variable == '_rhofirn.nc':
                rhofirn_old[i] = var_old_ds['rhofirn'][:]
            if variable == '_slwc.nc':
                slwc_old[i] = var_old_ds['slwc'][:]
            if variable == '_snic.nc':
                snic_old[i] = var_old_ds['snic'][:]
            if variable == '_snowc.nc':
                snowc_old[i] = var_old_ds['snowc'][:]
            if variable == '_T_ice.nc':
                T_ice_old[i] = var_old_ds['T_ice'][:]

    diff = False
    if density_opt != density_old:
        diff = True
    if compaction_opt != compaction_old:
        diff = True
    if rhofirn_opt != rhofirn_old:
        diff = True
    if slwc_opt != slwc_old:
        diff = True
    if snic_opt != snic_old:
        diff = True
    if snowc_opt != snowc_old:
        diff = True
    if T_ice_opt != T_ice_old:
        diff = True
    if diff == True:
        print("Differences exist.")
    else:
        print("No differences. Test OK.")

if __name__ == '__main__':
    compare_ncfiles()
    #compare_old_opt_firn()
