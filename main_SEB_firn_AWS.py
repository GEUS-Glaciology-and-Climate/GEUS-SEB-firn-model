# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

Runs the GEUS SEB-firn model forced by gap-filled PROMICE level 3 AWS data,
either inverting the surface temperature from the surface energy balance
(Tsurf_mode = 'inverted', standard model) or prescribing it from the observed
upward longwave radiation (Tsurf_mode = 'observed').

Before running:
    python tools/prepare_aws_forcing.py         -> input/weather data/AWS_L3/
    python -m tools.prepare_aws_initial_state   -> input/initial state/AWS_L3/
"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import lib.io as io
from lib.initialization import ImportConst
from lib.seb_smb_model import GEUS_model, observed_surface_temperature
import os
import time
import traceback


def run_SEB_firn_AWS(
    station: str,
    Tsurf_mode: str = 'inverted',
    output_path: str = './output/AWS_L3/',
    forcing_folder: str = './input/weather data/AWS_L3/',
    initial_state_folder: str = './input/initial state/AWS_L3/',
    Tsurf_obs_melt_tol: float = 0.5,
    GF_formulation: str = 'original',
    output_var: str = "all",
):
    start_time = time.time()
    c = ImportConst()
    c.station = station
    c.spin_up = False
    c.verbose = 1
    c.Tsurf_mode = Tsurf_mode
    c.Tsurf_obs_melt_tol = Tsurf_obs_melt_tol
    c.GF_formulation = GF_formulation

    c.surface_input_path = forcing_folder + station + '_forcing.csv'
    c.surface_input_driver = 'AWS_L3'
    c.output_path = output_path

    c.use_spin_up_init = True
    c.spin_up_path = initial_state_folder
    c.initial_state_pkl = initial_state_folder + station + '_initial_state.pkl'
    ini = pd.read_csv(initial_state_folder + station + '_initial_scalars.csv', index_col=0).value
    c.snowthick_ini = ini.snowthick_ini
    c.Tdeep = ini.Tdeep

    c.num_lay = 100
    c.lim_new_lay = 0.05
    c.z_max = 70
    c.RunName = f"{station}_Tsurf_{Tsurf_mode}"
    if Tsurf_mode == 'observed' and Tsurf_obs_melt_tol != 0.5:
        c.RunName += f"_tol{Tsurf_obs_melt_tol}"
    if GF_formulation != 'original':
        c.RunName += f"_GF{GF_formulation}"
    os.makedirs(os.path.join(c.output_path, c.RunName), exist_ok=True)

    df_in, c = io.load_surface_input_data(c)
    c.delta_time = pd.to_timedelta('1h').total_seconds()
    c.zdtime = c.delta_time
    if c.altitude < 1500:
        c.new_bottom_lay = 30

    print(c.RunName, 'Tdeep', c.Tdeep, 'start/end', df_in.index[0], df_in.index[-1])
    df_out = pd.DataFrame(index=df_in.index)
    df_out.index.name = 'time'

    try:
        (
            df_out["L"], df_out["LHF"], df_out["SHF"], df_out["theta_2m"],
            df_out["q_2m"], df_out["ws_10m"], df_out["Re"], df_out["melt_mweq"],
            df_out["sublimation_mweq"], df_out["SRin"], df_out["SRout"],
            df_out["LRin"], df_out["LRout_mdl"],
            snowc, snic, slwc, T_ice,
            zrfrz, rhofirn, zsupimp, dgrain,
            df_out['zrogl'], df_out['Tsurf'],
            grndc, grndd, pgrndcapc, df_out['GF'],
            dH_comp, df_out["snowbkt"], compaction, df_out['snowthick']
        ) = GEUS_model(df_in, c)
    except Exception as e:
        print(c.station, e); traceback.print_exc()
        with open(f"{c.output_path}/{c.RunName}/error.txt", "a+") as text_file:
            text_file.write(f"{c.RunName}\n{e}\n{traceback.format_exc()}\n")
        return
    print(c.RunName, 'SEB firn model took %0.03f sec' % (time.time() - start_time))

    c.write(c.output_path + "/" + c.RunName + "/constants.csv")
    thickness_act = snowc * c.rho_water / rhofirn + snic * c.rho_water / c.rho_ice
    depth_act = np.cumsum(thickness_act, 0)
    density_bulk = (snowc + snic) / (snowc / rhofirn + snic / c.rho_ice)

    depth_lowest_layer = np.insert(np.diff(depth_act[-1, :]), 0, 0)
    depth_lowest_layer[np.abs(depth_lowest_layer) <= 3] = 0
    df_out['surface_height'] = (depth_act[-1, :] - depth_act[-1, 0]) - depth_lowest_layer.cumsum()
    df_out['snowfall_mweq'] = df_in.Snowfallmweq
    df_out['rainfall_mweq'] = df_in.Rainfallmweq
    df_out['refreezing_mweq'] = zrfrz.sum(axis=0)
    df_out = df_out.rename(columns={'zrogl': 'runoff_mweq'})
    df_out['smb_mweq'] = df_out.snowfall_mweq + df_out.rainfall_mweq \
        - df_out.runoff_mweq - df_out.sublimation_mweq
    df_out['Tsurf_obs'] = observed_surface_temperature(df_in, c)
    # meltflux in W m-2, back-calculated from melt
    df_out['meltflux'] = df_out.melt_mweq * c.L_fus * c.rho_water / c.zdtime

    io.write_1d_netcdf(df_out, c)
    vars_map = {
        "snic": snic, "snowc": snowc, "slwc": slwc, "rhofirn": rhofirn,
        "density_bulk": density_bulk, "T_ice": T_ice, "rfrz": zrfrz,
        "dgrain": dgrain, "compaction": compaction,
    }
    for name, data in vars_map.items():
        if output_var == "all" or name in output_var:
            io.write_2d_netcdf(data, name, depth_act, df_in.index, c)
    print(c.RunName, 'done in %0.01f sec' % (time.time() - start_time))


if __name__ == "__main__":
    for station in ['TAS_A', 'KAN_U', 'KAN_M', 'QAS_M', 'QAS_U']:
        for Tsurf_mode in ['inverted', 'observed']:
            run_SEB_firn_AWS(station, Tsurf_mode)
        # sensitivity to the margin below 0 degC at which an observed surface
        # is considered melting
        for tol in [0.0, 1.0]:
            run_SEB_firn_AWS(station, 'observed', Tsurf_obs_melt_tol=tol)
