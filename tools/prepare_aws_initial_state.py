# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

Builds the initial subsurface state for the AWS-driven runs at TAS_A and KAN_U.

The state is saved in the same pickle format as the spin-up final state
([snowc, snic, slwc, T_ice, rhofirn, dgrain, Tsurf, grndc, grndd, snowbkt])
and is read by lib.initialization.IniVar through c.initial_state_pkl.

KAN_U (percolation zone, ice slabs):
    layer structure (snow/ice content, density, water content, grain size)
    taken from the last time step of the CARRA-driven run (2025-08-31), which
    carries the ice slab built up by 1990-2025 melt. The temperature is then
    replaced by the AWS thermistor string (daily mean of the first day), and
    below the deepest thermistor by the observed 10 m temperature.

TAS_A, KAN_M, QAS_M, QAS_U (ablation zone):
    bare ice column (the stations recorded < 3 cm of snow at the start), on the
    default model grid, with temperatures from the thermistor string as above
    (first day with working thermistors). At QAS_M, the upper thermistors
    had melted out in 2025-2026 (no depth) and the lower ones were not
    recording, so the first usable profile is the one of 2026-09-01, used as a
    proxy for the end-of-summer profile of 2025.
"""
import os
import pickle
import numpy as np
import pandas as pd
import xarray as xr
from lib.initialization import ImportConst

forcing_folder = './input/weather data/AWS_L3/'
output_folder = './input/initial state/AWS_L3/'


def observed_temperature_profile(df_forcing, depth_m, t_surface=None):
    '''Interpolates the thermistor temperatures of the first day of the forcing
    with at least 3 working thermistors onto depth_m (layer mid-depths, m).
    Below the deepest working thermistor, the 10 m temperature is used.
    Returns temperature in degC.'''
    t_cols = [c for c in df_forcing.columns if c.startswith('t_i_') and c[4:].isdigit()]
    n_valid = (df_forcing[t_cols].notnull().values
               & df_forcing[['d_' + c for c in t_cols]].notnull().values).sum(axis=1)
    t_start = df_forcing.index[np.argmax(n_valid >= 3)]
    if t_start != df_forcing.index[0]:
        print('Thermistor profile taken from', t_start)
    first_day = df_forcing.loc[t_start:t_start + pd.Timedelta('1D')]
    temp = first_day[t_cols].mean()
    depth = first_day[['d_' + c for c in t_cols]].mean().values
    temp = temp.values
    ok = ~np.isnan(temp) & ~np.isnan(depth) & (depth > 0)
    temp, depth = temp[ok], depth[ok]
    order = np.argsort(depth)
    temp, depth = temp[order], depth[order]

    t10 = first_day.t_i_10m.mean()
    if np.isnan(t10):
        t10 = temp[-1]
    # constant (observed) value below the deepest thermistor
    depth = np.append(depth, [max(depth[-1] + 1, 10), 1000])
    temp = np.append(temp, [t10, t10])
    # surface: extrapolate the top thermistor, or interpolate linearly
    # towards t_surface (degC) when given
    if t_surface is not None:
        depth = np.insert(depth, 0, 0.)
        temp = np.insert(temp, 0, t_surface)
    T = np.interp(depth_m, depth, temp)
    return np.minimum(T, 0), t10


def default_grid_weq(c):
    # same layering as lib.initialization.InitializationSubsurface
    depth_mod_weq = np.insert(
        np.arange(1, c.num_lay + 1) ** 4 / c.num_lay ** 4 * c.z_max, 0, 0)
    thickness = np.diff(depth_mod_weq)
    tmp = np.maximum(0, c.lim_new_lay - thickness)
    return thickness + tmp - np.flip(tmp)


def init_ablation_site(station, c):
    df = pd.read_csv(forcing_folder + station + '_forcing.csv', index_col=0, parse_dates=True)
    thickness_weq = default_grid_weq(c)
    snic = thickness_weq.copy()
    snowc = np.zeros_like(snic)
    slwc = np.zeros_like(snic)
    rhofirn = np.full_like(snic, 900.)
    dgrain = np.full_like(snic, 2.)
    depth_mid = np.cumsum(thickness_weq * c.rho_water / c.rho_ice) \
        - thickness_weq * c.rho_water / c.rho_ice / 2
    T_degC, t10 = observed_temperature_profile(df, depth_mid)
    return dict(snowc=snowc, snic=snic, slwc=slwc, T_ice=T_degC + c.T_0,
                rhofirn=rhofirn, dgrain=dgrain, depth_mid=depth_mid,
                snowbkt=0., snowthick=max(0., df.snow_height.iloc[:24].mean()),
                Tdeep=t10 + c.T_0)


def init_from_model_run(station, c, run_folder, t_surface=None):
    df = pd.read_csv(forcing_folder + station + '_forcing.csv', index_col=0, parse_dates=True)
    start = df.index[0]
    state = {}
    for var in ['snowc', 'snic', 'slwc', 'rhofirn', 'dgrain', 'T_ice']:
        with xr.open_dataset(f'{run_folder}/{station}_{var}.nc') as ds:
            ds_t = ds.sel(time=start, method='nearest')
            print(var, 'taken from', ds_t.time.values)
            state[var] = ds_t[var].values.astype('float64')
    with xr.open_dataset(f'{run_folder}/{station}_surface.nc') as ds:
        state['snowbkt'] = float(ds.snowbkt.sel(time=start, method='nearest'))

    thickness_m = state['snowc'] * c.rho_water / state['rhofirn'] \
        + state['snic'] * c.rho_water / c.rho_ice
    depth_mid = np.cumsum(thickness_m) - thickness_m / 2
    T_degC, t10 = observed_temperature_profile(df, depth_mid, t_surface)
    state['T_ice_model'] = state['T_ice']
    state['T_ice'] = T_degC + c.T_0
    state['depth_mid'] = depth_mid
    # seasonal snow on top of the previous summer surface is negligible end of August
    state['snowthick'] = 0.
    state["Tdeep"] = t10 + c.T_0
    return state


def save_state(station, state, c):
    os.makedirs(output_folder, exist_ok=True)
    T_ice = state['T_ice']
    with open(output_folder + station + '_initial_state.pkl', 'wb') as f:
        pickle.dump([state['snowc'], state['snic'], state['slwc'], T_ice,
                     state['rhofirn'], state['dgrain'], T_ice[0],
                     T_ice.copy(),  # grndc, as in IniVar
                     np.zeros_like(T_ice),  # grndd
                     state['snowbkt']], f)
    # scalars needed by the run
    pd.Series({'snowthick_ini': state['snowthick'], 'Tdeep': state['Tdeep']}).to_csv(
        output_folder + station + '_initial_scalars.csv', header=['value'])

    density_bulk = (state['snowc'] + state['snic']) / (
        state['snowc'] / state['rhofirn'] + state['snic'] / c.rho_ice)
    pd.DataFrame({'depth_mid_m': state['depth_mid'],
                  'density_bulk': density_bulk,
                  'T_ice_degC': T_ice - c.T_0,
                  'snowc': state['snowc'], 'snic': state['snic'],
                  'slwc': state['slwc']}).to_csv(
        output_folder + station + '_initial_profile.csv', index=False, float_format='%.5g')


def plot_state(station, state, c):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    density_bulk = (state['snowc'] + state['snic']) / (
        state['snowc'] / state['rhofirn'] + state['snic'] / c.rho_ice)
    fig, ax = plt.subplots(1, 2, sharey=True, figsize=(8, 6))
    ax[0].step(density_bulk, state['depth_mid'], where='mid')
    ax[0].set_xlabel('bulk density (kg m-3)')
    ax[1].plot(state['T_ice'] - c.T_0, state['depth_mid'], marker='.', label='initial state')
    if 'T_ice_model' in state:
        ax[1].plot(state['T_ice_model'] - c.T_0, state['depth_mid'], ls='--',
                   label='CARRA-driven run (replaced)')
        ax[1].legend()
    ax[1].set_xlabel('temperature (degC)')
    ax[0].set_ylim(20, 0)
    ax[0].set_ylabel('depth (m)')
    fig.suptitle(station + ' initial state')
    fig.savefig(output_folder + station + '_initial_state.png', dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    c = ImportConst()
    c.num_lay = 100
    c.z_max = 70
    c.lim_new_lay = 0.05

    for station in ['TAS_A', 'KAN_M', 'QAS_M', 'QAS_U']:
        state = init_ablation_site(station, c)
        save_state(station, state, c)
        plot_state(station, state, c)

    state = init_from_model_run('KAN_U', c, './output/2026-03-16/KAN_U_100_layers_3h')
    save_state('KAN_U', state, c)
    plot_state('KAN_U', state, c)
