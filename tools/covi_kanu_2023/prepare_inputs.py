# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

Prepares common inputs for running both the GEUS SEB-firn model and the
model of Covi et al. (2022, https://doi.org/10.1017/jog.2022.54,
https://github.com/fcovi/1D_SEB_model) at KAN_U for the 2023 melt season.

Both models get:
    - the same gap-filled AWS forcing (tools/prepare_aws_forcing.py)
    - the same initial firn state: density from the CARRA-driven GEUS run on
      the start date, temperature from the AWS thermistor string

Run from the GEUS-SEB-firn-model folder:
    python -m tools.covi_kanu_2023.prepare_inputs
"""
import os
import numpy as np
import pandas as pd
import xarray as xr
import tools.prepare_aws_forcing as paf
import tools.prepare_aws_initial_state as pis
from lib.initialization import ImportConst

station = 'KAN_U'
start, end = '2023-05-01', '2023-09-30 23:00'
carra_run = './output/2026-03-16/KAN_U_100_layers_3h'
covi_path = '../Covi_1D_SEB_model'
covi_firn_template = covi_path + '/input/firn/SiteJ2019_firninit_5cm35m.txt'

forcing_folder = './input/weather data/AWS_L3_KAN_U_2023/'
initial_state_folder = './input/initial state/AWS_L3_KAN_U_2023/'
covi_input_folder = covi_path + '/input/KAN_U_2023/'


def prepare_forcing():
    paf.output_folder = forcing_folder
    out = paf.prepare_forcing(station, start, end)

    # relative humidity (with regard to water, as used by Covi et al.)
    df = paf.read_l3(station, start, end)
    df_nb = paf.read_l3('KAN_M', start, end)
    rh, rh_origin = paf.gap_fill(df.rh_u, [df_nb.rh_u])
    out['RelativeHumidity'] = rh.clip(0, 100)
    out['RelativeHumidity_origin'] = rh_origin
    out.to_csv(forcing_folder + station + '_forcing.csv', float_format='%.6g')
    return out


def write_covi_met(out):
    met = pd.DataFrame(index=out.index)
    met.index.name = 'TIMESTAMP'
    met['Year'] = met.index.year
    met['JD'] = met.index.dayofyear
    met['Time'] = met.index.hour
    met['Tair'] = out.AirTemperatureC
    met['RH'] = out.RelativeHumidity
    met['Wind'] = out.WindSpeedms
    met['SWin'] = out.ShortwaveRadiationDownWm2
    met['SWout'] = out.ShortwaveRadiationUpWm2
    met['LWin'] = out.LongwaveRadiationDownWm2
    # the Covi model has no fallback for missing LWout: short gaps interpolated
    met['LWout'] = out.LongwaveRadiationUpWm2.interpolate(method='time').ffill().bfill()
    met['Net'] = met.SWin - met.SWout + met.LWin - met.LWout
    met['Prec'] = (out.Snowfallmweq + out.Rainfallmweq) * 1000  # mm w.e. per hour
    os.makedirs(covi_input_folder, exist_ok=True)
    with open(covi_input_folder + 'KAN_U2023_metinit.txt', 'w') as f:
        f.write('"KAN_U AWS (PROMICE L3, gap-filled) climate data to be used as model input"\n')
        met.to_csv(f, float_format='%.6f')
    print('LWout gaps interpolated:', int(out.LongwaveRadiationUpWm2.isnull().sum()), 'hours')


def surface_temperature_first_day(df_forcing):
    # mean surface temperature of the first day from observed LWout (emissivity 1).
    # The top thermistor is ~3 m deep at the start: the temperature between the
    # surface and that thermistor is interpolated linearly
    lw = df_forcing.LongwaveRadiationUpWm2.iloc[:24]
    return float(np.minimum((lw / 5.67e-8) ** 0.25 - 273.15, 0).mean())


def initial_state(c, t_surface):
    # GEUS: layers of the CARRA-driven run + observed temperature
    pis.forcing_folder = forcing_folder
    pis.output_folder = initial_state_folder
    state = pis.init_from_model_run(station, c, carra_run, t_surface)
    pis.save_state(station, state, c)
    pis.plot_state(station, state, c)
    return state


def write_covi_firn(state, df_forcing, c, t_surface):
    # Covi grid (5 cm to 1 m, 10 cm to 10 m, 50 cm to 35 m)
    grid = pd.read_csv(covi_firn_template, sep=' ', header=1)
    density_bulk = (state['snowc'] + state['snic']) / (
        state['snowc'] / state['rhofirn'] + state['snic'] / c.rho_ice)
    thickness_m = state['snowc'] * c.rho_water / state['rhofirn'] \
        + state['snic'] * c.rho_water / c.rho_ice
    depth_top = np.cumsum(thickness_m) - thickness_m
    # density: GEUS layer containing the middle of each Covi layer
    ind = np.searchsorted(depth_top, grid.layerdepth.values, side='right') - 1
    grid['layerdensity'] = np.minimum(density_bulk[ind], 900)
    T, _ = pis.observed_temperature_profile(df_forcing, grid.layerdepth.values, t_surface)
    grid['layertemperature'] = T
    with open(covi_input_folder + 'KAN_U2023_firninit_5cm35m.txt', 'w') as f:
        f.write('"KAN_U snow model initialization file: density from CARRA-driven '
                'GEUS run, temperature from AWS thermistors"\n')
        grid.to_csv(f, sep=' ', index=False, float_format='%.4f')
    return grid


if __name__ == "__main__":
    c = ImportConst()
    c.num_lay = 100
    c.z_max = 70
    c.lim_new_lay = 0.05

    out = prepare_forcing()
    write_covi_met(out)
    t_surface = surface_temperature_first_day(out)
    print('surface temperature for initialisation: %.2f degC' % t_surface)
    state = initial_state(c, t_surface)
    grid = write_covi_firn(state, out, c, t_surface)
    print(grid.iloc[::12].to_string())
