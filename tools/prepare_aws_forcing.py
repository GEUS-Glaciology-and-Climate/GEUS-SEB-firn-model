# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

Prepares gap-free hourly model forcing from PROMICE/GC-Net level 3 AWS data.

Gap-filling, applied variable by variable, in this order:
    1. gaps of up to `max_interp_gap` hours are linearly interpolated
    2. tilt-corrected shortwave radiation is filled from the station's own
       uncorrected radiometer through a zero-intercept regression
    3. longer gaps are filled from neighbour AWSs through linear regressions
       fitted on the hours where both stations have data (zero intercept for
       shortwave radiation); when several neighbours are available, their
       predictions are averaged
    4. what remains is filled with the mean diurnal cycle of the station
       itself, computed over a centred window of +/- `clim_window` days

Upward shortwave radiation is not gap-filled directly: daily albedo is
gap-filled (neighbour albedo + mean offset, then interpolation) and the
missing usr is recomputed as albedo x gap-filled dsr.

Snowfall is derived from the positive daily increments of the median-smoothed
surface height (z_surf_combined) converted to water equivalent with the fresh snow density.
Rainfall is the corrected rain gauge record, set to 0 where missing.

Observed upward longwave radiation is kept as is (not gap-filled) so that the
model can use the observed surface temperature when it is available.

Each output variable has an `<var>_origin` column:
    0 = observed, 1 = interpolated, 2 = neighbour station, 3 = diurnal climatology,
    4 = station's own uncorrected radiometer
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

path_l3 = "C:/Users/bav/GitHub/PROMICE data/thredds-data/level_3_sites/csv/hour/"
output_folder = "./input/weather data/AWS_L3/"

neighbours = {
    'TAS_A': ['TAS_L'],
    'KAN_U': ['KAN_M'],
    'KAN_M': ['KAN_U', 'KAN_L'],
    'QAS_M': ['QAS_U', 'QAS_L'],
    'QAS_U': ['QAS_M', 'QAS_L'],
}

# smoothing window applied to surface height before deriving snowfall.
# At the ablation sites, 1-day increments give 2-3 times more snowfall than
# the maximum observed snow depth can hold (e.g. 2.7 m w.e. for a 1.6 m deep
# snowpack at TAS_A), the SR50 noise and wind drift being counted as snowfall.
# 5-day smoothing gives values consistent with the observed snow depth.
snowfall_smoothing_default = '5D'
snowfall_smoothing = {'KAN_U': '1D'}

# model input name: L3 variable name
var_map = {
    'AirTemperatureC': 't_u',
    'SpecificHumiditykgkg': 'qh_u',
    'WindSpeedms': 'wspd_u',
    'AirPressurehPa': 'p_u',
    'ShortwaveRadiationDownWm2': 'dsr_cor',
    'LongwaveRadiationDownWm2': 'dlr',
}

# neighbour variable used as predictor when different from the station variable
# (the uncorrected radiation has better coverage than the tilt-corrected one,
# and the regression absorbs the tilt bias)
neighbour_var_map = {'dsr_cor': 'dsr', 'usr_cor': 'usr'}

# station's own uncorrected radiometer used to fill the tilt-corrected values
self_proxy_var_map = {'dsr_cor': 'dsr'}

ORIGIN_OBS, ORIGIN_INTERP, ORIGIN_NEIGHBOUR, ORIGIN_CLIM, ORIGIN_RAW = 0, 1, 2, 3, 4


def read_l3(station, start, end):
    df = pd.read_csv(path_l3 + station + '_hour.csv', index_col=0, parse_dates=True)
    df = df.loc[~df.index.duplicated()]
    return df.reindex(pd.date_range(start, end, freq='h'))


def interpolate_short_gaps(s, max_gap):
    # only fills NaN runs that are at most max_gap long and bounded by data
    isnull = s.isnull()
    run_id = (isnull != isnull.shift()).cumsum()
    run_len = isnull.groupby(run_id).transform('sum')
    filled = s.interpolate(method='time', limit_area='inside')
    return s.where(~(isnull & (run_len <= max_gap)), filled)


def regression_prediction(s, s_nb, zero_intercept=False, min_overlap=500):
    # prediction of s from s_nb, NaN where s_nb is missing or the overlap is too short
    both = s.notnull() & s_nb.notnull()
    if both.sum() < min_overlap:
        return pd.Series(np.nan, index=s.index)
    x, y = s_nb[both].values, s[both].values
    if zero_intercept:
        a, b = np.sum(x * y) / np.sum(x * x), 0.
    else:
        a, b = np.polyfit(x, y, 1)
    return a * s_nb + b


def fill_from_neighbour(s, s_nb_list, zero_intercept=False):
    # average of the predictions from all neighbours
    pred = pd.concat([regression_prediction(s, s_nb, zero_intercept)
                      for s_nb in s_nb_list], axis=1).mean(axis=1)
    return s.fillna(pred)


def fill_diurnal_climatology(s, clim_window):
    # mean value at each hour of the day over a centred window of +/- clim_window days
    tmp = s.to_frame('v')
    tmp['date'] = tmp.index.floor('D')
    tmp['hour'] = tmp.index.hour
    piv = tmp.pivot_table(index='date', columns='hour', values='v', dropna=False)
    piv = piv.reindex(pd.date_range(tmp.date.min(), tmp.date.max(), freq='D'))
    out = s.copy()
    for w in [clim_window, 3 * clim_window, 365]:
        clim = piv.rolling(2 * w + 1, center=True, min_periods=1).mean()
        clim_h = clim.stack(future_stack=True).reindex(
            list(zip(tmp.date, tmp.hour))).values
        out = out.where(out.notnull(), pd.Series(clim_h, index=s.index))
        if out.notnull().all():
            break
    return out


def gap_fill(s, s_nb_list=(), s_self_proxy=None, max_interp_gap=6, clim_window=7,
             zero_intercept=False):
    origin = pd.Series(np.where(s.notnull(), ORIGIN_OBS, -1), index=s.index)

    s = interpolate_short_gaps(s, max_interp_gap)
    origin[(origin == -1) & s.notnull()] = ORIGIN_INTERP

    if s_self_proxy is not None:
        s = s.fillna(regression_prediction(s, s_self_proxy, zero_intercept))
        origin[(origin == -1) & s.notnull()] = ORIGIN_RAW

    if len(s_nb_list) > 0:
        s = fill_from_neighbour(s, s_nb_list, zero_intercept=zero_intercept)
        origin[(origin == -1) & s.notnull()] = ORIGIN_NEIGHBOUR

    s = fill_diurnal_climatology(s, clim_window)
    origin[(origin == -1) & s.notnull()] = ORIGIN_CLIM
    return s, origin


def daily_albedo(df, dsr='dsr_cor', usr='usr_cor'):
    ok = df[dsr].notnull() & df[usr].notnull() & (df[dsr] > 50)
    day = df.loc[ok, [dsr, usr]].resample('D').sum()
    n_hours = ok.resample('D').sum()
    alb = (day[usr] / day[dsr]).where(n_hours >= 4)
    return alb.where((alb > 0.1) & (alb < 0.95))


def gap_fill_albedo(df, df_nb_list):
    alb = daily_albedo(df)
    alb = alb.reindex(pd.date_range(df.index[0].floor('D'), df.index[-1].floor('D'), freq='D'))
    origin = pd.Series(np.where(alb.notnull(), ORIGIN_OBS, -1), index=alb.index)

    # station's own uncorrected radiometers
    alb = alb.fillna(daily_albedo(df, 'dsr', 'usr').reindex(alb.index))
    origin[(origin == -1) & alb.notnull()] = ORIGIN_RAW

    for df_nb in df_nb_list:
        alb_nb = daily_albedo(df_nb, 'dsr', 'usr').reindex(alb.index)
        # for each gap, the neighbour albedo is shifted by the mean offset
        # between the two stations over the 30 days preceding the gap
        isnull = alb.isnull()
        run_id = (isnull != isnull.shift()).cumsum()
        for _, gap in alb[isnull].groupby(run_id[isnull]):
            if len(gap) <= 5:
                continue
            before = alb.loc[gap.index[0] - pd.Timedelta('30D'):gap.index[0]]
            offset = (before - alb_nb.reindex(before.index)).mean()
            if np.isnan(offset):
                after = alb.loc[gap.index[-1]:gap.index[-1] + pd.Timedelta('30D')]
                offset = (after - alb_nb.reindex(after.index)).mean()
            if not np.isnan(offset):
                alb.loc[gap.index] = (alb_nb.loc[gap.index] + offset).clip(0.1, 0.95)
        origin[(origin == -1) & alb.notnull()] = ORIGIN_NEIGHBOUR

    alb = alb.interpolate(method='time').ffill().bfill()
    origin[origin == -1] = ORIGIN_INTERP
    return alb.reindex(df.index, method='ffill'), origin.reindex(df.index, method='ffill')


def snowfall_from_surface_height(z_surf, t_air, rho_snow=315, min_increment=0.01,
                                 max_increment=0.4, max_t_air=2, window='1D'):
    # positive daily increments of the (median-smoothed) surface height = snowfall
    # days without surface height observations get no snowfall.
    # Rises larger than max_increment (m/day) or on days warmer than max_t_air
    # (degC, daily mean) are considered artefacts (station maintenance, sensor
    # resets) and discarded.
    z = z_surf.rolling(window, center=True, min_periods=6).median()
    dz = z.resample('D').mean().diff().fillna(0)
    t_daily = t_air.resample('D').mean().reindex(dz.index)
    dz = dz.where((dz > min_increment) & (dz < max_increment) & (t_daily < max_t_air), 0)
    sf_daily = dz * rho_snow / 1000  # m w.e.
    return sf_daily.reindex(z_surf.index, method='ffill') / 24


def prepare_forcing(station, start, end, plot=True):
    df = read_l3(station, start, end)
    df_nb_list = [read_l3(nb, start, end) for nb in neighbours.get(station, [])]

    out = pd.DataFrame(index=df.index)
    for var, l3_var in var_map.items():
        nb_var = neighbour_var_map.get(l3_var, l3_var)
        s_nb_list = [df_nb[nb_var] for df_nb in df_nb_list if nb_var in df_nb]
        proxy = self_proxy_var_map.get(l3_var, None)
        out[var], out[var + '_origin'] = gap_fill(
            df[l3_var], s_nb_list,
            s_self_proxy=df[proxy] if proxy in df else None,
            zero_intercept=(l3_var == 'dsr_cor'))
        print(f'{station} {var}: ' + ', '.join(
            f'{k}={(out[var + "_origin"] == v).mean() * 100:.1f}%' for k, v in
            zip(['obs', 'interp', 'neighbour', 'clim', 'own raw'], range(5))))

    out['SpecificHumiditykgkg'] = out.SpecificHumiditykgkg / 1000  # g/kg -> kg/kg
    out['WindSpeedms'] = out.WindSpeedms.clip(lower=0.1)
    out['ShortwaveRadiationDownWm2'] = out.ShortwaveRadiationDownWm2.clip(lower=0)

    # upward shortwave
    out['Albedo'], out['Albedo_origin'] = gap_fill_albedo(df, df_nb_list)
    usr = interpolate_short_gaps(df.usr_cor, 6)
    out['ShortwaveRadiationUpWm2_origin'] = np.where(
        df.usr_cor.notnull(), ORIGIN_OBS, np.where(usr.notnull(), ORIGIN_INTERP, out.Albedo_origin))
    out['ShortwaveRadiationUpWm2'] = usr.fillna(out.Albedo * out.ShortwaveRadiationDownWm2)
    out['ShortwaveRadiationUpWm2'] = np.minimum(
        out.ShortwaveRadiationUpWm2.clip(lower=0), 0.95 * out.ShortwaveRadiationDownWm2)

    # observed upward longwave, kept with its gaps for the observed Tsurf option
    out['LongwaveRadiationUpWm2'] = df.ulr

    # instrument heights (from boom height, PROMICE geometry)
    z_boom = df.z_boom_cor_u.interpolate(method='time').ffill().bfill()
    z_boom = z_boom.rolling(24, center=True, min_periods=1).median().clip(lower=0.5)
    out['HeightTemperaturem'] = z_boom - 0.1
    out['HeightHumiditym'] = z_boom - 0.1
    out['HeightWindSpeedm'] = z_boom + 0.4

    # precipitation
    out['Snowfallmweq'] = snowfall_from_surface_height(
        df.z_surf_combined, out.AirTemperatureC, window=snowfall_smoothing.get(station, snowfall_smoothing_default))
    rain = df.rainfall_cor_u if 'rainfall_cor_u' in df else pd.Series(np.nan, index=df.index)
    out['Rainfallmweq'] = (rain.fillna(0) / 1000).clip(lower=0)
    out['Rainfallmweq_origin'] = np.where(rain.notnull(), ORIGIN_OBS, ORIGIN_CLIM)

    # observations used for initialisation and evaluation
    out['z_surf_combined'] = df.z_surf_combined
    if 'z_ice_surf' in df:
        out['z_ice_surf'] = df.z_ice_surf
    out['snow_height'] = df.snow_height
    out['t_surf_l3'] = df.t_surf
    out['t_i_10m'] = df.t_i_10m
    for col in [c for c in df.columns if c.startswith('t_i_') or c.startswith('d_t_i_')]:
        if col != 't_i_10m':
            out[col] = df[col]

    out['latitude'] = df.lat.mean()
    out['longitude'] = df.lon.mean()
    out['altitude'] = df.alt.mean()
    out.index.name = 'time'

    os.makedirs(output_folder, exist_ok=True)
    out.to_csv(output_folder + station + '_forcing.csv', float_format='%.6g')

    if plot:
        plot_forcing(out, station)
    return out


def plot_forcing(out, station):
    var_list = list(var_map.keys()) + ['ShortwaveRadiationUpWm2', 'Albedo',
                                       'Snowfallmweq', 'Rainfallmweq']
    colors = {ORIGIN_OBS: 'k', ORIGIN_INTERP: 'tab:blue',
              ORIGIN_NEIGHBOUR: 'tab:orange', ORIGIN_CLIM: 'tab:red', ORIGIN_RAW: 'tab:green'}
    labels = {ORIGIN_OBS: 'observed', ORIGIN_INTERP: 'interpolated',
              ORIGIN_NEIGHBOUR: 'neighbour AWS', ORIGIN_CLIM: 'diurnal climatology / zero',
              ORIGIN_RAW: 'own uncorrected radiometer'}
    fig, axs = plt.subplots(len(var_list), 1, sharex=True, figsize=(12, 2 * len(var_list)))
    for ax, var in zip(axs, var_list):
        s = out[var].resample('D').mean() if var not in ['Snowfallmweq', 'Rainfallmweq'] \
            else out[var].resample('D').sum()
        if var + '_origin' in out:
            org = out[var + '_origin'].resample('D').agg(lambda x: x.mode().iloc[0] if len(x) else np.nan)
            for o, col in colors.items():
                ax.plot(s.where(org == o), marker='.', ls='-', color=col, label=labels[o])
        else:
            ax.plot(s, color='k', label='observed')
        ax.set_ylabel(var.replace('Radiation', '').replace('Wm2', ''), fontsize=8)
        ax.grid()
    axs[0].legend(ncol=5, loc='lower left', bbox_to_anchor=(0, 1.05))
    axs[0].set_title(station + ' gap-filled forcing (daily values)', pad=25)
    fig.tight_layout()
    fig.savefig(output_folder + station + '_forcing.png', dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    # start of the TAS_A record in the 2025-2026 balance year
    prepare_forcing('TAS_A', '2025-09-09', '2026-10-01 00:00')
    # KAN_U starts when the CARRA-driven run used for initialisation ends
    prepare_forcing('KAN_U', '2025-09-01', '2026-10-01 00:00')
    prepare_forcing('KAN_M', '2025-09-01', '2026-10-01 00:00')
    # before 2025-09-26, the QAS_M SR50 surface is offset by 2.3 m from the
    # (bare) ice surface given by the pressure transducer
    prepare_forcing('QAS_M', '2025-09-26', '2026-10-01 00:00')
    prepare_forcing('QAS_U', '2025-09-01', '2026-10-01 00:00')
