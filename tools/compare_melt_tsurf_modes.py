# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

Compares the melt calculated by the GEUS SEB-firn model when the surface
temperature is
    - inverted from the surface energy balance (standard model), or
    - prescribed from the observed upward longwave radiation,
for the runs made by main_SEB_firn_AWS.py.

Outputs (in <output_path>/comparison/):
    melt_comparison_summary.csv   seasonal totals and hourly statistics
    melt_tolerance_sensitivity.csv  observed-Tsurf runs with other melt
                                  tolerances (if they were run)
    <station>_melt_comparison.png cumulative/daily melt, surface temperature,
                                  surface height vs observations, melt contingency
    <station>_seb_difference.png  daily mean difference of the SEB components
"""
import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

output_path = './output/AWS_L3/'
forcing_folder = './input/weather data/AWS_L3/'
modes = ['inverted', 'observed']
colors = {'inverted': 'tab:blue', 'observed': 'tab:red'}
melt_season = ('2026-05-01', '2026-09-30 23:00')


def load_run(station, mode):
    with xr.open_dataset(f'{output_path}/{station}_Tsurf_{mode}/{station}_surface.nc') as ds:
        df = ds.to_dataframe()
    df['Tsurf_C'] = df.Tsurf - 273.15
    df['Tsurf_obs_C'] = df.Tsurf_obs - 273.15
    df['SRnet'] = df.SRin - df.SRout
    df['LRnet'] = df.LRin - df.LRout_mdl
    return df


def summary_stats(station, runs, df_obs):
    inv, obs = runs['inverted'], runs['observed']
    sl = slice(*melt_season)
    m_inv = inv.melt_mweq.loc[sl] * 1000  # mm w.e. per hour
    m_obs = obs.melt_mweq.loc[sl] * 1000
    melt_inv, melt_obs = m_inv > 0.01, m_obs > 0.01
    has_obs_T = obs.Tsurf_obs.loc[sl].notnull()
    out = {
        'station': station,
        'period': f'{melt_season[0]} to {melt_season[1][:10]}',
        'melt_inverted_mmwe': m_inv.sum(),
        'melt_observed_mmwe': m_obs.sum(),
        'melt_diff_obs_minus_inv_mmwe': m_obs.sum() - m_inv.sum(),
        'melt_diff_percent': 100 * (m_obs.sum() - m_inv.sum()) / m_inv.sum(),
        'runoff_inverted_mmwe': inv.runoff_mweq.loc[sl].sum() * 1000,
        'runoff_observed_mmwe': obs.runoff_mweq.loc[sl].sum() * 1000,
        'refreezing_inverted_mmwe': inv.refreezing_mweq.loc[sl].sum() * 1000,
        'refreezing_observed_mmwe': obs.refreezing_mweq.loc[sl].sum() * 1000,
        'sublimation_inverted_mmwe': inv.sublimation_mweq.loc[sl].sum() * 1000,
        'sublimation_observed_mmwe': obs.sublimation_mweq.loc[sl].sum() * 1000,
        'melt_hours_inverted': int(melt_inv.sum()),
        'melt_hours_observed': int(melt_obs.sum()),
        'melt_hours_both': int((melt_inv & melt_obs).sum()),
        'melt_hours_inverted_only': int((melt_inv & ~melt_obs).sum()),
        'melt_hours_observed_only': int((~melt_inv & melt_obs).sum()),
        'melt_inverted_during_inverted_only_hours_mmwe': m_inv[melt_inv & ~melt_obs].sum(),
        'melt_observed_during_observed_only_hours_mmwe': m_obs[~melt_inv & melt_obs].sum(),
        'hours_with_observed_Tsurf_percent': 100 * has_obs_T.mean(),
        # surface temperature of the inverted run vs observed (where observed)
        'Tsurf_inv_minus_obs_bias_K': (inv.Tsurf - obs.Tsurf_obs).loc[sl].mean(),
        'Tsurf_inv_minus_obs_rmse_K': np.sqrt(((inv.Tsurf - obs.Tsurf_obs).loc[sl] ** 2).mean()),
    }
    # independent check: modelled vs observed surface lowering over the melt season
    z = df_obs.z_surf_combined.loc[sl].dropna()
    if len(z) > 48:
        z_start = z.iloc[:24].mean()
        z_end = z.iloc[-24:].mean()
        out['obs_surface_height_change_m'] = z_end - z_start
        for mode in modes:
            h = runs[mode].surface_height.loc[z.index]
            out[f'mdl_surface_height_change_{mode}_m'] = h.iloc[-24:].mean() - h.iloc[:24].mean()
    return out


def plot_station(station, runs, df_obs):
    inv, obs = runs['inverted'], runs['observed']
    sl = slice(*melt_season)
    fig, axs = plt.subplots(5, 1, figsize=(12, 16), sharex=True)

    ax = axs[0]
    for mode in modes:
        m = runs[mode].melt_mweq.loc[sl].cumsum() * 1000
        ax.plot(m, color=colors[mode], label=f'{mode} Tsurf: {m.iloc[-1]:.0f} mm w.e.')
    ax.set_ylabel('cumulative melt\n(mm w.e.)')
    ax.legend(loc='upper left')
    ax.set_title(f'{station}: melt from inverted vs observed surface temperature')

    ax = axs[1]
    for mode in modes:
        ax.plot(runs[mode].melt_mweq.loc[sl].resample('D').sum() * 1000,
                color=colors[mode], drawstyle='steps-mid', label=mode)
    ax.set_ylabel('daily melt\n(mm w.e.)')
    ax.legend(loc='upper left')

    ax = axs[2]
    ax.plot(inv.Tsurf_C.loc[sl].resample('D').mean(), color=colors['inverted'], label='inverted (model)')
    ax.plot(obs.Tsurf_obs_C.loc[sl].resample('D').mean(), color=colors['observed'], label='observed (LRout)')
    ax.plot(inv.Tsurf_C.loc[sl].resample('D').max(), color=colors['inverted'], ls=':', lw=0.8)
    ax.plot(obs.Tsurf_obs_C.loc[sl].resample('D').max(), color=colors['observed'], ls=':', lw=0.8)
    ax.set_ylabel('surface temperature\n(°C, daily mean; dotted: max)')
    ax.legend(loc='lower left')

    ax = axs[3]
    # hourly melt contingency, daily counts
    m_inv = inv.melt_mweq.loc[sl] > 1e-5
    m_obs = obs.melt_mweq.loc[sl] > 1e-5
    cnt = pd.DataFrame({
        'both': (m_inv & m_obs), 'inverted only': (m_inv & ~m_obs),
        'observed only': (~m_inv & m_obs)}).resample('D').sum()
    bottom = np.zeros(len(cnt))
    for col, c in zip(cnt.columns, ['tab:purple', colors['inverted'], colors['observed']]):
        ax.bar(cnt.index, cnt[col], bottom=bottom, color=c, width=1, label=col)
        bottom += cnt[col].values
    ax.set_ylabel('melt hours\nper day')
    ax.legend(loc='upper left', ncol=3)

    ax = axs[4]
    z = df_obs.z_surf_combined.loc[sl]
    if z.notnull().sum() > 48:
        t0 = z.first_valid_index()
        ax.plot(z - z.loc[t0:].iloc[:24].mean(), color='k', label='observed (z_surf_combined)')
        for mode in modes:
            h = runs[mode].surface_height.loc[sl]
            ax.plot(h - h.loc[t0:].iloc[:24].mean(), color=colors[mode], label=f'model, {mode} Tsurf')
    ax.set_ylabel('surface height\nchange (m)')
    ax.legend(loc='lower left')
    for ax in axs:
        ax.grid(alpha=0.4)
    fig.tight_layout()
    fig.savefig(f'{output_path}/comparison/{station}_melt_comparison.png', dpi=130)
    plt.close(fig)

    # SEB component differences and Tsurf scatter
    fig, axs = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={'width_ratios': [2.5, 1]})
    ax = axs[0]
    for var, c in zip(['LRnet', 'SHF', 'LHF', 'GF', 'meltflux'],
                      ['tab:orange', 'tab:red', 'tab:green', 'tab:brown', 'k']):
        d = (obs[var] - inv[var]).loc[sl].resample('D').mean()
        ax.plot(d, color=c, label=var)
    ax.axhline(0, color='grey', lw=0.8)
    ax.set_ylabel('observed-Tsurf minus inverted-Tsurf run\n(daily mean, W m-2)')
    ax.legend(ncol=5, loc='upper left')
    ax.grid(alpha=0.4)
    ax.set_title(station)

    ax = axs[1]
    both = pd.concat([inv.Tsurf_C, obs.Tsurf_obs_C], axis=1, keys=['inv', 'obs']).loc[sl].dropna()
    ax.hexbin(both.obs, both.inv, gridsize=50, bins='log', cmap='viridis', mincnt=1)
    lim = [min(both.min()), 0.5]
    ax.plot(lim, lim, 'k--', lw=0.8)
    ax.set_xlabel('observed Tsurf (°C)')
    ax.set_ylabel('inverted Tsurf (°C)')
    ax.set_title('hourly surface temperature, melt season')
    fig.tight_layout()
    fig.savefig(f'{output_path}/comparison/{station}_seb_difference.png', dpi=130)
    plt.close(fig)


def tolerance_sensitivity(station):
    # all observed-Tsurf runs: default (0.5 K) and the ones named *_tol<value>
    rows = []
    sl = slice(*melt_season)
    for folder in sorted(glob.glob(f'{output_path}/{station}_Tsurf_observed*')):
        name = os.path.basename(folder)
        if '_GF' in name:
            continue  # other model configurations
        tol = float(name.split('_tol')[1]) if '_tol' in name else 0.5
        with xr.open_dataset(f'{folder}/{station}_surface.nc') as ds:
            df = ds.to_dataframe().loc[sl]
        rows.append({'station': station, 'Tsurf_obs_melt_tol_K': tol,
                     'melt_mmwe': df.melt_mweq.sum() * 1000,
                     'refreezing_mmwe': df.refreezing_mweq.sum() * 1000,
                     'runoff_mmwe': df.runoff_mweq.sum() * 1000,
                     'melt_hours': int((df.melt_mweq > 1e-5).sum())})
    inv = load_run(station, 'inverted').loc[sl]
    rows.append({'station': station, 'Tsurf_obs_melt_tol_K': np.nan,
                 'melt_mmwe': inv.melt_mweq.sum() * 1000,
                 'refreezing_mmwe': inv.refreezing_mweq.sum() * 1000,
                 'runoff_mmwe': inv.runoff_mweq.sum() * 1000,
                 'melt_hours': int((inv.melt_mweq > 1e-5).sum())})
    df = pd.DataFrame(rows)
    df['run'] = np.where(df.Tsurf_obs_melt_tol_K.isnull(), 'inverted', 'observed')
    return df.sort_values(['run', 'Tsurf_obs_melt_tol_K'])


def main(station_list=('TAS_A', 'KAN_U', 'KAN_M', 'QAS_M', 'QAS_U')):
    os.makedirs(f'{output_path}/comparison', exist_ok=True)
    summary = []
    for station in station_list:
        runs = {mode: load_run(station, mode) for mode in modes}
        df_obs = pd.read_csv(forcing_folder + station + '_forcing.csv',
                             index_col=0, parse_dates=True)
        summary.append(summary_stats(station, runs, df_obs))
        plot_station(station, runs, df_obs)
    summary = pd.DataFrame(summary).set_index('station').T
    summary.to_csv(f'{output_path}/comparison/melt_comparison_summary.csv')
    print(summary.to_string())

    sens = pd.concat([tolerance_sensitivity(s) for s in station_list])
    sens.to_csv(f'{output_path}/comparison/melt_tolerance_sensitivity.csv',
                index=False, float_format='%.1f')
    print(sens.to_string(index=False, float_format='%.1f'))
    return summary


if __name__ == "__main__":
    main()
