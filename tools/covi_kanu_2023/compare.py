# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

Compares, at KAN_U for May-September 2023, the melt and surface energy
budget of the GEUS SEB-firn model and of the model of Covi et al. (2022),
each run with an inverted (skin layer) and an observed surface temperature,
on the same forcing and initial firn state (see prepare_inputs.py).

Outputs in ./output/AWS_L3_KAN_U_2023/comparison/
"""
import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

geus_path = './output/AWS_L3_KAN_U_2023/'
covi_path = '../Covi_1D_SEB_model/output/'
forcing_file = './input/weather data/AWS_L3_KAN_U_2023/KAN_U_forcing.csv'
out_path = geus_path + 'comparison/'
sigma = 5.67e-8

runs = {
    'Covi skin layer': dict(model='covi', mode='skin', color='tab:red', ls='-'),
    'Covi observed Tsurf': dict(model='covi', mode='observed', color='tab:red', ls='--'),
    'GEUS inverted': dict(model='geus', mode='inverted', color='tab:blue', ls='-'),
    'GEUS observed Tsurf': dict(model='geus', mode='observed', color='tab:blue', ls='--'),
}


def load_covi(mode):
    df = pd.read_csv(f'{covi_path}/surface_KAN_U2023_{mode}.txt', index_col=0, parse_dates=True)
    # 4-min substeps -> hourly: fluxes averaged, melt summed (mm w.e.)
    flux = df[['Tsurf', 'SWnet', 'LWnet', 'QH', 'QL', 'QG', 'QR', 'Q', 'Qres']].resample('h').mean()
    flux['melt_mm'] = df.surfmelt.resample('h').sum()
    flux['sh'] = df.sh.resample('h').last()
    return flux


def load_geus(mode):
    with xr.open_dataset(f'{geus_path}/KAN_U_Tsurf_{mode}/KAN_U_surface.nc') as ds:
        df = ds.to_dataframe()
    with xr.open_dataset(f'{geus_path}/KAN_U_Tsurf_{mode}/KAN_U_density_bulk.nc') as ds:
        rho_top = ds.density_bulk.isel(level=0).to_series()
    out = pd.DataFrame(index=df.index)
    out['Tsurf'] = df.Tsurf - 273.15
    out['SWnet'] = df.SRin - df.SRout
    out['LWnet'] = df.LRin - df.LRout_mdl
    out['QH'] = df.SHF
    out['QL'] = df.LHF
    out['Q'] = df.meltflux
    out['melt_mm'] = df.melt_mweq * 1000
    # SEB ground heat flux as closure term (exact when the energy budget is
    # closed, i.e. in the inverted run)
    out['QG'] = out.Q - (out.SWnet + out.LWnet + out.QH + out.QL)
    # surface height computed as in Covi et al.: snowfall / 350 - melt / top layer density
    out['sh'] = (df.snowfall_mweq * 1000 / 350 - out.melt_mm / rho_top).cumsum()
    return out


def main():
    os.makedirs(out_path, exist_ok=True)
    forcing = pd.read_csv(forcing_file, index_col=0, parse_dates=True)
    # observed surface temperature with emissivity 1, as in Covi et al.
    Ts_obs = np.minimum((forcing.LongwaveRadiationUpWm2 / sigma) ** 0.25 - 273.15, 0)
    data = {name: (load_covi(r['mode']) if r['model'] == 'covi' else load_geus(r['mode']))
            for name, r in runs.items()}

    jja = slice('2023-06-01', '2023-08-31 23:00')
    rows = []
    for name, d in data.items():
        dT = (d.Tsurf - Ts_obs).loc[jja]
        rows.append({
            'run': name,
            'melt_May_Sep_mmwe': d.melt_mm.sum(),
            'melt_hours': int((d.melt_mm > 0.01).sum()),
            'Tsurf_bias_JJA_K': dT.mean(),
            'Tsurf_rmse_JJA_K': np.sqrt((dT ** 2).mean()),
            'hours_Tsurf_0_while_obs_below_-1C': int(((d.Tsurf >= -0.01) & (Ts_obs < -1)).loc[jja].sum()),
            'QH_JJA': d.QH.loc[jja].mean(), 'QL_JJA': d.QL.loc[jja].mean(),
            'QG_JJA': d.QG.loc[jja].mean(), 'LWnet_JJA': d.LWnet.loc[jja].mean(),
            'Qmelt_JJA': d.Q.loc[jja].clip(lower=0).mean(),
        })
    summary = pd.DataFrame(rows).set_index('run')
    summary.to_csv(out_path + 'summary.csv', float_format='%.2f')
    print(summary.round(2).to_string())

    # ---- figure
    fig, axs = plt.subplots(4, 1, figsize=(12, 15), sharex=True)
    for name, d in data.items():
        r = runs[name]
        axs[0].plot(d.melt_mm.cumsum(), color=r['color'], ls=r['ls'],
                    label=f'{name}: {d.melt_mm.sum():.0f} mm w.e.')
    axs[0].set_ylabel('cumulative surface melt\n(mm w.e.)')
    axs[0].legend(loc='upper left')
    axs[0].set_title('KAN_U 2023: GEUS model vs Covi et al. (2022) model, same forcing and initial state')

    ax = axs[1]
    ax.plot(Ts_obs.resample('D').mean(), color='k', lw=2, label='observed (LWout, emissivity 1)')
    for name in ['Covi skin layer', 'GEUS inverted']:
        r = runs[name]
        ax.plot(data[name].Tsurf.resample('D').mean(), color=r['color'], label=name)
    ax.set_ylabel('daily mean surface\ntemperature (°C)')
    ax.legend(loc='lower center', ncol=3)

    ax = axs[2]
    z = forcing.z_surf_combined
    z = z - z.iloc[:24].mean()
    ax.plot(z, color='k', lw=1, label='SR50')
    for name, d in data.items():
        r = runs[name]
        ax.plot(d.sh - d.sh.iloc[0], color=r['color'], ls=r['ls'], label=name)
    ax.set_ylabel('relative surface\nheight (m)')
    ax.legend(loc='lower left', ncol=2)

    ax = axs[3]
    for name in ['Covi skin layer', 'GEUS inverted']:
        r = runs[name]
        ax.plot(data[name].QG.rolling('15D', center=True).mean(), color=r['color'],
                label=f'QG {name}')
        ax.plot(data[name].QH.rolling('15D', center=True).mean(), color=r['color'], ls=':',
                label=f'QH {name}')
    ax.axhline(0, color='grey', lw=0.8)
    ax.set_ylabel('15-day mean flux\n(W m-2, + towards surface)')
    ax.legend(loc='upper left', ncol=2)
    for ax in axs:
        ax.grid(alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path + 'KAN_U_2023_GEUS_vs_Covi.png', dpi=130)
    plt.close(fig)

    # zoom on a 10-day period of July, as Fig. 12 of Covi et al.
    fig, ax = plt.subplots(figsize=(12, 4))
    sl = slice('2023-07-07', '2023-07-17')
    ax.plot(forcing.AirTemperatureC.loc[sl], color='tab:green', ls='--', label='air temperature')
    ax.plot(Ts_obs.loc[sl], color='k', lw=2, label='observed Tsurf')
    for name in ['Covi skin layer', 'GEUS inverted']:
        ax.plot(data[name].Tsurf.loc[sl], color=runs[name]['color'], label=name)
    ax.set_ylabel('temperature (°C)')
    ax.legend(ncol=4)
    ax.grid(alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path + 'KAN_U_2023_Tsurf_July.png', dpi=130)
    plt.close(fig)
    return summary


if __name__ == "__main__":
    main()
