"""Make the initial density and temperature profiles of DYE-2 (DY2) from
observations:

input/initial state/DY2_initial_density_bulk.csv  (depth,density_bulk)
    1 cm resolution, depth of the bottom of each cm (m), density (kg m-3):
    - 0 to 0.98 m: GEUS snow pit at DYE-2, 2022-06-15
    - 0.98 to 9.56 m: GEUS firn core (thermistor hole) at DYE-2, 2022-06-15.
      Samples flagged in dens_note (bad quality, broken pieces) are skipped
      and gaps are filled by linear interpolation.
    - below: the previous DY2 profile.
    Snow pit and core from Vandecrux et al. (2023), GEUS snow and firn data
    in Greenland, https://doi.org/10.22008/FK2/9QEOWZ

input/initial state/DY2_initial_T_ice.csv  (depth,T_ice)
    depth (m), temperature (degC): thermistor string at DYE-2 on
    2019-01-01 (Covi et al., 2022, doi:10.18739/A2BN9X444), as compiled in
    SUMup (data_formatted.csv).

Usage:
    python tools/make_DY2_initial_profiles.py --pits <snow_pit_sheets/csv folder>
        --covi <SUMup .../temperature data/Covi/data_formatted.csv>
"""
import argparse
import os

import numpy as np
import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FOLDER = os.path.join(REPO, "input", "initial state")
PIT = "2022_GC-Net_GEUS_DYE-2_Snow_Pit.csv"
CORE = "2022_GC-Net_GEUS_DYE-2_Therm_Core.csv"
TEMPERATURE_DATE = "2019-01-01"


def read_density(path):
    """Density samples (top and bottom depth in m, density in kg m-3)."""
    df = pd.read_csv(path, comment="#")
    df = df.loc[df.dens_note.isna()]
    rho = df[["dens_a", "dens_b"]].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    out = pd.DataFrame({"top": df.depth_top_dens / 100, "bottom": df.depth_bot_dens / 100,
                        "density": rho}).dropna()
    return out.loc[out.bottom > out.top].sort_values("top")


def on_cm_grid(samples, depth):
    """Density of the samples at the cm whose bottom is at depth; NaN outside."""
    out = np.full(len(depth), np.nan)
    for _, s in samples.iterrows():
        out[(depth > s.top + 1e-9) & (depth <= s.bottom + 1e-9)] = s.density
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pits", required=True, help="folder of snow pit csv files")
    parser.add_argument("--covi", required=True, help="Covi data_formatted.csv")
    args = parser.parse_args()

    # density
    density_file = os.path.join(FOLDER, "DY2_initial_density_bulk.csv")
    previous = pd.read_csv(density_file)
    depth = previous.depth.values
    pit = read_density(os.path.join(args.pits, PIT))
    core = read_density(os.path.join(args.pits, CORE))
    pit_bottom = pit.bottom.max()
    core = core.loc[core.bottom > pit_bottom]
    core.loc[core.top < pit_bottom, "top"] = pit_bottom
    rho = on_cm_grid(pit, depth)
    core_rho = on_cm_grid(core, depth)
    rho = np.where(np.isnan(rho), core_rho, rho)
    observed = depth <= core.bottom.max() + 1e-9
    # gaps between samples: linear interpolation
    rho[observed] = pd.Series(rho[observed]).interpolate(limit_direction="both").values
    rho[~observed] = previous.density_bulk.values[~observed]
    pd.DataFrame({"depth": depth, "density_bulk": np.round(rho, 2)}).to_csv(
        density_file, index=False)
    print(f"{density_file}: observed to {core.bottom.max():.2f} m, "
          f"mean density top 5 m {np.nanmean(rho[depth <= 5]):.0f} kg m-3")

    # temperature
    covi = pd.read_csv(args.covi, usecols=["timestamp", "name", "depth", "temperature"])
    profile = covi.loc[(covi.name == "DYE-2") & (covi.timestamp == TEMPERATURE_DATE)]
    profile = profile.dropna().sort_values("depth")
    temperature_file = os.path.join(FOLDER, "DY2_initial_T_ice.csv")
    pd.DataFrame({"depth": profile.depth.round(3).values,
                  "T_ice": profile.temperature.values}).to_csv(temperature_file, index=False)
    print(f"{temperature_file}: {len(profile)} depths, {profile.depth.min():.2f} to "
          f"{profile.depth.max():.2f} m, {profile.temperature.min():.1f} to "
          f"{profile.temperature.max():.1f} degC")


if __name__ == "__main__":
    main()
