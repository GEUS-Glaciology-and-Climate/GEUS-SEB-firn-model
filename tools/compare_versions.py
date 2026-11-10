# -*- coding: utf-8 -*-
"""
Created on %(date)s
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
    import matplotlib.pyplot as plt
"""
import argparse
import os
import re

# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

# default comparison
DEFAULT_NEW = "output/new version"
DEFAULT_BASE = "output/HH precipitation"
DEFAULT_STATIONS = ["all"]
# 2D variables given per layer in m w.e. (or m): shown per m of layer
# (label, unit, factor)
PER_LAYER_AMOUNT = {"slwc": ("Liquid water", "kg m⁻³", 1000),
                    "snic": ("Ice", "kg m⁻³", 1000),
                    "snowc": ("Snow", "kg m⁻³", 1000),
                    "rfrz": ("Refreezing", "kg m⁻³ per step", 1000),
                    "compaction": ("Compaction", "m m⁻¹ per step", 1)}
UNITS_2D = {"T_ice": ("Temperature", "K"), "rhofirn": ("Firn density", "kg m⁻³"),
            "density_bulk": ("Bulk density", "kg m⁻³"), "dgrain": ("Grain size", "mm")}

RUN_FOLDER = re.compile(r"^(?P<station>.+)_(?P<num_lay>\d+)_layers_(?P<freq>[^_]+)$")
ANNUAL_VARS = {  # variable: label
    "melt_mweq": "Melt",
    "runoff_mweq": "Runoff",
    "refreezing_mweq": "Refreezing",
    "sublimation_mweq": "Sublimation (+ = deposition)",
    "snowfall_mweq": "Snowfall",
    "rainfall_mweq": "Rainfall",
    "smb": "SMB",
}
DEPTHS_T = [1, 5, 10]  # m, firn temperature depths
RHO_WATER, RHO_ICE = 999.8395, 900  # as in the model


def list_runs(folder):
    """{station: run folder path} of the runs in a version folder."""
    runs = {}
    for name in sorted(os.listdir(folder)):
        m = RUN_FOLDER.match(name)
        path = os.path.join(folder, name)
        if m and os.path.isfile(os.path.join(path, f"{m['station']}_surface.nc")):
            runs[m["station"]] = path
    return runs


def emissivity(run):
    """Surface emissivity used in the run (constants.csv), default 0.98."""
    try:
        c = pd.read_csv(os.path.join(run, "constants.csv"), index_col=0).iloc[:, 0]
        return float(c["em"])
    except Exception:
        return 0.98


def load_surface(run, station):
    with xr.open_dataset(os.path.join(run, f"{station}_surface.nc")) as ds:
        df = ds.to_dataframe()
    df = df[~df.index.duplicated(keep="first")]
    # SMB with the same convention for both versions (see module docstring)
    df["smb"] = (df.snowfall_mweq + df.rainfall_mweq - df.runoff_mweq
                 + df.sublimation_mweq)
    # surface temperature from the modelled outgoing longwave radiation
    em = emissivity(run)
    df["Tsurf"] = ((df.LRout_mdl - (1 - em) * df.LRin) / em / 5.67e-8) ** 0.25 - 273.15
    return df


def annual_sums(df):
    """Sums per calendar year, keeping only years covered at >= 95 %."""
    steps = df.index.to_series().diff().median()
    full = pd.Timedelta("365D") / steps
    g = df[list(ANNUAL_VARS)].groupby(df.index.year)
    out = g.sum()
    return out[g.size() >= 0.95 * full]


def load_subsurface(run, station, time, step=8):
    """Firn temperature at DEPTHS_T, mean bulk density of the top 10 m and
    total liquid water, every `step`-th time step of `time`."""
    def read(var):
        with xr.open_dataset(os.path.join(run, f"{station}_{var}.nc")) as ds:
            da = ds[var].sel(time=time[::step])
            return da.transpose("time", "level").values.astype(float)
    try:
        snowc, snic, rhofirn = read("snowc"), read("snic"), read("rhofirn")
        T, slwc = read("T_ice") - 273.15, read("slwc")
    except (FileNotFoundError, KeyError):
        return None
    # physical layer thickness (m) and mid-depth, without liquid water
    thick = snowc * RHO_WATER / np.where(rhofirn > 0, rhofirn, RHO_ICE) + snic * RHO_WATER / RHO_ICE
    bottom = np.cumsum(thick, axis=1)
    mid = bottom - thick / 2
    out = pd.DataFrame(index=time[::step])
    for d in DEPTHS_T:
        out[f"T_{d}m"] = [np.interp(d, mid[i], T[i], left=np.nan, right=np.nan)
                          for i in range(len(out))]
    # mass (kg m-2) of snow and ice in the top 10 m
    frac = np.clip((10 - (bottom - thick)) / np.where(thick > 0, thick, 1), 0, 1)
    mass = ((snowc + snic) * RHO_WATER * frac).sum(axis=1)
    out["rho_10m"] = np.where(bottom[:, -1] >= 10, mass / 10, np.nan)
    out["lwc"] = slwc.sum(axis=1)
    return out


def plot_surface(station, dfs, names, out):
    fig, ax = plt.subplots(4, 1, figsize=(10, 13), sharex=False)
    colors = {"melt_mweq": "tab:red", "runoff_mweq": "tab:blue", "refreezing_mweq": "tab:green"}
    styles = ["-", "--"]

    for df, name, ls in zip(dfs, names, styles):
        ax[0].plot(df.index, df.smb.cumsum(), ls=ls, label=name)
    ax[0].set_ylabel("Cumulative SMB\n(m w.e.)")
    ax[0].legend()

    ann = [annual_sums(df) for df in dfs]
    for var, col in colors.items():
        for a, name, ls in zip(ann, names, styles):
            ax[1].plot(a.index, a[var], ls=ls, color=col, marker=".",
                       label=f"{ANNUAL_VARS[var]}, {name}")
    ax[1].set_ylabel("Annual total\n(m w.e.)")
    ax[1].legend(fontsize=8, ncol=3)

    diff = (ann[0] - ann[1]).dropna()
    w = 0.4
    ax[2].bar(diff.index - w / 2, diff.smb, w, label="SMB")
    ax[2].bar(diff.index + w / 2, diff.sublimation_mweq, w, label="Sublimation")
    ax[2].axhline(0, color="k", lw=0.5)
    ax[2].set_ylabel(f"Annual difference\n{names[0]} − {names[1]}\n(m w.e.)")
    ax[2].legend(fontsize=8)

    ts = [df.Tsurf.resample("MS").mean() for df in dfs]
    ax[3].plot(ts[0].index, ts[0] - ts[1], color="k")
    ax[3].axhline(0, color="k", lw=0.5)
    ax[3].set_ylabel(f"Monthly surface T\n{names[0]} − {names[1]} (K)")

    fig.suptitle(station)
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{station}_surface.png"), dpi=120)
    plt.close(fig)


def plot_subsurface(station, subs, names, out):
    fig, ax = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    styles = ["-", "--"]
    for sub, name, ls in zip(subs, names, styles):
        for d, col in zip(DEPTHS_T, ["tab:orange", "tab:purple", "tab:blue"]):
            ax[0].plot(sub.index, sub[f"T_{d}m"], ls=ls, color=col, lw=0.8,
                       label=f"{d} m, {name}")
        ax[1].plot(sub.index, sub.rho_10m, ls=ls, label=name)
        ax[2].plot(sub.index, sub.lwc, ls=ls, lw=0.8, label=name)
    ax[0].set_ylabel("Firn temperature (°C)")
    ax[0].legend(fontsize=8, ncol=3)
    ax[1].set_ylabel("Mean density of\nthe top 10 m (kg m⁻³)")
    ax[1].legend(fontsize=8)
    ax[2].set_ylabel("Total liquid water\n(m w.e.)")
    fig.suptitle(station)
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{station}_subsurface.png"), dpi=120)
    plt.close(fig)


def plot_summary(summary, names, out):
    variables = ["melt_mweq", "runoff_mweq", "refreezing_mweq",
                 "sublimation_mweq", "smb", "T_10m"]
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    for ax, var in zip(axes.flat, variables):
        x, y = summary[(var, names[1])], summary[(var, names[0])]
        ax.scatter(x, y, s=12)
        for st in summary.index:
            ax.annotate(st, (x[st], y[st]), fontsize=6)
        lim = [np.nanmin([x.min(), y.min()]), np.nanmax([x.max(), y.max()])]
        ax.plot(lim, lim, color="k", lw=0.5)
        ax.set_xlabel(names[1])
        ax.set_ylabel(names[0])
        unit = "°C" if var == "T_10m" else "m w.e. yr⁻¹"
        ax.set_title(f"{ANNUAL_VARS.get(var, 'Firn temperature at 10 m')} ({unit})")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "summary_scatter.png"), dpi=120)
    plt.close(fig)


def is_amount(var):
    """Variables summed per month (amounts per time step)."""
    return var.endswith("mweq") or var == "smb"


def plot_all_surface_variables(station, dfs, names, out):
    """One panel per variable present in both surface files."""
    variables = [v for v in dfs[0].columns if v in dfs[1].columns
                 and pd.api.types.is_numeric_dtype(dfs[0][v])]
    ncol = 3
    nrow = int(np.ceil(len(variables) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(16, 2.3 * nrow), sharex=True,
                             squeeze=False)
    for ax, var in zip(axes.flat, variables):
        for df, name, ls in zip(dfs, names, ["-", "--"]):
            monthly = df[var].resample("MS")
            monthly = monthly.sum() if is_amount(var) else monthly.mean()
            ax.plot(monthly.index, monthly, ls=ls, lw=0.8, label=name)
        ax.set_title(f"{var} (monthly {'sum' if is_amount(var) else 'mean'})", fontsize=9)
        ax.tick_params(labelsize=8)
    for ax in list(axes.flat)[len(variables):]:
        ax.set_visible(False)
    axes.flat[0].legend(fontsize=8)
    fig.suptitle(station)
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{station}_all_surface_variables.png"), dpi=110)
    plt.close(fig)


def list_2d_variables(run, station):
    """2D variables written for a run: <station>_<var>.nc with a variable
    <var> on (level, time) and the layer depth."""
    out = []
    for f in os.listdir(run):
        if not (f.startswith(station + "_") and f.endswith(".nc")):
            continue
        var = f[len(station) + 1:-3]
        if var == "surface":
            continue
        try:
            with xr.open_dataset(os.path.join(run, f)) as ds:
                if var in ds and set(ds[var].dims) == {"level", "time"} and "depth" in ds:
                    out.append(var)
        except Exception:
            pass
    return sorted(out)


def section_on_depth_grid(run, station, var, time, grid, step=8):
    """Values of a 2D variable on a regular depth grid (m), every `step`-th
    time step: each grid point takes the value of the layer it falls in.
    Amounts per layer (PER_LAYER_AMOUNT) are divided by the layer thickness."""
    with xr.open_dataset(os.path.join(run, f"{station}_{var}.nc")) as ds:
        ds = ds.sel(time=time[::step])
        values = ds[var].transpose("level", "time").values.astype(float)
        bottom = ds.depth.transpose("level", "time").values.astype(float)  # layer bottom (m)
    if var in PER_LAYER_AMOUNT:
        top = np.vstack([np.zeros((1, bottom.shape[1])), bottom[:-1]])
        thick = bottom - top
        values = PER_LAYER_AMOUNT[var][2] * values / np.where(thick > 0, thick, np.nan)
    out = np.full((len(grid), bottom.shape[1]), np.nan)
    for t in range(bottom.shape[1]):
        idx = np.searchsorted(bottom[:, t], grid)
        valid = idx < bottom.shape[0]
        out[valid, t] = values[idx[valid], t]
    return out


def plot_2d_variable(station, var, runs, names, time, out, max_depth, step=8):
    """Depth-time sections of both versions and their difference."""
    grid = np.linspace(0, max_depth, 201)[1:] - max_depth / 400
    sections = [section_on_depth_grid(r, station, var, time, grid, step) for r in runs]
    if all(np.all(np.isnan(s)) for s in sections):
        return
    t = time[::step]
    label, unit = (PER_LAYER_AMOUNT[var][:2] if var in PER_LAYER_AMOUNT
                   else UNITS_2D.get(var, (var, "")))
    vmin, vmax = np.nanpercentile(np.concatenate([s.ravel() for s in sections]), [1, 99])
    diff = sections[0] - sections[1]
    dmax = np.nanpercentile(np.abs(diff), 99)
    dmax = dmax if dmax > 0 else 1e-12

    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
    for ax, sec, title in zip(axes[:2], sections, names):
        im = ax.pcolormesh(t, grid, sec, vmin=vmin, vmax=vmax, cmap="viridis",
                           shading="auto")
        ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, label=f"{label} ({unit})")
    im = axes[2].pcolormesh(t, grid, diff, vmin=-dmax, vmax=dmax, cmap="RdBu_r",
                            shading="auto")
    axes[2].set_title(f"{names[0]} − {names[1]}", fontsize=10)
    fig.colorbar(im, ax=axes[2], label=f"Δ {label} ({unit})")
    for ax in axes:
        ax.set_ylim(max_depth, 0)
        ax.set_ylabel("Depth (m)")
    fig.suptitle(f"{station}: {var}")
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{station}_2D_{var}.png"), dpi=110)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--new", default=DEFAULT_NEW,
                   help=f"folder of the new version runs (default: {DEFAULT_NEW})")
    p.add_argument("--base", default=DEFAULT_BASE,
                   help=f"folder of the reference runs (default: {DEFAULT_BASE})")
    p.add_argument("--names", nargs=2, default=None,
                   help="labels of the new and base versions (default: folder names)")
    p.add_argument("--stations", nargs="*", default=DEFAULT_STATIONS,
                   help="stations to compare, or 'all' for all stations in both "
                        f"folders (default: {' '.join(DEFAULT_STATIONS)})")
    p.add_argument("--max-depth", type=float, default=20,
                   help="depth of the 2D sections (m, default 20)")
    p.add_argument("--out", help="output folder (default: output/comparison <new> vs <base>)")
    p.add_argument("--no-subsurface", action="store_true",
                   help="skip the subsurface figures (faster)")
    args = p.parse_args()
    # folders are relative to the repository root, wherever the script is run from
    os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

    names = args.names or [os.path.basename(os.path.normpath(f)) for f in (args.new, args.base)]
    out = args.out or os.path.join("output", f"comparison {names[0]} vs {names[1]}")
    os.makedirs(out, exist_ok=True)

    runs_new, runs_base = list_runs(args.new), list_runs(args.base)
    if not args.stations or [s.lower() for s in args.stations] == ["all"]:
        stations = sorted(set(runs_new) & set(runs_base))
    else:
        stations = args.stations
    missing = sorted(set(runs_new) ^ set(runs_base))
    if missing and stations == sorted(set(runs_new) & set(runs_base)):
        print("Only in one version (skipped):", ", ".join(missing))

    rows = {}
    for station in stations:
        if station not in runs_new or station not in runs_base:
            print(f"{station}: not in both versions, skipped")
            continue
        try:
            dfs = [load_surface(runs_new[station], station),
                   load_surface(runs_base[station], station)]
        except Exception as e:
            print(f"{station}: cannot read surface output ({e}), skipped")
            continue
        common = dfs[0].index.intersection(dfs[1].index)
        if len(common) == 0:
            print(f"{station}: no common period, skipped")
            continue
        dfs = [df.loc[common] for df in dfs]
        print(f"{station}: {common[0]:%Y-%m-%d} to {common[-1]:%Y-%m-%d}")
        plot_surface(station, dfs, names, out)
        plot_all_surface_variables(station, dfs, names, out)

        row = {}
        for df, name in zip(dfs, names):
            mean = annual_sums(df).mean()
            for var in ANNUAL_VARS:
                row[(var, name)] = mean[var]
        if not args.no_subsurface:
            subs = [load_subsurface(runs_new[station], station, common),
                    load_subsurface(runs_base[station], station, common)]
            if all(s is not None for s in subs):
                plot_subsurface(station, subs, names, out)
                for sub, name in zip(subs, names):
                    row[("T_10m", name)] = sub.T_10m.mean()
            runs = [runs_new[station], runs_base[station]]
            common_2d = sorted(set(list_2d_variables(runs[0], station))
                               & set(list_2d_variables(runs[1], station)))
            for var in common_2d:
                try:
                    plot_2d_variable(station, var, runs, names, common, out, args.max_depth)
                except Exception as e:
                    print(f"{station}: {var} not plotted ({e})")
            print(f"{station}: 2D variables:", ", ".join(common_2d))
        rows[station] = row

    if not rows:
        print("Nothing to compare")
        return
    summary = pd.DataFrame(rows).T
    summary.columns = pd.MultiIndex.from_tuples(summary.columns)
    for var in summary.columns.get_level_values(0).unique():
        summary[(var, "difference")] = summary[(var, names[0])] - summary[(var, names[1])]
    summary = summary.sort_index(axis=1, level=0, sort_remaining=False)
    summary.to_csv(os.path.join(out, "summary.csv"))
    if ("T_10m", names[0]) not in summary:
        summary[("T_10m", names[0])] = np.nan
        summary[("T_10m", names[1])] = np.nan
    plot_summary(summary, names, out)

    # markdown table of the differences
    labels = {**ANNUAL_VARS, "T_10m": "T 10 m (°C)"}
    cols = [v for v in list(ANNUAL_VARS) + ["T_10m"] if (v, "difference") in summary]
    lines = [f"Mean annual difference, {names[0]} − {names[1]} (m w.e. yr⁻¹ unless noted)", "",
             "| Station | " + " | ".join(labels[v] for v in cols) + " |",
             "|---" * (len(cols) + 1) + "|"]
    for st in summary.index:
        lines.append(f"| {st} | " + " | ".join(
            f"{summary.loc[st, (v, 'difference')]:+.3f}" for v in cols) + " |")
    with open(os.path.join(out, "summary.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("Results written to", os.path.abspath(out))


if __name__ == "__main__":
    main()
