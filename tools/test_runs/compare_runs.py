"""Compare two sets of test runs made by run_test_cases.py.

For every test case (<site>_<year>) found in either run, writes to --out:
    <site>_<year>_cumulative_melt_smb.png
    <site>_<year>_cumulative_runoff_refreezing.png
    <site>_<year>_lwc_profiles.png   total LWC, LWC profiles of both runs and
                                     their difference
    summary.md            table of annual totals and differences
    PR_description.md     template for the PR description, with the table
                          and the figures (if --image-url-prefix is given)

Usage:
    python compare_runs.py --base-runs <dir> --head-runs <dir> --out <dir>
        [--base-name main] [--head-name PR] [--image-url-prefix <url> ...]
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# variables compared, all in m w.e. per time step in <site>_surface.nc
VARIABLES = {
    "melt_mweq": "Melt",
    "runoff_mweq": "Runoff",
    "refreezing_mweq": "Refreezing",
    "sublimation_mweq": "Sublimation",
    "smb_mweq": "SMB",
}
FIGURES = {
    "cumulative_melt_smb": ["melt_mweq", "smb_mweq"],
    "cumulative_runoff_refreezing": ["runoff_mweq", "refreezing_mweq"],
}


def list_cases(*run_dirs):
    cases = set()
    for run_dir in run_dirs:
        for f in glob.glob(os.path.join(run_dir, "*", "*", "*_surface.nc")):
            cases.add(os.path.basename(os.path.dirname(os.path.dirname(f))))
    return sorted(cases)


def load_case(run_dir, case):
    """Return the surface output of one test case as a DataFrame, or None."""
    files = glob.glob(os.path.join(run_dir, case, "*", "*_surface.nc"))
    if not files:
        return None
    with xr.open_dataset(files[0]) as ds:
        variables = [v for v in VARIABLES if v in ds]
        return ds[variables].to_dataframe()


def plot_case(case, base, head, base_name, head_name, out_dir):
    """Write the figures of one test case and return their file names."""
    names = []
    for fig_name, variables in FIGURES.items():
        fig, axes = plt.subplots(len(variables), 1, figsize=(8, 6), sharex=True)
        for ax, var in zip(np.atleast_1d(axes), variables):
            for df, label, style in [(base, base_name, dict(color="0.4", lw=2.5)),
                                     (head, head_name, dict(color="tab:red", lw=1.2))]:
                if df is not None and var in df:
                    cum = df[var].cumsum()
                    ax.plot(cum.index, cum.values,
                            label=f"{label}: {cum.iloc[-1]:.3f}", **style)
            ax.set_ylabel(f"Cumulative {VARIABLES[var]}\n(m w.e.)")
            ax.legend(loc="best", fontsize=9)
            ax.grid(alpha=0.3)
        np.atleast_1d(axes)[0].set_title(case)
        fig.autofmt_xdate()
        fig.tight_layout()
        name = f"{case}_{fig_name}.png"
        fig.savefig(os.path.join(out_dir, name), dpi=100)
        plt.close(fig)
        names.append(name)
    return names


def load_lwc(run_dir, case):
    """Return the liquid water profiles (<site>_slwc.nc) of one test case, or None."""
    files = glob.glob(os.path.join(run_dir, case, "*", "*_slwc.nc"))
    if not files:
        return None
    with xr.open_dataset(files[0]) as ds:
        return ds[["slwc", "depth"]].load()


def lwc_on_depth_grid(ds, grid):
    """Liquid water content (kg m-3) of the layers interpolated, as a step
    function, on a regular depth grid (m) at each time step. The model layers
    move with the snow and are merged/split, so the two runs only share a grid
    after this interpolation."""
    slwc = ds.slwc.transpose("level", "time").values  # m w.e. per layer
    bottom = ds.depth.transpose("level", "time").values  # m, depth of layer bottom
    top = np.vstack([np.zeros((1, bottom.shape[1])), bottom[:-1]])
    density = 1000 * slwc / np.where(bottom - top > 0, bottom - top, np.nan)
    out = np.full((len(grid), bottom.shape[1]), np.nan)
    for t in range(bottom.shape[1]):
        idx = np.searchsorted(bottom[:, t], grid)
        valid = idx < bottom.shape[0]
        out[valid, t] = density[idx[valid], t]
    return out


def plot_lwc(case, base, head, base_name, head_name, out_dir):
    """Four-panel figure: total LWC of both runs, LWC profiles of each run,
    and their difference. Returns the file name, or None if no profile output."""
    if base is None or head is None:
        return None
    # depth range: deepest layer holding water in either run, at least 1 m
    wet = [ds.depth.where(ds.slwc > 1e-6).max().item() for ds in (base, head)]
    max_depth = min(np.nanmax(wet + [1.0]) * 1.2, 20.0)
    grid = np.arange(0, max_depth, max_depth / 200) + max_depth / 400
    lwc_base = lwc_on_depth_grid(base, grid)
    lwc_head = lwc_on_depth_grid(head, grid)
    delta = lwc_head - lwc_base

    fig, axes = plt.subplots(4, 1, figsize=(9, 11), sharex=True,
                             gridspec_kw=dict(height_ratios=[0.8, 1, 1, 1]))
    for ds, label, style in [(base, base_name, dict(color="0.4", lw=2.5)),
                             (head, head_name, dict(color="tab:red", lw=1.2))]:
        total = 1000 * ds.slwc.sum("level")  # mm w.e.
        axes[0].plot(ds.time, total, label=f"{label}: max {total.max().item():.1f} mm w.e.",
                     **style)
    axes[0].set_ylabel("Total LWC\n(mm w.e.)")
    axes[0].legend(loc="upper left", fontsize=9)
    axes[0].grid(alpha=0.3)
    axes[0].set_title(case)

    vmax = max(np.nanmax(lwc_base), np.nanmax(lwc_head), 1e-3)
    dmax = max(np.nanmax(np.abs(delta)), 1e-3)
    panels = [(axes[1], lwc_base, f"LWC {base_name}", "Blues", 0, vmax, "kg m$^{-3}$"),
              (axes[2], lwc_head, f"LWC {head_name}", "Blues", 0, vmax, "kg m$^{-3}$"),
              (axes[3], delta, f"LWC {head_name} $-$ {base_name}", "RdBu_r", -dmax, dmax,
               "kg m$^{-3}$")]
    for ax, data, title, cmap, vmin, vmx, unit in panels:
        im = ax.pcolormesh(base.time.values, grid, data, cmap=cmap, vmin=vmin, vmax=vmx,
                           shading="nearest", rasterized=True)
        ax.set_ylim(max_depth, 0)
        ax.set_ylabel("Depth (m)")
        ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, label=unit, pad=0.01)
    # keep the top panel aligned with the colorbar panels
    fig.colorbar(plt.cm.ScalarMappable(), ax=axes[0], pad=0.01).ax.set_visible(False)
    fig.autofmt_xdate()
    name = f"{case}_lwc_profiles.png"
    fig.savefig(os.path.join(out_dir, name), dpi=100, bbox_inches="tight")
    plt.close(fig)
    return name


def summary_table(cases, runs, base_name, head_name):
    lines = ["| Test case | Variable | " + base_name + " | " + head_name
             + " | Difference | Difference (%) |",
             "|---|---|---|---|---|---|"]
    identical = True
    for case in cases:
        base, head = runs[case]
        for var, label in VARIABLES.items():
            b = base[var].sum() if base is not None and var in base else np.nan
            h = head[var].sum() if head is not None and var in head else np.nan
            diff = h - b
            pct = f"{100 * diff / abs(b):+.1f}" if b else "n/a"
            if not (base is not None and head is not None and var in base and var in head
                    and np.array_equal(base[var].values, head[var].values)):
                identical = False
            lines.append(f"| {case} | {label} | {b:.3f} | {h:.3f} | {diff:+.3f} | {pct} |")
    return lines, identical


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-runs", required=True)
    parser.add_argument("--head-runs", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--base-name", default="main")
    parser.add_argument("--head-name", default="PR")
    parser.add_argument("--image-url-prefix", action="append", default=[],
                        help="URL of the --out folder on GitHub, used to embed "
                             "the figures in PR_description.md. Can be given "
                             "several times (e.g. PR branch and main): each "
                             "figure is then linked from every location, with "
                             "an empty alt text so that links that do not "
                             "resolve (yet, or any more) display nothing")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    cases = list_cases(args.base_runs, args.head_runs)
    runs = {case: (load_case(args.base_runs, case), load_case(args.head_runs, case))
            for case in cases}

    figures = {case: plot_case(case, *runs[case], args.base_name, args.head_name, args.out)
               for case in cases}
    for case in cases:
        name = plot_lwc(case, load_lwc(args.base_runs, case), load_lwc(args.head_runs, case),
                        args.base_name, args.head_name, args.out)
        if name:
            figures[case].append(name)

    table, identical = summary_table(cases, runs, args.base_name, args.head_name)
    missing = [f"- No output for **{name}** in {case} (run failed?)"
               for case in cases
               for name, df in zip([args.base_name, args.head_name], runs[case])
               if df is None]
    if not cases:
        missing = ["- No test case output found in either run."]

    header = [f"Totals over each test case in m w.e. ({args.head_name} vs {args.base_name}). "
              "Constant Tdeep = 263 K, no spin-up (initial profiles from input/initial state).",
              ""]
    verdict = (["**All compared variables are identical.**", ""] if identical and cases
               else [])
    summary = header + verdict + table + [""] + missing
    with open(os.path.join(args.out, "summary.md"), "w") as f:
        f.write("\n".join(summary) + "\n")

    description = ["## Changes", "", "<!-- describe the bug fix and why it is needed -->", "",
                   "## Effect on the test runs", ""] + summary + [""]
    if args.image_url_prefix:
        for case in cases:
            description += [f"### {case}", ""]
            for name in figures[case]:
                description += [f'<img src="{prefix}/{name}?raw=true" alt="">'
                                for prefix in args.image_url_prefix]
            description.append("")
    with open(os.path.join(args.out, "PR_description.md"), "w") as f:
        f.write("\n".join(description) + "\n")

    print("\n".join(summary))


if __name__ == "__main__":
    main()
