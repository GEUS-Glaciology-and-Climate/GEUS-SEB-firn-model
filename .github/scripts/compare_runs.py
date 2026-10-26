"""Compare melt and SMB between two sets of test runs.

Reads the <station>_surface.nc files written by run_test_cases.py for a
reference run (e.g. main) and a new run (e.g. the PR branch) and writes a
markdown report with annual and total sums per station.

Usage:
    python compare_runs.py --ref <dir> --new <dir> --report <file.md>
        [--ref-name main] [--new-name PR]
"""
import argparse
import glob
import os

import xarray as xr

VARIABLES = {"melt_mweq": "Melt", "smb_mweq": "SMB"}


def load_annual_sums(run_dir):
    """Return {station: DataFrame of annual sums (m w.e.)} for one run."""
    sums = {}
    for f in sorted(glob.glob(os.path.join(run_dir, "*", "*_surface.nc"))):
        with xr.open_dataset(f) as ds:
            station = ds.attrs["station"]
            df = ds[list(VARIABLES)].to_dataframe()
        sums[station] = df.resample("YS").sum()
        sums[station].index = sums[station].index.year
    return sums


def fmt(x):
    return f"{x:.3f}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", required=True)
    parser.add_argument("--new", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--ref-name", default="main")
    parser.add_argument("--new-name", default="PR")
    args = parser.parse_args()

    ref = load_annual_sums(args.ref)
    new = load_annual_sums(args.new)
    stations = sorted(set(ref) | set(new))

    lines = [f"## Test runs: {args.new_name} vs {args.ref_name}", ""]
    if not ref:
        lines += [f"No output found for **{args.ref_name}** (run failed?).", ""]
    if not new:
        lines += [f"No output found for **{args.new_name}** (run failed?).", ""]

    # summary over the whole period
    lines += ["### Total over the test period (m w.e.)", "",
              "| Station | Variable | " + args.ref_name + " | " + args.new_name
              + " | Difference | Difference (%) |",
              "|---|---|---|---|---|---|"]
    for station in stations:
        for var, label in VARIABLES.items():
            r = ref[station][var].sum() if station in ref else float("nan")
            n = new[station][var].sum() if station in new else float("nan")
            pct = 100 * (n - r) / abs(r) if r not in (0,) else float("nan")
            lines.append(f"| {station} | {label} | {fmt(r)} | {fmt(n)} | "
                         f"{fmt(n - r)} | {pct:.1f} |")
    lines.append("")

    # annual detail per station
    for station in stations:
        if station not in ref or station not in new:
            continue
        lines += [f"<details><summary>{station}: annual sums (m w.e.)</summary>",
                  "",
                  "| Year | Melt " + args.ref_name + " | Melt " + args.new_name
                  + " | SMB " + args.ref_name + " | SMB " + args.new_name + " |",
                  "|---|---|---|---|---|"]
        for year in ref[station].index.union(new[station].index):
            r = ref[station].reindex([year]).iloc[0]
            n = new[station].reindex([year]).iloc[0]
            lines.append(f"| {year} | {fmt(r.melt_mweq)} | {fmt(n.melt_mweq)} | "
                         f"{fmt(r.smb_mweq)} | {fmt(n.smb_mweq)} |")
        lines += ["", "</details>", ""]

    report = "\n".join(lines)
    with open(args.report, "w") as f:
        f.write(report + "\n")
    print(report)


if __name__ == "__main__":
    main()
