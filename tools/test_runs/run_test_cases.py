"""Run one version of the GEUS SEB-firn model on all test cases.

Every file test_data/<site>_<year>.nc is one test case: a one-year forcing
for one site, in the same format as CARRA_at_AWS.nc. Each case is run
without spin-up (initial profiles from input/initial state) and with a
constant deep firn temperature, so that two model versions are run in
exactly the same conditions.

Only the run_SEB_firn() arguments shared by all model versions are used, so
this script can run older versions (e.g. main) as well as new ones.

Usage:
    python run_test_cases.py --model-dir <model folder> --data-dir <test_data> --out <dir>

Output: <out>/<site>_<year>/<site>_100_layers_3h/<site>_surface.nc and <site>_slwc.nc
Exit code 1 if any test case failed.
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr

# constant deep firn temperature used for all test runs (K)
TDEEP_K = 263.0


def write_constant_T10m(path):
    """Write a T10m_prediction.nc covering Greenland with a uniform value so
    that main_SEB_firn.py assigns Tdeep = TDEEP_K at every site."""
    ds = xr.Dataset(
        {"T10m": (("time", "latitude", "longitude"),
                  np.full((2, 2, 2), TDEEP_K - 273.15, dtype="float32"))},
        coords={"time": pd.to_datetime(["1980-01-01", "1990-12-01"]),
                "latitude": [58.0, 84.0],
                "longitude": [-75.0, -10.0]},
    )
    ds.to_netcdf(path)


def list_test_cases(data_dir):
    """Return [(case name, site, path)] for all files in data_dir."""
    cases = []
    for path in sorted(glob.glob(os.path.join(data_dir, "*.nc"))):
        name = os.path.basename(path)[:-3]
        with xr.open_dataset(path) as ds:
            site = str(np.atleast_1d(ds.stid.values)[0])
        cases.append((name, site, path))
    return cases


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    data_dir = os.path.abspath(args.data_dir)
    out_dir = os.path.abspath(args.out)
    os.makedirs(out_dir, exist_ok=True)

    # the model uses paths relative to its root folder
    os.chdir(os.path.abspath(args.model_dir))
    sys.path.insert(0, os.getcwd())
    write_constant_T10m("input/T10m_prediction.nc")

    from main_SEB_firn import run_SEB_firn

    failed = []
    for name, site, path in list_test_cases(data_dir):
        print(f"===== running {name} =====", flush=True)
        case_dir = os.path.join(out_dir, name)
        run_SEB_firn(
            station=site,
            surface_input_path=path,
            output_path=case_dir,
            spin_up_path=os.path.join(case_dir, "no_spin_up"),
            output_var=["slwc"],  # surface output + liquid water profiles
            plot=False,
            silent=True,
        )
        # run_SEB_firn catches model errors and only writes error.txt
        if not glob.glob(os.path.join(case_dir, "*", f"{site}_surface.nc")):
            failed.append(name)

    if failed:
        print("Test cases failed:", failed)
        sys.exit(1)


if __name__ == "__main__":
    main()
