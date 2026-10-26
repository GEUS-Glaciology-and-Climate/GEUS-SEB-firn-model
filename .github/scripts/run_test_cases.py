"""Run the GEUS SEB-firn model on the test forcing files.

Used by the PR comparison workflow to run the same test cases on two
checkouts of the model (e.g. main and the PR branch). Only the
run_SEB_firn() arguments shared by all model versions are used.

Usage:
    python run_test_cases.py --model-dir <checkout> --data-dir <test_data> --out <dir>
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
    that main_SEB_firn.py assigns Tdeep = TDEEP_K at every station."""
    ds = xr.Dataset(
        {"T10m": (("time", "latitude", "longitude"),
                  np.full((2, 2, 2), TDEEP_K - 273.15, dtype="float32"))},
        coords={"time": pd.to_datetime(["1980-01-01", "1990-12-01"]),
                "latitude": [58.0, 84.0],
                "longitude": [-75.0, -10.0]},
    )
    ds.to_netcdf(path)


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

    stations = sorted(os.path.basename(f)[:-3]
                      for f in glob.glob(os.path.join(data_dir, "*.nc")))
    failed = []
    for station in stations:
        print(f"===== running {station} =====", flush=True)
        run_SEB_firn(
            station=station,
            surface_input_path=os.path.join(data_dir, f"{station}.nc"),
            output_path=out_dir,
            spin_up_path=os.path.join(out_dir, "no_spin_up"),
            plot=False,
            silent=True,
        )
        # run_SEB_firn catches model errors and only writes error.txt
        if not glob.glob(os.path.join(out_dir, f"{station}_*", f"{station}_surface.nc")):
            failed.append(station)

    if failed:
        print("Runs failed for:", failed)
        sys.exit(1)


if __name__ == "__main__":
    main()
