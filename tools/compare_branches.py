"""Compare the model output of the current working tree with the main branch.

Run it locally (e.g. by double-clicking tools/compare_branches.bat). It:
  1. exports the base branch (default origin/main) and the current working
     tree, including uncommitted changes, to a temporary folder,
  2. runs both model versions on the test cases in test_data/
     (tools/test_runs/run_test_cases.py),
  3. writes figures, a summary table and a PR description template to
     docs/PR_results/<date>_<name>/ (tools/test_runs/compare_runs.py).
Commit that folder with the PR: the figures are then shown in the PR
description when PR_description.md is pasted into it.

Usage:
    python tools/compare_branches.py --name <short_pr_name>
        [--base origin/main] [--out <dir>] [--keep-temp]
"""
import argparse
import datetime
import io
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SCRIPT_DIR = os.path.join(REPO, "tools", "test_runs")


def git(*args):
    return subprocess.run(["git", *args], cwd=REPO, check=True,
                          capture_output=True, text=True).stdout.strip()


def export_ref(ref, dest):
    """Write the files of a git ref (branch, commit) to dest."""
    archive = subprocess.run(["git", "archive", "--format=tar", ref], cwd=REPO,
                             check=True, capture_output=True).stdout
    with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
        tar.extractall(dest)


def export_working_tree(dest):
    """Copy the tracked and untracked (not ignored) files of the working tree,
    including uncommitted changes, to dest."""
    files = git("ls-files", "--cached", "--others", "--exclude-standard", "-z").split("\0")
    for f in files:
        src = os.path.join(REPO, f)
        if f and os.path.isfile(src):  # skips files deleted in the working tree
            os.makedirs(os.path.dirname(os.path.join(dest, f)), exist_ok=True)
            shutil.copy2(src, os.path.join(dest, f))


def github_url_of(folder):
    """URL of a repository folder on the current branch on GitHub, or None."""
    try:
        remote = git("remote", "get-url", "origin")
        branch = git("branch", "--show-current")
    except subprocess.CalledProcessError:
        return None
    match = re.search(r"github\.com[:/](.+?)(\.git)?$", remote)
    if not match or not branch:
        return None
    rel = os.path.relpath(folder, REPO).replace(os.sep, "/")
    return f"https://github.com/{match.group(1)}/blob/{branch}/{rel}"


def run_version(label, model_dir, data_dir, runs_dir, log_dir):
    print(f"--- running test cases with {label} (log: {log_dir}/{label}.log)", flush=True)
    with open(os.path.join(log_dir, f"{label}.log"), "w") as log:
        result = subprocess.run(
            [sys.executable, os.path.join(SCRIPT_DIR, "run_test_cases.py"),
             "--model-dir", model_dir, "--data-dir", data_dir, "--out", runs_dir],
            stdout=log, stderr=subprocess.STDOUT,
            env={**os.environ, "PYTHONWARNINGS": "ignore", "MPLBACKEND": "Agg"})
    print(f"    {label}: {'ok' if result.returncode == 0 else 'FAILED'}", flush=True)
    return result.returncode == 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", help="short name of the PR, e.g. melting-fix")
    parser.add_argument("--base", default="origin/main", help="git ref to compare with")
    parser.add_argument("--out", help="output folder (default docs/PR_results/<date>_<name>)")
    parser.add_argument("--no-fetch", action="store_true", help="do not run git fetch")
    parser.add_argument("--keep-temp", action="store_true", help="keep the model runs")
    args = parser.parse_args()

    name = args.name or input("Short name of the PR (e.g. melting-fix): ")
    name = re.sub(r"[^A-Za-z0-9_.-]+", "-", name.strip()).strip("-") or "unnamed"
    out = os.path.abspath(args.out or os.path.join(
        REPO, "docs", "PR_results", f"{datetime.date.today():%Y-%m-%d}_{name}"))
    head_name = git("branch", "--show-current") or "working tree"

    if not args.no_fetch:
        try:
            git("fetch", "origin")
        except subprocess.CalledProcessError as e:
            print("Warning: git fetch failed, using the local copy of", args.base, "\n", e.stderr)

    tmp = tempfile.mkdtemp(prefix="seb_firn_compare_")
    try:
        base_dir, head_dir = os.path.join(tmp, "base"), os.path.join(tmp, "head")
        print(f"Exporting {args.base} and the working tree of {head_name} to {tmp}")
        export_ref(args.base, base_dir)
        export_working_tree(head_dir)

        # the test cases are taken from the working tree, so that they also
        # exist when the base branch does not have them yet
        data_dir = os.path.join(head_dir, "test_data")
        base_ok = run_version("base", base_dir, data_dir, os.path.join(tmp, "runs_base"), tmp)
        head_ok = run_version("head", head_dir, data_dir, os.path.join(tmp, "runs_head"), tmp)

        os.makedirs(out, exist_ok=True)
        for label, ok in [("base", base_ok), ("head", head_ok)]:
            if not ok:  # keep the log of failed runs next to the results
                shutil.copy(os.path.join(tmp, f"{label}.log"),
                            os.path.join(out, f"run_{label}.log"))
        cmd = [sys.executable, os.path.join(SCRIPT_DIR, "compare_runs.py"),
               "--base-runs", os.path.join(tmp, "runs_base"),
               "--head-runs", os.path.join(tmp, "runs_head"),
               "--out", out,
               "--base-name", args.base.replace("origin/", ""),
               "--head-name", head_name]
        url = github_url_of(out)
        if url:
            cmd += ["--image-url-prefix", url]
        subprocess.run(cmd, check=True)
    finally:
        if args.keep_temp:
            print("Model runs kept in", tmp)
        else:
            shutil.rmtree(tmp, ignore_errors=True)

    print(f"\nResults written to {out}")
    print("Commit that folder with the PR and paste PR_description.md into the PR description.")
    if not base_ok:
        print("Warning: the test cases failed with", args.base)
    if not head_ok:
        print("ERROR: the test cases failed with", head_name)
        sys.exit(1)


if __name__ == "__main__":
    main()
