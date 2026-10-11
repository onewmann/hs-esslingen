"""Run time of estimate_velocity: first call and repeated calls, by image size.

    python experiments/timing.py [--quick] [--out DIR] [--replot] [--no-octave] [--matlab]
    python experiments/timing.py --sizes 512 2048 --repeats 5 --out DIR

Times ``estimate_velocity(bp)`` with the defaults (L = 4, pad_factor 2) on
simulated n x n pulsed images, n = 128, 256, 512 and 1024 (--quick: up to
512, fewer runs). Every first call runs in a fresh process: it includes the
bin weights of the projection, which are computed once per image size and
then cached, and the set-up of the FFT, as a single call of the command line
or a script does. One of these processes then repeats the call, as a live
loop over camera frames does. The table gives the median of each. The
Python process also reports the time to import fftvel with NumPy.

The same images are timed in Octave (octave-cli or octave on the PATH; skip
with --no-octave) and, with --matlab, in MATLAB (matlab -batch). The calls
run one after the other, never in parallel, so that they do not compete for
the CPU; close other programs for a clean measurement.

Writes timing.csv and environment.json (platform, CPU count, versions) to
the output folder: experiments/results for a full run, the folder that the
README quotes; experiments/results_quick with --quick; DIR with --out.
--replot prints the table of timing.csv again. There is no figure.
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import platform
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# common (and through it fftvel) is imported inside the functions, so that the
# worker process can time the import of fftvel itself.
HERE = Path(__file__).resolve().parent
REPO = HERE.parent
SRC = REPO / "python" / "src"

CSV = ["timing.csv"]
SIZES = (128, 256, 512, 1024)
QUICK_SIZES = (128, 256, 512)
FIRST_RUNS = {128: 5, 256: 5, 512: 3}          # fresh processes per size (3 for larger sizes)
REPEATS = {128: 20, 256: 10, 512: 5}           # repeated calls per size (3 for larger sizes)
QUICK_FIRST_RUNS, QUICK_REPEATS = 2, 3
VELOCITY = (3.0, -2.0)                         # px per step; the run time hardly depends on it
TIMEOUT = 3600                                 # s per process

# M-code run by Octave and MATLAB: load the image, time the first call and the repeats.
M_CODE = ("addpath('{matlab}'); d = load('{mat}'); bp = d.bp; "
          "t0 = tic; fftvel.estimate_velocity(bp); first = toc(t0); "
          "reps = zeros(1, {repeats}); "
          "for k = 1:{repeats}, t0 = tic; fftvel.estimate_velocity(bp); reps(k) = toc(t0); end; "
          "fprintf('VERSION %s\\n', version); "
          "fprintf('TIMING %.9g', first); fprintf(' %.9g', reps); fprintf('\\n');")


def worker(path: str, repeats: int) -> None:
    """Fresh Python process: time the import, the first call and ``repeats`` further calls."""
    t0 = time.perf_counter()
    sys.path.insert(0, str(SRC))
    import numpy as np

    import fftvel
    t_import = time.perf_counter() - t0
    bp = np.load(path)
    t0 = time.perf_counter()
    fftvel.estimate_velocity(bp)
    first = time.perf_counter() - t0
    reps = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fftvel.estimate_velocity(bp)
        reps.append(time.perf_counter() - t0)
    print(json.dumps({"import_s": t_import, "first_s": first, "repeated_s": reps}))


def _python_run(npy: Path, repeats: int) -> dict:
    r = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker", str(npy), str(repeats)],
                       capture_output=True, text=True, timeout=TIMEOUT)
    if r.returncode != 0:
        raise RuntimeError(f"Python worker failed:\n{r.stderr}")
    return json.loads(r.stdout.strip().splitlines()[-1])


def m_command(tool: str):
    """Command line that runs M-code in Octave or MATLAB, or None if the program is not on the PATH."""
    if tool == "octave":
        exe = shutil.which("octave-cli") or shutil.which("octave")
        return None if exe is None else [exe, "--no-gui", "--quiet", "--no-window-system", "--eval"]
    exe = shutil.which("matlab")
    return None if exe is None else [exe, "-batch"]


def _m_run(cmd, mat: Path, repeats: int) -> dict:
    code = M_CODE.format(matlab=(REPO / "matlab").as_posix().replace("'", "''"),
                         mat=mat.as_posix().replace("'", "''"), repeats=repeats)
    r = subprocess.run([*cmd, code], capture_output=True, text=True, timeout=TIMEOUT)
    lines = r.stdout.splitlines()
    timing = [ln.split()[1:] for ln in lines if ln.startswith("TIMING ")]
    if r.returncode != 0 or not timing:
        raise RuntimeError(f"{cmd[0]} failed (exit {r.returncode}):\n{r.stdout}\n{r.stderr}")
    version = next((ln[8:].strip() for ln in lines if ln.startswith("VERSION ")), "")
    values = [float(x) for x in timing[0]]
    return {"first_s": values[0], "repeated_s": values[1:], "version": version}


def measure(implementation: str, n: int, files: dict, first_runs: int, repeats: int, cmd=None) -> dict:
    """CSV row: medians over ``first_runs`` fresh processes and ``repeats`` repeated calls."""
    from fftvel.core import DEFAULT_L, DEFAULT_PAD_FACTOR

    runs = []
    for k in range(first_runs):
        reps = repeats if k == 0 else 0
        if implementation == "python":
            runs.append(_python_run(files["npy"], reps))
        else:
            runs.append(_m_run(cmd, files["mat"], reps))
    first = [r["first_s"] for r in runs]
    rep = runs[0]["repeated_s"]
    row = dict(implementation=implementation, n=n, npad=DEFAULT_PAD_FACTOR * n, L=DEFAULT_L,
               first_runs=len(first),
               first_median_s=statistics.median(first), first_min_s=min(first),
               repeated_runs=len(rep),
               repeated_median_s=statistics.median(rep) if rep else float("nan"),
               repeated_min_s=min(rep) if rep else float("nan"))
    if implementation == "python":
        row["import_median_s"] = statistics.median(r["import_s"] for r in runs)
    else:
        row["version"] = runs[0]["version"]
    return row


def environment(versions: dict) -> dict:
    import numpy as np

    import fftvel
    env = dict(date_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
               platform=platform.platform(), machine=platform.machine(), processor=platform.processor(),
               cpu_count=os.cpu_count(), python=platform.python_version(),
               python_implementation=platform.python_implementation(), numpy=np.__version__,
               fftvel=fftvel.__version__)
    try:
        import scipy
        env["scipy"] = scipy.__version__
    except ImportError:
        pass
    env.update(versions)
    return env


def _log(msg: str) -> None:
    print(msg, flush=True)


def run(quick: bool, sizes=None, first_runs=None, repeats=None, tools=("octave",), log=_log):
    """Time all sizes in Python and in the M-code tools found; write timing.csv and environment.json."""
    import numpy as np
    from scipy.io import savemat

    import common
    from fftvel import simulate_pulsed_image
    from fftvel.core import DEFAULT_L

    sizes = sizes or (QUICK_SIZES if quick else SIZES)
    cmds = {}
    for tool in tools:
        cmd = m_command(tool)
        if cmd is None:
            log(f"{tool}: not found on the PATH, not timed")
        else:
            cmds[tool] = cmd
    rows = []
    versions = {}
    with tempfile.TemporaryDirectory() as tmp:
        for n in sizes:
            F = first_runs or (QUICK_FIRST_RUNS if quick else FIRST_RUNS.get(n, 3))
            R = repeats or (QUICK_REPEATS if quick else REPEATS.get(n, 3))
            bp = simulate_pulsed_image(n, VELOCITY, DEFAULT_L, rng=n)
            files = {"npy": Path(tmp) / f"bp{n}.npy", "mat": Path(tmp) / f"bp{n}.mat"}
            np.save(files["npy"], bp)
            savemat(files["mat"], {"bp": bp})
            for impl in ("python", *cmds):
                log(f"{impl:7s} n = {n:4d}: {F} first calls, {R} repeated calls")
                row = measure(impl, n, files, F, R, cmds.get(impl))
                if impl != "python":
                    versions[f"{impl}_version"] = row["version"]
                rows.append(row)
    common.write_csv("timing.csv", rows)
    common.RESULTS.mkdir(parents=True, exist_ok=True)
    with open(common.RESULTS / "environment.json", "w") as f:
        json.dump(environment(versions), f, indent=2)
        f.write("\n")
    summary(rows)
    return rows


def summary(rows):
    """Print the table."""
    print("\nRun time of estimate_velocity (L = 4, pad_factor 2), median in s:")
    print(f"  {'':7s} {'n':>5s} {'npad':>5s}  {'first call':>14s}  {'repeated':>14s}  {'import':>6s}")
    for r in rows:
        imp = r.get("import_median_s", "")
        imp = f"{imp:6.2f}" if isinstance(imp, float) else f"{'':6s}"
        print(f"  {r['implementation']:7s} {int(r['n']):5d} {int(r['npad']):5d}  "
              f"{r['first_median_s']:8.3f} ({int(r['first_runs']):2d})   "
              f"{r['repeated_median_s']:8.3f} ({int(r['repeated_runs']):2d})   {imp}")
    print("  (runs in brackets; first call: fresh process, import of fftvel not included)")


def replot():
    import common
    summary(common.read_csv("timing.csv"))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--quick", action="store_true",
                   help="sizes up to 512, fewer runs; writes to experiments/results_quick unless --out is given")
    p.add_argument("--out", type=Path, metavar="DIR",
                   help="folder for timing.csv and environment.json (default: experiments/results; with "
                        "--quick experiments/results_quick)")
    p.add_argument("--replot", action="store_true", help="print the table of timing.csv in the output folder")
    p.add_argument("--sizes", type=int, nargs="+", metavar="N", help="image sizes in px (default: 128 256 512 1024)")
    p.add_argument("--first-runs", type=int, metavar="F", help="fresh processes per size (first calls)")
    p.add_argument("--repeats", type=int, metavar="R", help="repeated calls per size")
    p.add_argument("--no-octave", action="store_true", help="do not time Octave")
    p.add_argument("--matlab", action="store_true", help="also time MATLAB (matlab -batch, R2019a or newer)")
    a = p.parse_args(argv)
    if a.sizes and min(a.sizes) < 16:
        p.error("--sizes must be at least 16")
    for name in ("first_runs", "repeats"):
        if getattr(a, name) is not None and getattr(a, name) < 1:
            p.error(f"--{name.replace('_', '-')} must be at least 1")

    sys.path.insert(0, str(HERE))
    import common
    out = (a.out if a.out is not None else common.QUICK if a.quick else common.COMMITTED).resolve()
    if a.quick and out == common.COMMITTED:
        p.error("--quick must not write to experiments/results, which holds the full-run results; "
                "leave out --out or choose another folder")
    common.set_output(out)
    print(f"CSV files: {common.RESULTS}")
    if a.replot:
        missing = common.missing_csv(CSV)
        if missing:
            sys.exit(f"--replot: {', '.join(missing)} not found in {common.RESULTS}")
        replot()
        return
    tools = () if a.no_octave else ("octave",)
    if a.matlab:
        tools += ("matlab",)
    run(a.quick, a.sizes, a.first_runs, a.repeats, tools)


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "--worker":
        worker(sys.argv[2], int(sys.argv[3]))
    else:
        main()
