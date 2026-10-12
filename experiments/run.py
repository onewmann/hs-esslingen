"""Accuracy, noise robustness, validity check, test sequences, camera
recordings and before/after comparison; on request also the threshold
calibration, the baselines, the model violations and the run time.

    python experiments/run.py                 # default parts, about 15 min on 4 cores
    python experiments/run.py --all           # all ten parts, about 70 min on 4 cores
    python experiments/run.py --quick         # smaller grids, for a check
    python experiments/run.py --quick --all
    python experiments/run.py --only calibration violations
    python experiments/run.py --only recordings --recordings PATH
    python experiments/run.py --replot --only noise
    python experiments/noise.py --quick       # one part alone

Parts run without --only (the default set): accuracy, noise, before_after,
pipeline, testsequences and recordings. --all adds the four parts that run
only on request (or name them with --only):

* calibration: the validity threshold (calibrate_threshold.py, about 20 min
  on 4 cores),
* baselines: simple single-image baselines on the images of the noise study
  (baselines.py, 7 to 10 min on 4 cores),
* violations: an independent forward model of the camera and the factors
  that break the estimator (model_violations.py, about 20 min on 4 cores),
* timing: run time of estimate_velocity in Python and, if octave is on the
  PATH, in Octave (timing.py, about 8 min; it runs last and alone, so that
  nothing else competes for the CPU).

A full run writes CSV files to experiments/results/ and figures to
docs/figures/, the files the README quotes and shows. --quick writes both to
experiments/results_quick/ (figures in its figures/ folder), which git
ignores, so a check never overwrites the committed results. --out DIR writes
the CSV files to DIR and the figures to DIR/figures instead. --replot redraws
the figures of the selected parts from the CSV files in that folder (the
pipeline figure has no CSV file and is computed again; timing has no figure
and prints its table).

``--recordings`` is the folder with the camera recordings vid1.mat to
vid4.mat (default: data/recordings, see data/README.md). The recordings part
analyses the frame-to-frame motion of vid2 to vid4 and has no smaller
variant for --quick; the calibration part takes single frames of all four.

Each part is a module that also runs alone: accuracy.py, noise.py,
history/before_after.py, pipeline_figure.py, testsequences.py,
recordings_2025.py (part recordings), calibrate_threshold.py (part
calibration), baselines.py, model_violations.py (part violations) and
timing.py. ``python experiments/noise.py ARGS`` is ``python
experiments/run.py --only noise ARGS``; the last four have options of their
own (timing.py: --matlab, --no-octave, --sizes). common.py holds the paths,
the plot style, the simulated case and the CSV files. Your own
measurements (multi-flash images, or frame stacks in average mode) are
evaluated by real_measurements.py, which is not a part of run.py (see
data/measurements/README.md).

Needs scipy and matplotlib 3.6 or newer (python -m pip install -e "./python[plot,data]").
"""

from __future__ import annotations

import argparse
import os
import sys
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import accuracy  # noqa: E402
import baselines  # noqa: E402
import calibrate_threshold  # noqa: E402
import common  # noqa: E402
import model_violations  # noqa: E402
import noise  # noqa: E402
import pipeline_figure  # noqa: E402
import recordings_2025  # noqa: E402
import testsequences  # noqa: E402
import timing  # noqa: E402
from history import before_after  # noqa: E402

PARTS = ["accuracy", "noise", "before_after", "pipeline", "testsequences", "recordings"]   # the default set
EXTRA = ["calibration", "baselines", "violations", "timing"]           # only with --only or --all
ALL_PARTS = PARTS + EXTRA
STUDIES = dict(zip(ALL_PARTS, (accuracy, noise, before_after, pipeline_figure, testsequences, recordings_2025,
                               calibrate_threshold, baselines, model_violations, timing), strict=True))
SIMULATIONS = ["accuracy", "noise", "before_after", "calibration", "baselines", "violations"]   # share the pool
SUMMARISED = ["accuracy", "noise", "before_after"]      # summary also for the CSV files of an earlier run
RECORDING_FILES = {"recordings": recordings_2025.FILES, "calibration": calibrate_threshold.FILES}


def replot(parts, requested):
    """Redraw the figures of the given parts from the CSV files in the output folder."""
    for part in ALL_PARTS:
        if part not in parts:
            continue
        missing = common.missing_csv(STUDIES[part].CSV)
        if missing:
            if part in requested:
                sys.exit(f"--replot {part}: {', '.join(missing)} not found in {common.RESULTS}")
            print(f"{part}: skipped, {', '.join(missing)} not found in {common.RESULTS}")
            continue
        STUDIES[part].replot()


def main(argv=None, only=None, description=__doc__):
    """Command line of the experiments; ``only`` fixes the parts when a study module runs alone."""
    offered = only or ALL_PARTS
    p = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--quick", action="store_true",
                   help="smaller grids; writes to experiments/results_quick unless --out is given")
    p.add_argument("--out", type=Path, metavar="DIR",
                   help="folder for the CSV files, figures go to DIR/figures (default: experiments/results "
                        "and docs/figures; with --quick experiments/results_quick)")
    needed = sorted({f for part in offered if part in RECORDING_FILES for f in RECORDING_FILES[part]})
    if needed:
        p.add_argument("--recordings", type=Path, metavar="PATH",
                       help=f"folder with {', '.join(needed)} (default: data/recordings)")
    if set(offered) & set(SIMULATIONS):
        p.add_argument("--jobs", type=int, default=os.cpu_count() or 1,
                       help="worker processes for the simulations (default: number of CPUs)")
    p.add_argument("--replot", action="store_true",
                   help="redraw the figures from the CSV files in the output folder")
    if only is None:
        which = p.add_mutually_exclusive_group()
        which.add_argument("--only", nargs="+", choices=ALL_PARTS, metavar="PART",
                           help="run (or with --replot redraw) only these parts: " + ", ".join(ALL_PARTS)
                                + " (default: " + ", ".join(PARTS) + ")")
        which.add_argument("--all", action="store_true",
                           help="run (or with --replot redraw) all parts, also " + ", ".join(EXTRA))
    a = p.parse_args(argv)
    chosen = only or a.only
    requested = set(chosen or [])
    parts = set(ALL_PARTS if getattr(a, "all", False) else chosen or PARTS)
    if getattr(a, "jobs", 1) < 1:
        p.error("--jobs must be at least 1")
    given = getattr(a, "recordings", None)
    folder = given if given is not None else common.RECORDINGS
    if not a.replot:
        uses = [part for part in RECORDING_FILES if part in parts]
        if given is not None and not uses:
            p.error("--recordings is given, but --only leaves out " + " and ".join(RECORDING_FILES))
        missing = [f for f in sorted({f for part in uses for f in RECORDING_FILES[part]})
                   if not (folder / f).is_file()]
        if missing:
            p.error(f"--recordings {folder}: {', '.join(missing)} not found")
    out = (a.out if a.out is not None else common.QUICK if a.quick else common.COMMITTED).resolve()
    if a.quick and out == common.COMMITTED:
        p.error("--quick must not write to experiments/results, which holds the full-run results; "
                "leave out --out or choose another folder")
    if a.quick and out / "figures" == common.README_FIGURES:
        p.error("--quick must not write to docs/figures, which holds the figures of the full run; "
                "leave out --out or choose another folder")
    common.set_output(out)
    import matplotlib
    matplotlib.use("Agg")
    common.FIGURES.mkdir(parents=True, exist_ok=True)
    print(f"CSV files: {common.RESULTS}\nfigures:   {common.FIGURES}")
    if a.replot:
        replot(parts, requested)
        return
    sim = parts & set(SIMULATIONS)
    if sim:
        with Pool(a.jobs, initializer=common.init_worker) as pool:
            if sim & set(SUMMARISED):
                acc = accuracy.run(pool, a.quick) if "accuracy" in sim else common.read_if("accuracy.csv")
                if "noise" in sim:
                    noi, stat = noise.run(pool, a.quick)
                else:
                    noi, stat = common.read_if("noise.csv"), common.read_if("static.csv")
                ba = before_after.run(pool, a.quick) if "before_after" in sim else common.read_if("before_after.csv")
                # the numbers the README quotes, also for parts of an earlier run in the same folder
                if acc is not None:
                    accuracy.summary(acc)
                noise.summary(noi, stat)
                if ba is not None:
                    before_after.summary(ba)
            if "calibration" in sim:
                calibrate_threshold.run(pool, a.quick, folder.resolve())     # prints its summary
            if "baselines" in sim:
                baselines.summary(baselines.run(pool, a.quick))
            if "violations" in sim:
                model_violations.summary(model_violations.run(pool, a.quick), check=False)   # run printed the check
    if "pipeline" in parts:
        pipeline_figure.run()
    if "testsequences" in parts:
        testsequences.summary(testsequences.run())
    if "recordings" in parts:
        recordings_2025.run(folder)
    if "timing" in parts:
        timing.run(a.quick)         # last, after the worker pool has closed


if __name__ == "__main__":
    main()
