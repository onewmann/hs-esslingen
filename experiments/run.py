"""Accuracy, noise robustness, validity check, test sequences, camera
recordings and before/after comparison.

    python experiments/run.py                 # about 15 min on 4 cores
    python experiments/run.py --quick         # smaller grids, for a check
    python experiments/run.py --only recordings --recordings PATH
    python experiments/run.py --replot --only noise
    python experiments/noise.py --quick       # one part alone

A full run writes CSV files to experiments/results/ and figures to
docs/figures/, the files the README quotes and shows. --quick writes both to
experiments/results_quick/ (figures in its figures/ folder), which git
ignores, so a check never overwrites the committed results. --out DIR writes
the CSV files to DIR and the figures to DIR/figures instead. --replot redraws
the figures of the selected parts from the CSV files in that folder (the
pipeline figure has no CSV file and is computed again).

``--recordings`` is the folder with the camera recordings vid2.mat, vid3.mat
and vid4.mat of the 2025 student project (default: data/recordings, see
data/README.md); the recordings part analyses their frame-to-frame motion
and has no smaller variant for --quick.

Each part is a module that also runs alone: accuracy.py, noise.py,
history/before_after.py, pipeline_figure.py, test_sequences.py (part
testsequences) and recordings_2025.py (part recordings).
``python experiments/noise.py ARGS`` is ``python experiments/run.py --only
noise ARGS``. common.py holds the paths, the plot style, the simulated case
and the CSV files.

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
import common  # noqa: E402
import noise  # noqa: E402
import pipeline_figure  # noqa: E402
import recordings_2025  # noqa: E402
import test_sequences  # noqa: E402
from history import before_after  # noqa: E402

PARTS = ["accuracy", "noise", "before_after", "pipeline", "testsequences", "recordings"]
STUDIES = dict(zip(PARTS, (accuracy, noise, before_after, pipeline_figure, test_sequences,
                           recordings_2025)))
SIMULATIONS = ["accuracy", "noise", "before_after"]      # parts that share the worker pool


def replot(parts, requested):
    """Redraw the figures of the given parts from the CSV files in the output folder."""
    for part in PARTS:
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
    offered = only or PARTS
    p = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--quick", action="store_true",
                   help="smaller grids; writes to experiments/results_quick unless --out is given")
    p.add_argument("--out", type=Path, metavar="DIR",
                   help="folder for the CSV files, figures go to DIR/figures (default: experiments/results "
                        "and docs/figures; with --quick experiments/results_quick)")
    if "recordings" in offered:
        p.add_argument("--recordings", type=Path, metavar="PATH",
                       help="folder with vid2.mat, vid3.mat and vid4.mat (default: data/recordings)")
    if set(offered) & set(SIMULATIONS):
        p.add_argument("--jobs", type=int, default=os.cpu_count(),
                       help="worker processes for the simulations (default: number of CPUs)")
    p.add_argument("--replot", action="store_true",
                   help="redraw the figures from the CSV files in the output folder")
    if only is None:
        p.add_argument("--only", nargs="+", choices=PARTS,
                       help="run (or with --replot redraw) only these parts (default: all)")
    a = p.parse_args(argv)
    chosen = only or a.only
    requested = set(chosen or [])
    parts = set(chosen or PARTS)
    given = getattr(a, "recordings", None)
    folder = given if given is not None else common.RECORDINGS
    if not a.replot:
        if given is not None and "recordings" not in parts:
            p.error("--recordings is given, but --only leaves out recordings")
        if "recordings" in parts:
            missing = [n for n in recordings_2025.FILES if not (folder / n).is_file()]
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
    if "pipeline" in parts:
        pipeline_figure.run()
    if "testsequences" in parts:
        test_sequences.summary(test_sequences.run())
    if "recordings" in parts:
        recordings_2025.run(folder)


if __name__ == "__main__":
    main()
