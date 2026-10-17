"""Compare the estimated directions of the moving graph-paper video with the grid lines and the steps.

    python data/measurements/2026-10-17_phone_grid/directions.py [folder]

Reads results.csv (from experiments/real_measurements.py) and grid_lines.csv
(from extract_frames.py) in folder, by default the folder of this file. For the valid windows with even
steps it prints how far the estimate is from the mean step (the reference) and
from the grid lines, and the direction of the motion relative to the lines.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def wrap(a):
    """Angle difference in deg, wrapped to [-90, 90) (directions modulo 180)."""
    return (np.asarray(a, float) + 90) % 180 - 90


def main(folder: Path = HERE, tol=2.0):
    with open(folder / "results.csv", newline="", encoding="utf-8") as f:
        res = {r["file"]: r for r in csv.DictReader(f)}
    with open(folder / "grid_lines.csv", newline="", encoding="utf-8") as f:
        lines = {r["file"]: float(r["line_angle_deg"]) for r in csv.DictReader(f)}
    win = [r for name, r in res.items() if name.startswith("frames_moving") and r["valid"] == "1"
           and r["uneven_steps"] == "0"]
    est = np.array([float(r["angle_deg"]) for r in win])
    step = np.array([float(r["ref_angle_deg"]) for r in win])
    grid = np.array([lines[r["file"]] for r in win])
    off_step = wrap(est - step)
    off = np.abs(off_step) > tol
    near = np.abs(wrap(est - grid)) < np.abs(wrap(step - grid))
    rel = wrap(step - grid)
    print(f"{len(win)} valid windows with even steps; grid lines at {grid.min():.1f} to {grid.max():.1f} deg")
    print(f"motion relative to the lines: {rel.min():.1f} to {rel.max():.1f} deg")
    print(f"estimate more than {tol:g} deg from the step: {off.sum()}; of these, closer to the lines than the step: "
          f"{(off & near).sum()}, largest distance from the lines {np.abs(wrap(est - grid))[off].max():.2f} deg")
    print(f"all windows: estimate closer to the lines than the step in {near.sum()} of {len(win)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]) if len(sys.argv) > 1 else HERE))
