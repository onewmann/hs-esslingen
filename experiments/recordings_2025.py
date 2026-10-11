"""Camera recordings of the 2025 student project (data/recordings): the
frame-to-frame motion in vid2.mat, vid3.mat and vid4.mat by phase
correlation, and the single-image estimate of two 9-frame averages.

    python experiments/recordings_2025.py [--recordings PATH] [--out DIR] [--replot]

does the same as ``python experiments/run.py --only recordings ...`` (see
run.py). Writes recordings.csv and recordings.png. ``--quick`` only selects
the output folder; this study has no smaller variant.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

import common
from common import AQUA, BLUE, INK2, ORANGE
from fftvel import estimate_velocity, phase_correlation

CSV = ["recordings.csv"]
FILES = [f"vid{k}.mat" for k in (2, 3, 4)]


def run(folder: Path):
    from scipy.io import loadmat

    rows = []
    for k in (2, 3, 4):
        rec = loadmat(folder / f"vid{k}.mat")["recording"][:, :, 0, :].astype(float)
        H, W, F = rec.shape
        c = rec[H // 2 - 512:H // 2 + 512, W // 2 - 512:W // 2 + 512, :]
        for i in range(F - 1):
            dx, dy, _ = phase_correlation(c[:, :, i], c[:, :, i + 1])
            rows.append(dict(recording=f"vid{k}", step=i + 1, dx=dx, dy=dy))
        for s in (0, 1):
            e = estimate_velocity(c[:, :, s:s + 9].mean(axis=2), L=4)
            rows.append(dict(recording=f"vid{k}", window=f"{s + 1}-{s + 9}", est_speed=e.speed,
                             quality=e.quality, isotropic=int(e.isotropic), valid=int(e.valid)))
    common.write_csv("recordings.csv", rows)
    plot(rows)
    return rows


def plot(rows):
    fig, (ax,) = common.figure(1, 6.4)
    for name, col in zip(("vid2", "vid3", "vid4"), (BLUE, ORANGE, AQUA)):
        sel = [r for r in rows if r["recording"] == name and r.get("step", "") != ""]
        ax.plot([r["step"] for r in sel], [np.hypot(r["dx"], r["dy"]) for r in sel], "-o",
                color=col, lw=2, ms=4.5, label=name)
    common.style(ax, "frame pair", "displacement between frames (px)")
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK2, loc="upper left")
    common.title(fig, "Original recordings: steps are not equal",
                 "phase correlation of consecutive frames, central 1024 x 1024 px")
    common.save(fig, "recordings.png")


def replot():
    plot(common.read_csv("recordings.csv"))


def main(argv=None):
    import run as cli
    cli.main(argv, only=["recordings"], description=__doc__)


if __name__ == "__main__":
    main()
