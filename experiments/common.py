"""Shared parts of the experiments: paths, output folders, plot style, the
simulated test case and the CSV files.

The study modules write through ``common.RESULTS`` and ``common.FIGURES``,
which ``run.py`` sets once with :func:`set_output` before any study runs.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "python" / "src"))

from fftvel import angle_difference, estimate_velocity, simulate_pulsed_image  # noqa: E402
from fftvel.core import DEFAULT_PAD_FACTOR, RFRAC  # noqa: E402

COMMITTED = HERE / "results"           # full-run results quoted in the README
QUICK = HERE / "results_quick"         # --quick results, ignored by git
README_FIGURES = REPO / "docs" / "figures"
RECORDINGS = REPO / "data" / "recordings"
RESULTS = COMMITTED                    # set by set_output()
FIGURES = README_FIGURES               # set by set_output()
N = 256

# colours of the reference palette (categorical slots 1-3, text, surface)
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, SURFACE, MUTED = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb", "#a3a29d"


def set_output(out: Path) -> None:
    """Write the CSV files to ``out`` and the figures to out/figures.

    For the committed results (experiments/results) the figures go to
    docs/figures, which the README shows.
    """
    global RESULTS, FIGURES
    RESULTS = out
    FIGURES = README_FIGURES if out == COMMITTED else out / "figures"


# Simulated cases

_texture = None


def init_worker():
    """Pool initializer: load the camera texture once in each worker process."""
    global _texture
    import matplotlib.image as mpimg
    _texture = mpimg.imread(REPO / "data" / "surface_texture.png").astype(float)


def simulate_case(case):
    """Pulsed image of one case ``(texture, speed, angle, L, noise, seed)``."""
    texture, speed, angle, L, noise, seed = case
    a = np.deg2rad(angle)
    v = (speed * np.cos(a), -speed * np.sin(a))
    tex = _texture if texture == "surface" else None
    return simulate_pulsed_image(N, v, L, texture=tex, noise=noise, rng=seed)


def evaluate_case(case, bp):
    """CSV row with the case and the estimate for its pulsed image ``bp``."""
    texture, speed, angle, L, noise, seed = case
    e = estimate_velocity(bp, L=L)
    row = dict(texture=texture, speed=speed, angle=angle, L=L, noise=noise, seed=seed,
               est_angle=e.angle_deg, est_speed=e.speed, geba=e.geba, quality=e.quality,
               quality_perp=e.quality_perp, isotropic=int(e.isotropic), valid=int(e.valid),
               quality_blur=e.quality_blur, blur_length=e.blur_length, blur_like=int(e.blur_like))
    if speed > 0:
        row["angle_err"] = abs(angle_difference(e.angle_deg, angle))
        row["speed_err"] = abs(e.speed / speed - 1)
        row["correct"] = int(row["angle_err"] <= 2 and row["speed_err"] <= 0.05)
    return row


def blur_limit(n=N, pad_factor=DEFAULT_PAD_FACTOR):
    """Speed (px per step) below which estimate_velocity can call a result blur-like: npad / max(profile_j).

    About 2.24 px for 256 px images (512 / 229); see ``Estimate.blur_like``.
    """
    npad = pad_factor * n
    return npad / math.floor(RFRAC * (npad // 2 - 1))


def run_case(case):
    """One simulated pulsed image -> CSV row with the estimate."""
    return evaluate_case(case, simulate_case(case))


def run_cases(cases, pool, worker=run_case):
    """Rows of all cases in their order; ``worker`` must be a top-level function (Windows spawn)."""
    return pool.map(worker, cases, chunksize=1)


# CSV files

def write_csv(name, rows):
    RESULTS.mkdir(parents=True, exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(RESULTS / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def read_csv(name):
    out = []
    with open(RESULTS / name, newline="") as f:
        for r in csv.DictReader(f):
            row = {}
            for k, v in r.items():
                try:
                    row[k] = float(v)
                except ValueError:
                    row[k] = v
            out.append(row)
    return out


def read_if(name):
    """Rows of an earlier run in the same output folder, or None if there is none."""
    return read_csv(name) if (RESULTS / name).is_file() else None


def missing_csv(names):
    """Those of the given CSV files that are not in the output folder."""
    return [n for n in names if not (RESULTS / n).is_file()]


# Plot style

def style(ax, xlabel, ylabel):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.set_xlabel(xlabel, color=INK2, fontsize=10)
    ax.set_ylabel(ylabel, color=INK2, fontsize=10)


def figure(ncols=1, width=6.4):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, ncols, figsize=(width, 4.2), layout="constrained")
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.86))     # room for title and subtitle
    fig.patch.set_facecolor(SURFACE)
    return fig, np.atleast_1d(axes)


def title(fig, text, sub=None):
    fig.text(0.012, 0.965, text, ha="left", va="center", color=INK, fontsize=12, fontweight="bold")
    if sub:
        fig.text(0.012, 0.905, sub, ha="left", va="center", color=INK2, fontsize=9)


def plain_log_ticks(ax, axis="y"):
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
    a = ax.yaxis if axis == "y" else ax.xaxis
    a.set_major_locator(LogLocator(subs=(1, 2, 5)))
    a.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    a.set_minor_formatter(NullFormatter())


def save(fig, name, dpi=150, facecolor=SURFACE):
    """Write ``fig`` to the figure folder as ``name``."""
    fig.savefig(FIGURES / name, dpi=dpi, facecolor=facecolor)
