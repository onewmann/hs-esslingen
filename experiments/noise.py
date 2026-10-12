"""Sensor noise and the validity check: random velocities at several noise
levels, and images without motion.

    python experiments/noise.py [--quick] [--out DIR] [--replot] [--jobs N]

does the same as ``python experiments/run.py --only noise ...`` (see run.py).
Writes noise.csv, static.csv, noise.png and quality.png.
"""

from __future__ import annotations

import numpy as np

import common
from common import BLUE, INK, INK2, MUTED, ORANGE, SURFACE

CSV = ["noise.csv", "static.csv"]


def run(pool, quick):
    levels = [0, 0.5, 1, 2] if quick else [0, 0.25, 0.5, 0.75, 1, 1.5, 2]
    per_level = 12 if quick else 60
    rng = np.random.default_rng(2025)
    cases = []
    for texture in ("random", "surface"):
        for nl in levels:
            for _ in range(per_level):        # random velocity: speed log-uniform 1..32, any direction
                s = float(np.exp(rng.uniform(np.log(1), np.log(32))))
                a = float(rng.uniform(0, 180))
                cases.append((texture, s, a, 4, nl, 5000 + len(cases)))
    static = []
    for texture in ("random", "surface"):
        for nl in (0, 0.5, 1):
            for k in range(6 if quick else 20):
                static.append((texture, 0.0, 0, 4, nl, 9000 + len(static)))
    rows = common.run_cases(cases, pool)
    srows = common.run_cases(static, pool)
    common.write_csv("noise.csv", rows)
    common.write_csv("static.csv", srows)
    plot(rows, srows)
    return rows, srows


def plot(rows, srows):
    levels = sorted({r["noise"] for r in rows})
    fig, (ax1, ax2) = common.figure(2, 9.6)
    for texture, col, label in (("random", BLUE, "random texture"), ("surface", ORANGE, "camera texture")):
        corr, valid = [], []
        for nl in levels:
            sel = [r for r in rows if r["texture"] == texture and r["noise"] == nl]
            corr.append(100 * np.mean([r["correct"] for r in sel]))
            valid.append(100 * np.mean([r["valid"] for r in sel]))
        ax1.plot(levels, corr, "-o", color=col, lw=2, ms=4.5, label=label)
        ax2.plot(levels, valid, "-o", color=col, lw=2, ms=4.5, label=label)
    common.style(ax1, "noise / signal (standard deviation)", "correct estimates (%)")
    common.style(ax2, "noise / signal (standard deviation)", "marked valid (%)")
    for ax in (ax1, ax2):
        ax.set_ylim(0, 104)
    ax1.legend(frameon=False, fontsize=9, labelcolor=INK2, loc="lower left")
    val = [r for r in rows if r["valid"]]
    wrong_valid = int(sum(1 - r["correct"] for r in val))
    common.title(fig, "Sensor noise",
                 f"correct = direction within 2 deg and speed within 5 %; random velocities, 1 to 32 px per step; "
                 f"{wrong_valid} of {len(val)} valid estimates were wrong")
    common.save(fig, "noise.png")

    fig, (ax,) = common.figure(1, 6.4)
    bins = np.linspace(0, 1, 41)
    groups = [([r["quality"] for r in srows], MUTED, "no motion"),
              ([r["quality"] for r in rows if not r["correct"]], ORANGE, "moving, wrong estimate"),
              ([r["quality"] for r in rows if r["correct"]], BLUE, "moving, correct estimate")]
    ax.hist([g[0] for g in groups], bins=bins, stacked=True, color=[g[1] for g in groups],
            label=[g[2] for g in groups], edgecolor=SURFACE, linewidth=0.6)
    ax.axvline(0.6, color=INK, ls="--", lw=1.2)
    ax.text(0.61, ax.get_ylim()[1] * 0.92, "threshold 0.6", color=INK, fontsize=9)
    common.style(ax, "fit quality", "images")
    ax.legend(frameon=False, fontsize=9, labelcolor=INK2, loc="upper left")
    common.title(fig, "Fit quality separates usable from unusable images",
                 "simulated images with and without motion, all noise levels")
    common.save(fig, "quality.png")


def summary(rows=None, srows=None):
    """Print the numbers the README quotes, for the rows that are given."""
    if rows is not None:
        print("Noise (share correct / share correct among valid):")
        for nl in sorted({r["noise"] for r in rows}):
            sel = [r for r in rows if r["noise"] == nl]
            val = [r for r in sel if r["valid"]]
            print(f"  noise {nl:4.2f}: {100 * np.mean([r['correct'] for r in sel]):5.1f} % correct, "
                  f"{100 * len(val) / len(sel):5.1f} % valid, "
                  f"{100 * np.mean([r['correct'] for r in val]) if val else float('nan'):5.1f} % of valid correct")
        if rows and "blur_like" in rows[0]:          # not in the CSV files of runs before the blur check
            blur_summary(rows)
    if srows is not None:
        print(f"Static scenes marked valid: {int(sum(r['valid'] for r in srows))} of {len(srows)}, "
              f"max quality {max(r['quality'] for r in srows):.3f}")


def blur_summary(rows):
    """Print how many valid estimates are blur-like, so that ``reject_blur`` would make them not valid."""
    slow = common.blur_limit()
    val = [r for r in rows if r["valid"]]
    groups = (("correct", [r for r in val if r["correct"]]), ("wrong", [r for r in val if not r["correct"]]))
    parts = []
    for name, sel in groups:
        below = [r for r in sel if r["speed"] < slow]
        parts.append(f"{int(sum(r['blur_like'] for r in sel))} of {len(sel)} {name} "
                     f"({int(sum(r['blur_like'] for r in below))} of {len(below)} below {slow:.2f} px)")
    print("Blur-like valid estimates (not valid with reject_blur): " + ", ".join(parts))


def replot():
    rows, srows = common.read_csv("noise.csv"), common.read_csv("static.csv")
    plot(rows, srows)
    summary(rows, srows)


def main(argv=None):
    import run as cli
    cli.main(argv, only=["noise"], description=__doc__)


if __name__ == "__main__":
    main()
