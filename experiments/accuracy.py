"""Accuracy over the speed range: simulated images without noise.

    python experiments/accuracy.py [--quick] [--out DIR] [--replot] [--jobs N]

does the same as ``python experiments/run.py --only accuracy ...`` (see run.py).
Writes accuracy.csv and accuracy.png.
"""

from __future__ import annotations

import numpy as np

import common
from common import AQUA, BLUE, INK2, ORANGE

CSV = ["accuracy.csv"]


def run(pool, quick):
    speeds = [0.5, 1, 2, 4, 8, 16, 32, 64] if quick else \
        [0.3, 0.4, 0.5, 0.7, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]
    angles = range(0, 180, 30 if quick else 15)
    cases = []
    for texture, L, top in (("random", 4, 64), ("surface", 4, 64), ("random", 20, 12)):
        for s in speeds:
            if s > top:
                continue
            for a in angles:
                cases.append((texture, s, a, L, 0.0, 1000 + len(cases)))
    if not quick:                         # slower motion for 41 exposures
        for s in (0.1, 0.15, 0.2):
            for a in angles:
                cases.append(("random", s, a, 20, 0.0, 1000 + len(cases)))
    rows = common.run_cases(cases, pool)
    common.write_csv("accuracy.csv", rows)
    plot(rows)
    return rows


def plot(rows):
    fig, (ax1, ax2) = common.figure(2, 9.6)
    series = [("random", 4, BLUE, "random texture, 9 flashes"),
              ("surface", 4, ORANGE, "camera texture, 9 flashes"),
              ("random", 20, AQUA, "random texture, 41 flashes")]
    for texture, L, col, label in series:
        sel = [r for r in rows if r["texture"] == texture and r["L"] == L]
        sp = sorted(s for s in {r["speed"] for r in sel}
                    if all(r["valid"] for r in sel if r["speed"] == s))
        ae = [max(r["angle_err"] for r in sel if r["speed"] == s) for s in sp]
        se = [100 * max(r["speed_err"] for r in sel if r["speed"] == s) for s in sp]
        ax1.plot(sp, np.maximum(ae, 1e-3), "-o", color=col, lw=2, ms=4.5, label=label)
        ax2.plot(sp, np.maximum(se, 1e-3), "-o", color=col, lw=2, ms=4.5, label=label)
    for ax, ylab in ((ax1, "largest direction error (deg)"), (ax2, "largest speed error (%)")):
        ax.set_xscale("log")
        ax.set_yscale("log")
        common.plain_log_ticks(ax, "x")
        common.plain_log_ticks(ax, "y")
        common.style(ax, "true speed (px per flash)", ylab)
    ax1.legend(frameon=False, fontsize=9, labelcolor=INK2, loc="upper right")
    ndir = len({r["angle"] for r in rows})
    common.title(fig, "Accuracy over the speed range",
                 f"256 x 256 px, no noise, worst case over {ndir} directions; slower motion than shown "
                 f"is outside the measurable range and comes out 'not valid'")
    common.save(fig, "accuracy.png")


def summary(rows):
    """Print the numbers the README quotes."""
    print("\nAccuracy (no noise):")
    for texture, L in (("random", 4), ("surface", 4), ("random", 20)):
        sel = [r for r in rows if r["texture"] == texture and r["L"] == L]
        good = sorted({r["speed"] for r in sel
                       if all(x["angle_err"] <= 0.5 and x["speed_err"] <= 0.01
                              for x in sel if x["speed"] == r["speed"])})
        if good:
            print(f"  {texture:7s} L={L:2d}: within 0.5 deg and 1 % for {good[0]}..{good[-1]} px/step "
                  f"(median errors {np.median([r['angle_err'] for r in sel if r['speed'] in good]):.3f} deg, "
                  f"{100 * np.median([r['speed_err'] for r in sel if r['speed'] in good]):.3f} %)")


def replot():
    rows = common.read_csv("accuracy.csv")
    plot(rows)
    summary(rows)


def main(argv=None):
    import run as cli
    cli.main(argv, only=["accuracy"], description=__doc__)


if __name__ == "__main__":
    main()
