"""Speed estimate of the 2025 code and of the rework on the same simulated
pulsed images. The 2025 code is the estimator of the 2025 student project,
ported in original_2025.py.

    python experiments/history/before_after.py [--quick] [--out DIR] [--replot] [--jobs N]

does the same as ``python experiments/run.py --only before_after ...`` (see
run.py). Writes before_after.csv and before_after.png.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # experiments/, when run directly

import numpy as np  # noqa: E402

import common  # noqa: E402
from common import BLUE, INK2, MUTED, ORANGE, SURFACE  # noqa: E402

CSV = ["before_after.csv"]


def _case(case):
    """One simulated pulsed image -> CSV row with the estimate of the rework and of the 2025 code."""
    from history.original_2025 import original_estimate

    bp = common.simulate_case(case)
    row = common.evaluate_case(case, bp)
    scaled = (bp - bp.min()) / np.ptp(bp)
    beta, vmag = original_estimate(scaled)
    row["orig_angle"] = beta
    row["orig_speed"] = vmag
    return row


def run(pool, quick):
    speeds = [2, 8, 32] if quick else [1, 2, 4, 8, 16, 32]
    angles = [20, 110] if quick else [10, 40, 70, 100, 130, 160]
    cases = [("random", s, a, 4, 0.0, 7000 + i)
             for i, (s, a) in enumerate((s, a) for s in speeds for a in angles)]
    rows = common.run_cases(cases, pool, _case)
    common.write_csv("before_after.csv", rows)
    plot(rows)
    return rows


def plot(rows):
    fig, (ax,) = common.figure(1, 6.4)
    lim = [0.6, 120]
    ax.plot(lim, lim, color=MUTED, lw=1.2, zorder=1)
    ax.scatter([r["speed"] for r in rows], [r["orig_speed"] for r in rows], s=34, color=ORANGE,
               edgecolor=SURFACE, linewidth=1.2, zorder=2, label="2025 code")
    ax.scatter([r["speed"] for r in rows], [r["est_speed"] for r in rows], s=34, color=BLUE,
               edgecolor=SURFACE, linewidth=1.2, zorder=3, label="rework")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    common.plain_log_ticks(ax, "x")
    common.plain_log_ticks(ax, "y")
    common.style(ax, "true speed (px per exposure step)", "estimated speed")
    ax.legend(frameon=False, fontsize=9, labelcolor=INK2, loc="lower right")
    common.title(fig, "Speed estimate of the 2025 code and the rework",
                 f"same {len(rows)} simulated pulsed images (9 exposures, 256 x 256 px)")
    common.save(fig, "before_after.png")


def summary(rows):
    """Print the numbers the README quotes."""
    oe = [abs(r["orig_speed"] / r["speed"] - 1) for r in rows]
    ne = [abs(r["est_speed"] / r["speed"] - 1) for r in rows]
    print(f"Before/after, median speed error: 2025 code {100 * np.median(oe):.0f} %, "
          f"rework {100 * np.median(ne):.3f} %")


def replot():
    rows = common.read_csv("before_after.csv")
    plot(rows)
    summary(rows)


def main(argv=None):
    import run as cli
    cli.main(argv, only=["before_after"], description=__doc__)


if __name__ == "__main__":
    main()
