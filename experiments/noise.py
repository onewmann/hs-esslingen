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


TEXTURES = (("surface", "Camera texture"), ("random", "Random texture"))     # camera texture first


def threshold_at(n):
    """Validity threshold of estimate_velocity for n-by-n images; it depends on the size, not on the content."""
    from fftvel import estimate_velocity, simulate_pulsed_image
    return estimate_velocity(simulate_pulsed_image(n, (0, 0), 4, rng=0), L=4).threshold


def plot(rows, srows):
    levels = sorted({r["noise"] for r in rows})
    per_level = sorted({sum(1 for r in rows if r["texture"] == t and r["noise"] == nl)
                        for t, _ in TEXTURES for nl in levels})
    flashes = sorted({2 * int(r["L"]) + 1 for r in rows})
    fig, axes = common.figure(2, 9.6)
    for ax, (texture, name) in zip(axes, TEXTURES):
        corr, valid, wrong = [], [], []
        for nl in levels:
            sel = [r for r in rows if r["texture"] == texture and r["noise"] == nl]
            corr.append(100 * np.mean([r["correct"] for r in sel]))
            valid.append(100 * np.mean([r["valid"] for r in sel]))
            wrong.append(100 * np.mean([r["valid"] * (1 - r["correct"]) for r in sel]))
        ax.plot(levels, corr, "-o", color=BLUE, lw=2, ms=4.5, label="correct estimate")
        ax.plot(levels, valid, "--", color=INK2, lw=1.6, label="marked valid")
        ax.plot(levels, wrong, "-o", color=ORANGE, lw=2, ms=4.5, label="valid but wrong")
        sel = [r for r in rows if r["texture"] == texture]
        val = [r for r in sel if r["valid"]]
        ax.set_title(f"{name}: {len(sel)} images, {int(sum(1 - r['correct'] for r in val))} of {len(val)} "
                     "valid estimates wrong", color=INK2, fontsize=10, loc="left")
        common.style(ax, "noise / signal (standard deviation)", "images (%)" if ax is axes[0] else "")
        ax.set_ylim(-2, 104)
    axes[0].legend(frameon=False, fontsize=9, labelcolor=INK2, loc="center left")
    common.title(fig, "Sensor noise",
                 f"in simulation: {_join(per_level)} images per noise level and texture ({common.N} \u00d7 {common.N} px, "
                 f"{_join(flashes)} flashes, 1 to 32 px per flash); correct = within 2 deg and 5 %")
    common.save(fig, "noise.png")

    fig, (ax,) = common.figure(1, 9.6)
    bins = np.linspace(0, 1, 41)
    moving_wrong = [r["quality"] for r in rows if not r["correct"]]
    groups = [([r["quality"] for r in srows], MUTED, f"no motion ({len(srows)} images)"),
              (moving_wrong, ORANGE, f"moving, wrong estimate ({len(moving_wrong)} images)"),
              ([r["quality"] for r in rows if r["correct"]], BLUE,
               f"moving, correct estimate ({int(sum(r['correct'] for r in rows))} images)")]
    ax.hist([g[0] for g in groups], bins=bins, stacked=True, color=[g[1] for g in groups],
            label=[g[2] for g in groups], edgecolor=SURFACE, linewidth=0.6)
    thr = threshold_at(common.N)
    ax.axvline(thr, color=INK, ls="--", lw=1.2)
    ax.text(thr + 0.01, 0.96, f"threshold {thr:.2g}\n({common.N} px images)", color=INK, fontsize=9, va="top",
            transform=ax.get_xaxis_transform())
    common.style(ax, "fit quality", "images")
    ax.set_xlim(0, 1.02)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK2, loc="upper left")
    over_wrong = sum(q >= thr for q in moving_wrong)
    valid_wrong = int(sum(r["valid"] * (1 - r["correct"]) for r in rows))
    over_static = sum(r["quality"] >= thr for r in srows)
    common.title(fig, "Fit quality of moving images and of images without motion",
                 f"in simulation, {common.N} \u00d7 {common.N} px, all noise levels: {over_wrong} of {len(moving_wrong)} "
                 f"wrong estimates reach the threshold ({valid_wrong} marked valid), {over_static} of {len(srows)} "
                 "images without motion")
    common.save(fig, "quality.png")


def _join(values):
    """'60' for one value, '12 to 60' for several (the --quick grid has the same count everywhere)."""
    values = [f"{v:g}" for v in values]
    return values[0] if len(values) == 1 else f"{values[0]} to {values[-1]}"


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
