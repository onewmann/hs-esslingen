"""Accuracy, noise robustness, validity check and before/after comparison.

    python benchmark/run_benchmark.py                 # about 15 min on 4 cores
    python benchmark/run_benchmark.py --quick         # smaller grids, for a check
    python benchmark/run_benchmark.py --recordings ../2D-FFT/data

Writes CSV files to benchmark/results/ and figures to docs/figures/.
``--recordings`` points to the folder with the original camera recordings
(vid2.mat .. vid4.mat, not part of this repository) and adds the analysis of
their frame-to-frame motion.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "python" / "src"))
sys.path.insert(0, str(HERE))

from fftvel import estimate_velocity, phase_correlation, simulate_pulsed_image  # noqa: E402

RESULTS = HERE / "results"
FIGURES = REPO / "docs" / "figures"
N = 256

# colours of the reference palette (categorical slots 1-3, text, surface)
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, SURFACE, MUTED = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb", "#a3a29d"

_texture = None


def _init():
    global _texture
    import matplotlib.image as mpimg
    _texture = mpimg.imread(REPO / "data" / "surface_texture.png").astype(float)


def angle_error(a, b):
    return (a - b + 90) % 180 - 90


def _case(args):
    """One simulated pulsed image -> new estimate (and optionally the original)."""
    texture, speed, angle, L, noise, seed, with_original = args
    a = np.deg2rad(angle)
    v = (speed * np.cos(a), -speed * np.sin(a))
    tex = _texture if texture == "surface" else None
    bp = simulate_pulsed_image(N, v, L, texture=tex, noise=noise, rng=seed)
    e = estimate_velocity(bp, L=L)
    row = dict(texture=texture, speed=speed, angle=angle, L=L, noise=noise, seed=seed,
               est_angle=e.angle_deg, est_speed=e.speed, geba=e.geba, quality=e.quality,
               valid=int(e.valid))
    if speed > 0:
        row["angle_err"] = abs(angle_error(e.angle_deg, angle))
        row["speed_err"] = abs(e.speed / speed - 1)
        row["correct"] = int(row["angle_err"] <= 2 and row["speed_err"] <= 0.05)
    if with_original:
        from original_estimator import original_estimate
        scaled = (bp - bp.min()) / np.ptp(bp)
        beta, vmag = original_estimate(scaled)
        row["orig_angle"] = beta
        row["orig_speed"] = vmag
    return row


def _run(cases, pool):
    return pool.map(_case, cases, chunksize=1)


def _write(name, rows):
    RESULTS.mkdir(exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(RESULTS / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def _read(name):
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


def _style(ax, xlabel, ylabel):
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


def _figure(ncols=1, width=6.4):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, ncols, figsize=(width, 4.2), layout="constrained")
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.86))     # room for title and subtitle
    fig.patch.set_facecolor(SURFACE)
    return fig, np.atleast_1d(axes)


def _title(fig, text, sub=None):
    fig.text(0.012, 0.965, text, ha="left", va="center", color=INK, fontsize=12, fontweight="bold")
    if sub:
        fig.text(0.012, 0.905, sub, ha="left", va="center", color=INK2, fontsize=9)


def _plain_log_ticks(ax, axis="y"):
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
    a = ax.yaxis if axis == "y" else ax.xaxis
    a.set_major_locator(LogLocator(subs=(1, 2, 5)))
    a.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    a.set_minor_formatter(NullFormatter())


def accuracy(pool, quick):
    speeds = [0.5, 1, 2, 4, 8, 16, 32, 64] if quick else \
        [0.3, 0.4, 0.5, 0.7, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]
    angles = range(0, 180, 30 if quick else 15)
    cases = []
    for texture, L, top in (("random", 4, 64), ("surface", 4, 64), ("random", 20, 12)):
        for s in speeds:
            if s > top:
                continue
            for a in angles:
                cases.append((texture, s, a, L, 0.0, 1000 + len(cases), False))
    if not quick:                         # slower motion for 41 exposures
        for s in (0.1, 0.15, 0.2):
            for a in angles:
                cases.append(("random", s, a, 20, 0.0, 1000 + len(cases), False))
    rows = _run(cases, pool)
    _write("accuracy.csv", rows)
    plot_accuracy(rows)
    return rows


def plot_accuracy(rows):
    fig, (ax1, ax2) = _figure(2, 9.6)
    series = [("random", 4, BLUE, "random texture, 9 exposures"),
              ("surface", 4, ORANGE, "camera texture, 9 exposures"),
              ("random", 20, AQUA, "random texture, 41 exposures")]
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
        _plain_log_ticks(ax, "x")
        _plain_log_ticks(ax, "y")
        _style(ax, "true speed (px per exposure step)", ylab)
    ax1.legend(frameon=False, fontsize=9, labelcolor=INK2, loc="upper right")
    ndir = len({r["angle"] for r in rows})
    _title(fig, "Accuracy over the speed range",
           f"256 x 256 px, no noise, worst case over {ndir} directions; slower motion than shown "
           f"is outside the search range and comes out 'not valid'")
    fig.savefig(FIGURES / "accuracy.png", dpi=150, facecolor=SURFACE)


def noise(pool, quick):
    levels = [0, 0.5, 1, 2] if quick else [0, 0.25, 0.5, 0.75, 1, 1.5, 2]
    seeds = 1 if quick else 3
    cases = []
    for texture in ("random", "surface"):
        for nl in levels:
            for s in (2, 8, 24):
                for a in range(0, 180, 30):
                    for k in range(seeds):
                        cases.append((texture, s, a, 4, nl, 5000 + len(cases), False))
    static = []
    for texture in ("random", "surface"):
        for nl in (0, 0.5, 1):
            for k in range(6 if quick else 20):
                static.append((texture, 0.0, 0, 4, nl, 9000 + len(static), False))
    rows = _run(cases, pool)
    srows = _run(static, pool)
    _write("noise.csv", rows)
    _write("static.csv", srows)
    plot_noise(rows, srows)
    return rows, srows


def plot_noise(rows, srows):
    levels = sorted({r["noise"] for r in rows})
    fig, (ax1, ax2) = _figure(2, 9.6)
    for texture, col, label in (("random", BLUE, "random texture"), ("surface", ORANGE, "camera texture")):
        corr, valid = [], []
        for nl in levels:
            sel = [r for r in rows if r["texture"] == texture and r["noise"] == nl]
            corr.append(100 * np.mean([r["correct"] for r in sel]))
            valid.append(100 * np.mean([r["valid"] for r in sel]))
        ax1.plot(levels, corr, "-o", color=col, lw=2, ms=4.5, label=label)
        ax2.plot(levels, valid, "-o", color=col, lw=2, ms=4.5, label=label)
    _style(ax1, "noise / signal (standard deviation)", "correct estimates (%)")
    _style(ax2, "noise / signal (standard deviation)", "marked valid (%)")
    for ax in (ax1, ax2):
        ax.set_ylim(0, 104)
    ax1.legend(frameon=False, fontsize=9, labelcolor=INK2, loc="lower left")
    val = [r for r in rows if r["valid"]]
    wrong_valid = sum(1 - r["correct"] for r in val)
    _title(fig, "Sensor noise",
           f"correct = direction within 2 deg and speed within 5 %, speeds 2, 8 and 24 px per step; "
           f"{wrong_valid} of {len(val)} valid estimates were wrong")
    fig.savefig(FIGURES / "noise.png", dpi=150, facecolor=SURFACE)

    fig, (ax,) = _figure(1, 6.4)
    bins = np.linspace(0, 1, 41)
    groups = [([r["quality"] for r in srows], MUTED, "no motion"),
              ([r["quality"] for r in rows if not r["correct"]], ORANGE, "moving, wrong estimate"),
              ([r["quality"] for r in rows if r["correct"]], BLUE, "moving, correct estimate")]
    ax.hist([g[0] for g in groups], bins=bins, stacked=True, color=[g[1] for g in groups],
            label=[g[2] for g in groups], edgecolor=SURFACE, linewidth=0.6)
    ax.axvline(0.6, color=INK, ls="--", lw=1.2)
    ax.text(0.61, ax.get_ylim()[1] * 0.92, "threshold 0.6", color=INK, fontsize=9)
    _style(ax, "fit quality", "images")
    ax.legend(frameon=False, fontsize=9, labelcolor=INK2, loc="upper left")
    _title(fig, "Fit quality separates usable from unusable images",
           "simulated images with and without motion, all noise levels")
    fig.savefig(FIGURES / "quality.png", dpi=150, facecolor=SURFACE)


def before_after(pool, quick):
    speeds = [2, 8, 32] if quick else [1, 2, 4, 8, 16, 32]
    angles = [20, 110] if quick else [10, 40, 70, 100, 130, 160]
    cases = [("random", s, a, 4, 0.0, 7000 + i, True)
             for i, (s, a) in enumerate((s, a) for s in speeds for a in angles)]
    rows = _run(cases, pool)
    _write("before_after.csv", rows)
    plot_before_after(rows)
    return rows


def plot_before_after(rows):
    fig, (ax,) = _figure(1, 6.4)
    lim = [0.6, 120]
    ax.plot(lim, lim, color=MUTED, lw=1.2, zorder=1)
    ax.scatter([r["speed"] for r in rows], [r["orig_speed"] for r in rows], s=34, color=ORANGE,
               edgecolor=SURFACE, linewidth=1.2, zorder=2, label="original code")
    ax.scatter([r["speed"] for r in rows], [r["est_speed"] for r in rows], s=34, color=BLUE,
               edgecolor=SURFACE, linewidth=1.2, zorder=3, label="rewrite")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    _plain_log_ticks(ax, "x")
    _plain_log_ticks(ax, "y")
    _style(ax, "true speed (px per exposure step)", "estimated speed")
    ax.legend(frameon=False, fontsize=9, labelcolor=INK2, loc="upper left")
    _title(fig, "Speed estimate before and after the rewrite",
           f"same {len(rows)} simulated pulsed images (9 exposures, 256 x 256 px)")
    fig.savefig(FIGURES / "before_after.png", dpi=150, facecolor=SURFACE)


def pipeline_figure():
    import matplotlib.image as mpimg
    from fftvel.plotting import plot_estimate

    tex = mpimg.imread(REPO / "data" / "surface_texture.png").astype(float)
    bp, frames = simulate_pulsed_image(N, (5, -3), 4, texture=tex, rng=3, return_frames=True)
    e = estimate_velocity(bp, L=4, diagnostics=True)
    fig = plot_estimate(bp, e, frame=frames[:, :, 4])
    fig.patch.set_facecolor("white")
    fig.savefig(FIGURES / "pipeline.png", dpi=110, facecolor="white")


def testsequences():
    from scipy.io import loadmat

    data = loadmat(REPO / "data" / "testsequences.mat")
    rows = []
    fig, axes = _figure(3, 11.5)
    for ax, (name, title) in zip(axes, (("model", "model: (5, 5) px/frame"),
                                        ("tyre", "tyre: (4, 4) px/frame"),
                                        ("accelerating", "accelerating: vx 1.8-8.2, vy 3"))):
        s = data[name].astype(float) / 255
        nf = s.shape[2]
        shift = np.array([phase_correlation(s[:, :, f], s[:, :, f + 1])[:2] for f in range(nf - 1)])
        starts = np.arange(nf - 8)
        est, ref, valid = [], [], []
        for st in starts:
            last = st + 8
            border = int(np.ceil(np.abs(shift[:last].sum(axis=0)).max()))
            e = estimate_velocity(s[border:, border:, st:last + 1].mean(axis=2), L=4)
            r = float(np.hypot(*shift[st:last].mean(axis=0)))
            est.append(e.speed)
            valid.append(e.valid)
            ref.append(r)
            rows.append(dict(sequence=name, first_frame=st + 1, est_speed=e.speed, ref_speed=r,
                             est_angle=e.angle_deg, quality=e.quality, valid=int(e.valid)))
        est, ref, valid = np.array(est), np.array(ref), np.array(valid)
        ax.plot(starts + 1, ref, color=INK2, lw=2, label="phase correlation (reference)")
        ax.plot(starts[valid] + 1, est[valid], "o", color=BLUE, ms=6, mec=SURFACE, mew=1.2,
                label="single-image estimate")
        if (~valid).any():
            ax.plot(starts[~valid] + 1, est[~valid], "o", color=ORANGE, ms=6, mfc="none", label="not valid")
        ax.set_ylim(0, 10)
        ax.set_title(title, color=INK2, fontsize=10, loc="left")
        _style(ax, "first frame of the 9-frame window", "speed (px/frame)" if ax is axes[0] else "")
    axes[0].legend(frameon=False, fontsize=9, labelcolor=INK2, loc="lower left")
    _title(fig, "Test sequences with known motion",
           "every window of 9 frames, cropped to the area without zero-filled borders")
    fig.savefig(FIGURES / "testsequences.png", dpi=150, facecolor=SURFACE)
    _write("testsequences.csv", rows)
    return rows


def recordings(folder: Path):
    from scipy.io import loadmat

    rows = []
    fig, (ax,) = _figure(1, 6.4)
    for k, col in zip((2, 3, 4), (BLUE, ORANGE, AQUA)):
        rec = loadmat(folder / f"vid{k}.mat")["recording"][:, :, 0, :].astype(float)
        H, W, F = rec.shape
        c = rec[H // 2 - 512:H // 2 + 512, W // 2 - 512:W // 2 + 512, :]
        steps = []
        for i in range(F - 1):
            dx, dy, _ = phase_correlation(c[:, :, i], c[:, :, i + 1])
            steps.append(np.hypot(dx, dy))
            rows.append(dict(recording=f"vid{k}", step=i + 1, dx=dx, dy=dy))
        for s in (0, 1):
            e = estimate_velocity(c[:, :, s:s + 9].mean(axis=2), L=4)
            rows.append(dict(recording=f"vid{k}", window=f"{s + 1}-{s + 9}", est_speed=e.speed,
                             quality=e.quality, valid=int(e.valid)))
        ax.plot(range(1, F), steps, "-o", color=col, lw=2, ms=4.5, label=f"vid{k}")
    _style(ax, "frame pair", "displacement between frames (px)")
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK2, loc="upper left")
    _title(fig, "Original recordings: steps are not equal",
           "phase correlation of consecutive frames, central 1024 x 1024 px")
    fig.savefig(FIGURES / "recordings.png", dpi=150, facecolor=SURFACE)
    _write("recordings.csv", rows)
    return rows


def summary(acc, noi, stat, ba):
    print("\nAccuracy (no noise):")
    for texture, L in (("random", 4), ("surface", 4), ("random", 20)):
        sel = [r for r in acc if r["texture"] == texture and r["L"] == L]
        good = sorted({r["speed"] for r in sel
                       if all(x["angle_err"] <= 0.5 and x["speed_err"] <= 0.01
                              for x in sel if x["speed"] == r["speed"])})
        if good:
            print(f"  {texture:7s} L={L:2d}: within 0.5 deg and 1 % for {good[0]}..{good[-1]} px/step "
                  f"(median errors {np.median([r['angle_err'] for r in sel if r['speed'] in good]):.3f} deg, "
                  f"{100 * np.median([r['speed_err'] for r in sel if r['speed'] in good]):.3f} %)")
    print("Noise (share correct / share correct among valid):")
    for nl in sorted({r["noise"] for r in noi}):
        sel = [r for r in noi if r["noise"] == nl]
        val = [r for r in sel if r["valid"]]
        print(f"  noise {nl:4.2f}: {100 * np.mean([r['correct'] for r in sel]):5.1f} % correct, "
              f"{100 * len(val) / len(sel):5.1f} % valid, "
              f"{100 * np.mean([r['correct'] for r in val]) if val else float('nan'):5.1f} % of valid correct")
    print(f"Static scenes marked valid: {sum(r['valid'] for r in stat)} of {len(stat)}, "
          f"max quality {max(r['quality'] for r in stat):.3f}")
    oe = [abs(r["orig_speed"] / r["speed"] - 1) for r in ba]
    ne = [abs(r["est_speed"] / r["speed"] - 1) for r in ba]
    print(f"Before/after, median speed error: original {100 * np.median(oe):.0f} %, rewrite {100 * np.median(ne):.3f} %")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--recordings", type=Path)
    p.add_argument("--jobs", type=int, default=os.cpu_count())
    p.add_argument("--replot", action="store_true",
                   help="redraw the simulation figures from benchmark/results/*.csv")
    p.add_argument("--only", nargs="+", choices=["simulation", "pipeline", "testsequences", "recordings"],
                   help="run only these parts (default: all)")
    a = p.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    FIGURES.mkdir(parents=True, exist_ok=True)
    parts = set(a.only or ["simulation", "pipeline", "testsequences", "recordings"])
    if a.replot:
        acc, noi, stat, ba = (_read(n) for n in ("accuracy.csv", "noise.csv", "static.csv",
                                                  "before_after.csv"))
        plot_accuracy(acc)
        plot_noise(noi, stat)
        plot_before_after(ba)
        summary(acc, noi, stat, ba)
        return
    if "simulation" in parts:
        with Pool(a.jobs, initializer=_init) as pool:
            acc = accuracy(pool, a.quick)
            noi, stat = noise(pool, a.quick)
            ba = before_after(pool, a.quick)
        summary(acc, noi, stat, ba)
    if "pipeline" in parts:
        pipeline_figure()
    if "testsequences" in parts:
        ts = testsequences()
        for name in ("model", "tyre", "accelerating"):
            sel = [r for r in ts if r["sequence"] == name]
            val = [r for r in sel if r["valid"]]
            err = [abs(r["est_speed"] - r["ref_speed"]) for r in val]
            print(f"{name:13s} {len(val)}/{len(sel)} windows valid, max |error| {max(err):.3f} px/frame")
    if "recordings" in parts and a.recordings:
        recordings(a.recordings)


if __name__ == "__main__":
    main()
