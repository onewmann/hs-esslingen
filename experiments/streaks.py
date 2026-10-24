"""Static textures of streaks: how often a texture whose correlation length
along one axis is many times that across it passes the validity check while
nothing moves, and whether the same texture is measured correctly in motion.

    python experiments/streaks.py [--quick] [--out DIR] [--replot] [--jobs N]
    python experiments/run.py --only streaks

The isotropy check compares the fit along the direction found with the fit
across it. A texture that is itself anisotropic, such as grain, brushing or
fibres, has a spectrum that falls faster across the streaks than along them,
so the two fits differ as they would for motion, and the check cannot tell
the two apart.

Each static image is a 256 px crop of a 512 px random texture shaped like
fftvel.make_texture, but with two correlation lengths: ``across`` px along x
and ``along`` px along y, so the texture consists of streaks along y. Nothing
moves. The estimate uses the defaults (K = 9); ``valid_strict`` is the result
with ``reject_blur`` and ``clip_outliers`` set. A texture with the same
length in both directions is the control. The moving part shifts a 1024 px
texture of streaks by 1 to 8 px per step at four angles, as
simulate_pulsed_image does, and counts correct estimates (direction within
2 deg, speed within 5 %). Writes streaks.csv (one row per image, ``speed`` 0
for the static ones) and prints the counts; there is no figure. About 6 min
of CPU time, 2 min on 4 cores; the seeds fix every image. --quick takes 3
images per pair of lengths and two speeds.
"""

from __future__ import annotations

import argparse
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import common  # noqa: E402
from fftvel import angle_difference, estimate_velocity, simulate_pulsed_image  # noqa: E402
from fftvel.core import crop_square  # noqa: E402

CSV = ["streaks.csv"]
ACROSS = (1.0, 2.0, 4.0)          # correlation length across the streaks, px
ALONG = (8.0, 16.0, 32.0, 64.0)   # correlation length along the streaks, px
SEEDS = 12                        # static images per pair of lengths (--quick: 3)
MOVING = ((1.0, 32.0), (4.0, 64.0))               # pairs of lengths of the moving part
ANGLES = (0.0, 30.0, 90.0, 135.0)                 # deg, direction of the motion
SPEEDS = (1.0, 3.0, 4.5, 8.0)                     # px per step (--quick: 1 and 4.5)
SIDE, N, L = 512, 256, 4
SEED0 = 600000


def streak_texture(m, across, along, rng):
    """Random texture with correlation lengths ``across`` along x and ``along`` along y (unit std)."""
    f = np.fft.fftfreq(m)
    FX, FY = np.meshgrid(f, f)
    psd = (1 + (2 * np.pi) ** 2 * ((across * FX) ** 2 + (along * FY) ** 2)) ** -1.5
    t = np.real(np.fft.ifft2(np.fft.fft2(rng.standard_normal((m, m))) * np.sqrt(psd)))
    return (t - t.mean()) / t.std(ddof=1)


def cases(quick=False):
    """Static images for every pair of lengths (the pair (along, along) is the isotropic control), then the moving ones.

    The seeds count on in this order, so the full run's images do not depend on --quick.
    """
    out = []
    for along in ALONG:
        for across in sorted({*ACROSS, along}):
            for i in range(SEEDS):
                out.append(dict(across=across, along=along, speed=0.0, angle=0.0, seed=SEED0 + len(out), quick=i < 3))
    for across, along in MOVING:
        for angle in ANGLES:
            for speed in SPEEDS:
                out.append(dict(across=across, along=along, speed=speed, angle=angle, seed=SEED0 + len(out),
                                quick=speed in (1.0, 4.5)))
    return [{k: v for k, v in c.items() if k != "quick"} for c in out if c["quick"] or not quick]


def evaluate(case):
    """CSV row with the case and the estimate; a top-level function for the worker pool."""
    rng = np.random.default_rng(case["seed"])
    if case["speed"] == 0:
        img = crop_square(streak_texture(SIDE, case["across"], case["along"], rng), N)
    else:
        a = np.deg2rad(case["angle"])
        v = (case["speed"] * np.cos(a), -case["speed"] * np.sin(a))
        img = simulate_pulsed_image(N, v, L, texture=streak_texture(2 * SIDE, case["across"], case["along"], rng),
                                    rng=rng)
    e = estimate_velocity(img)
    strict = estimate_velocity(img, reject_blur=True, clip_outliers=True)
    row = dict(case, est_angle=e.angle_deg, est_speed=e.speed, quality=e.quality, threshold=e.threshold,
               quality_perp=e.quality_perp, isotropic=int(e.isotropic), in_range=int(e.in_range),
               quality_blur=e.quality_blur, blur_length=e.blur_length, blur_like=int(e.blur_like),
               valid=int(e.valid), valid_strict=int(strict.valid))
    if case["speed"] > 0:
        row["angle_err"] = abs(angle_difference(e.angle_deg, case["angle"]))
        row["speed_err"] = abs(e.speed / case["speed"] - 1)
        row["correct"] = int(row["angle_err"] <= 2 and row["speed_err"] <= 0.05)
    return row


def run(pool, quick):
    rows = common.run_cases(cases(quick), pool, evaluate)
    common.write_csv(CSV[0], rows)
    return rows


def summary(rows):
    static = [r for r in rows if r["speed"] == 0]
    per_pair = len([r for r in static if r["across"] == ACROSS[0] and r["along"] == ALONG[0]])
    print(f"static textures of streaks, {N} px, K = 9, {per_pair} images per pair of correlation lengths")
    print("   across  along  ratio  valid  strict  median speed of the valid  largest quality")
    for along in ALONG:
        for across in sorted({*ACROSS, along}):
            sel = [r for r in static if r["across"] == across and r["along"] == along]
            val = [r for r in sel if r["valid"]]
            med = f"{np.median([r['est_speed'] for r in val]):.2f} px" if val else "-"
            print(f"   {across:6.0f}  {along:5.0f}  {along / across:5.0f}  {len(val):2d}/{len(sel)}  "
                  f"{int(sum(r['valid_strict'] for r in sel)):6d}  {med:>26s}  {max(r['quality'] for r in sel):.3f}")
    aniso = [r for r in static if r["along"] / r["across"] >= 16]
    val = [r for r in aniso if r["valid"]]
    speeds = [r["est_speed"] for r in val] or [float("nan")]
    print(f"   ratio 16 or more: {len(val)} of {len(aniso)} valid ({int(sum(r['valid_strict'] for r in aniso))} with "
          f"reject_blur and clip_outliers), speeds {min(speeds):.2f} to {max(speeds):.2f} px per step "
          f"(median {np.median(speeds):.2f}), blur_like {int(sum(r['blur_like'] for r in val))}, "
          f"isotropic {int(sum(r['isotropic'] for r in aniso))} of {len(aniso)}")
    ctrl = [r for r in static if r["along"] == r["across"]]
    print(f"   isotropic controls: {int(sum(r['valid'] for r in ctrl))} of {len(ctrl)} valid, "
          f"largest quality {max(r['quality'] for r in ctrl):.3f}")
    moving = [r for r in rows if r["speed"] > 0]
    print(f"moving textures of streaks, {len(moving)} images: correct {int(sum(r['correct'] for r in moving))}, "
          f"valid {int(sum(r['valid'] for r in moving))}, wrong but valid "
          f"{int(sum(r['valid'] and not r['correct'] for r in moving))}, largest direction error "
          f"{max(r['angle_err'] for r in moving):.2f} deg, largest speed error "
          f"{100 * max(r['speed_err'] for r in moving):.2f} %")


def replot():
    summary(common.read_csv(CSV[0]))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--quick", action="store_true",
                   help="small sample; writes to experiments/results_quick unless --out is given")
    p.add_argument("--out", type=Path, metavar="DIR",
                   help="folder for the CSV file (default: experiments/results; with --quick "
                        "experiments/results_quick)")
    p.add_argument("--jobs", type=int, default=os.cpu_count(),
                   help="worker processes (default: number of CPUs)")
    p.add_argument("--replot", action="store_true",
                   help="only print the summary from the CSV file in the output folder")
    a = p.parse_args(argv)
    if a.jobs < 1:
        p.error("--jobs must be at least 1")
    out = (a.out if a.out is not None else common.QUICK if a.quick else common.COMMITTED).resolve()
    if a.quick and out == common.COMMITTED:
        p.error("--quick must not write to experiments/results, which holds the full-run results; "
                "leave out --out or choose another folder")
    common.set_output(out)
    print(f"CSV files: {common.RESULTS}")
    if a.replot:
        missing = common.missing_csv(CSV)
        if missing:
            sys.exit(f"--replot: {', '.join(missing)} not found in {common.RESULTS}")
        replot()
        return
    with Pool(a.jobs, initializer=common.init_worker) as pool:
        rows = run(pool, a.quick)
    summary(rows)


if __name__ == "__main__":
    main()
