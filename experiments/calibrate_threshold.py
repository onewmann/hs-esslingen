"""Calibration of the validity threshold: how high the fit quality of images
without motion gets, checked on held-out images and on single camera frames,
the share of correct estimates the threshold lets through, and the half-speed
failure.

    python experiments/calibrate_threshold.py [--quick] [--out DIR] [--replot] [--jobs N]
                                              [--recordings PATH]

An estimate is valid if ``quality >= max(0.6, tanh(c/sqrt(m - 5)))``, the
speed is in range and the spectrum is not isotropic; ``m`` is the number of
independent profile points and ``c = 9`` (``NULL_K`` in fftvel.core). Each
image gives the statistic ``k = atanh(quality) * sqrt(m - 5)``, so the
quality alone passes ``tanh(c/sqrt(m - 5))`` exactly when ``k >= c``. Parts:

1. calibration: 1008 images without motion (sizes 32, 48, 64, 96, 128, 192
   and 256 px, L = 1, 2, 4, 10, white noise, a static random texture and
   static crops of data/surface_texture.png with a little noise, 12 seeds
   each); the largest k and its 99th percentile.
2. held_out: new seeds, and textures the calibration did not see (other
   correlation lengths, power-law spectra, crops of the static room scene
   vid1.mat); the count of valid estimates with a Clopper-Pearson bound.
3. real: single frames of vid1 (static) and of vid2 to vid4 (moving surface,
   but one exposure per frame, so nothing repeats), crops of 256 and 512 px.
4. power: moving images of the camera texture with random velocity and L
   at noise 0, 0.25, 0.5 and 1; the share of correct estimates that pass
   for c = 8, 8.6, 9 and 10, and a table of the quality rule against GEBA.
5. half speed: for valid moving estimates, the fit score at half and at
   twice the estimated speed, correct against wrong estimates. Besides the
   images of part 4 (equal steps) this uses images with uneven steps
   (intervals alternating between 0.8 and 1.2 steps, or one flash missing),
   which the README names as the main source of twice the true speed.
6. blur: single exposures of the moving camera texture, blurred over 1 to
   12 px and without flashes, as a control for the single frames of part 3.
7. raw: the static crops of data/surface_texture.png of parts 1 and 2 again,
   at the same positions but without the added noise. The texture is itself
   a camera frame, so this is what a static surface gives; it shows how much
   the result of parts 1 and 2 depends on the added noise.

Writes threshold_calibration.csv, threshold_held_out.csv, threshold_real.csv,
threshold_power.csv, threshold_uneven.csv, threshold_blur.csv, threshold_raw.csv
(one row per image), threshold_roc.csv (accepted counts per rule and threshold)
and threshold.png. Rows of camera crops (vid1 in part 2, part 3) also hold
``saturated``, the share of pixels at 255.
A full run takes about 20 min on 4 cores. --quick runs a small sample and
writes to experiments/results_quick; --replot redraws the figure and prints
the numbers from the CSV files in the output folder. ``--recordings`` is the
folder with vid1.mat to vid4.mat (default: data/recordings). Every image is
fixed by the seed in its row (the crops of part 3 by frame and position, in
``where``), so a run gives the same CSV files on the same machine.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import common  # noqa: E402
from common import AQUA, BLUE, INK, INK2, MUTED, ORANGE, SURFACE  # noqa: E402
from fftvel import angle_difference, estimate_velocity, simulate_pulsed_image  # noqa: E402
from fftvel.core import DEFAULT_MIN_QUALITY, INDEPENDENT_PAD, NULL_K, NULL_OFFSET, SIM_MARGIN  # noqa: E402

CSV = ["threshold_calibration.csv", "threshold_held_out.csv", "threshold_real.csv", "threshold_power.csv",
       "threshold_uneven.csv", "threshold_blur.csv", "threshold_raw.csv"]
ROC_CSV = "threshold_roc.csv"
VIDEOS = ("vid1", "vid2", "vid3", "vid4")
FILES = [f"{v}.mat" for v in VIDEOS]

SIZES = (32, 48, 64, 96, 128, 192, 256)
LS = (1, 2, 4, 10)
NULL_KINDS = ("white", "random", "surface")
NEW_KINDS = ("corr_len", "power_law", "vid1")      # held-out textures the calibration did not see
CORR_LENS = (0.5, 4.0, 8.0)                        # px, random texture of make_texture (calibration: 1.5)
BETAS = (1.0, 2.0, 3.0)                            # power spectrum ~ |f|**-beta
NOISE_LEVELS = (0.0, 0.25, 0.5, 1.0)
CONSTANTS = (8.0, 8.6, 9.0, 10.0)                  # the 8.6 of the docstring, the constant 9 and two others
REAL_GRID = {256: (4, 7), 512: (2, 3)}             # roi -> rows x columns of non-overlapping crops
UNEVEN_KINDS = ("alternating", "dropped")
UNEVEN_SIZES = (64, 128, 256)
UNEVEN_NOISE = (0.0, 0.5)
BLUR_LENGTHS = (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0)    # px of motion during the single exposure
BLUR_SIZES = (128, 256)
BLUR_STEP = 0.05            # px between the positions that make up the continuous exposure
ALTERNATING = 0.2           # intervals alternate between 1 - 0.2 and 1 + 0.2 steps
STATIC_NOISE = 0.1          # white noise added to static crops and power-law textures, times their std
TOL_ANGLE, TOL_SPEED = 2.0, 0.05                   # correct: as in common.evaluate_case
DOUBLE_TOL = 0.1            # |log2(est/true) - 1| below this: twice the speed (half: + 1)
CONFIDENCE = 0.95           # of the one-sided Clopper-Pearson upper bounds
TEXTURE_SIDE = 1024         # side of data/surface_texture.png
SATURATED = 255             # value of a saturated pixel of the 8-bit recordings
SAT_SHARE = 0.1             # camera crops with at least this share of saturated pixels are counted apart

# seeds: one block per part; the seed of every image is stored in its row
SEED_CALIBRATION, SEED_HELD_OUT, SEED_POWER, SEED_UNEVEN, SEED_BLUR = 100000, 200000, 300000, 400000, 500000

# constants of the power curve in the figure, and thresholds of threshold_roc.csv
CURVE = np.round(np.arange(5.0, 12.0001, 0.1), 1)
K_THRESHOLDS = [float(t) for t in sorted({*np.round(np.arange(5.0, 12.0001, 0.5), 1), 8.6})]
GEBA_THRESHOLDS = [float(t) for t in np.round(np.arange(1.10, 2.0001, 0.05), 2)] + \
    [2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 5.0]
HALF_CUTS = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]        # reject a valid estimate if score_half/quality >= this
K_MAX_SHOWN = 14.0          # the figure clips k above this


# Cases

def calibration_cases(quick=False):
    """Images without motion of part 1: sizes x L x kinds x seeds."""
    sizes, ls, seeds = ((32, 64, 128), (1, 4), 2) if quick else (SIZES, LS, 12)
    cases = []
    for n in sizes:
        for L in ls:
            for kind in NULL_KINDS:
                for _ in range(seeds):
                    cases.append(dict(part="calibration", kind=kind, n=n, L=L,
                                      seed=SEED_CALIBRATION + len(cases)))
    return cases


def held_out_cases(recordings, quick=False):
    """Images without motion of part 2: new seeds for the calibration kinds, and new textures."""
    sizes, ls, seeds = ((32, 64, 128), (1, 4), 1) if quick else (SIZES, LS, 9)
    cases = []
    for n in sizes:
        for L in ls:
            for kind in NULL_KINDS + NEW_KINDS:
                for s in range(seeds):
                    c = dict(part="held_out", kind=kind, n=n, L=L, seed=SEED_HELD_OUT + len(cases))
                    if kind == "corr_len":
                        c["param"] = CORR_LENS[s % len(CORR_LENS)]
                    elif kind == "power_law":
                        c["param"] = BETAS[s % len(BETAS)]
                    elif kind == "vid1":
                        c["file"] = str(Path(recordings) / "vid1.mat")
                    cases.append(c)
    return cases


def real_cases(recordings, quick=False):
    """Single camera frames of part 3: every frame, a grid of crops per roi, L taken in turn."""
    frames = (0, 5) if quick else range(10)
    cases = []
    for video in VIDEOS:
        for f in frames:
            for roi, (nr, nc) in REAL_GRID.items():
                if quick:
                    nr, nc = (1, 2) if roi == 256 else (1, 1)
                for i in range(nr * nc):
                    cases.append(dict(part="real", kind=video, file=str(Path(recordings) / f"{video}.mat"),
                                      frame=f, n=roi, grid=(nr, nc), cell=divmod(i, nc), L=LS[(f + i) % len(LS)]))
    return cases


def _velocity(rng, n, K, span):
    """Random velocity: speed log-uniform from 4/K to n/4 (less if the texture is too small), any direction."""
    vmax = min(n / 4, (TEXTURE_SIDE - n - SIM_MARGIN - 4) / span)
    speed = float(np.exp(rng.uniform(np.log(4 / K), np.log(vmax))))
    angle = float(rng.uniform(0, 180))
    return speed, angle


def power_cases(quick=False):
    """Moving images of part 4: camera texture, random velocity and L, four noise levels."""
    sizes, per = ((32, 64, 128), 4) if quick else (SIZES, 48)
    cases = []
    for n in sizes:
        for noise in NOISE_LEVELS:
            for _ in range(per):
                seed = SEED_POWER + len(cases)
                rng = np.random.default_rng(seed)
                L = int(rng.integers(1, 11))
                speed, angle = _velocity(rng, n, 2 * L + 1, 2 * L)
                cases.append(dict(part="power", kind="equal", n=n, L=L, noise=noise, speed=speed,
                                  angle=angle, seed=seed))
    return cases


def uneven_cases(quick=False):
    """Moving images with uneven steps for part 5: alternating intervals, or one flash missing."""
    sizes, per = ((64,), 3) if quick else (UNEVEN_SIZES, 40)
    cases = []
    for kind in UNEVEN_KINDS:
        for n in sizes:
            for noise in UNEVEN_NOISE:
                for _ in range(per):
                    seed = SEED_UNEVEN + len(cases)
                    rng = np.random.default_rng(seed)
                    L = int(rng.integers(2, 11))
                    K = 2 * L + 1
                    speed, angle = _velocity(rng, n, K, K)
                    # alternating: K flashes, intervals 0.8 and 1.2 steps; dropped: K + 1 equal steps with
                    # one inner flash missing (param: its index), so K flashes with one doubled interval
                    param = ALTERNATING if kind == "alternating" else int(rng.integers(1, K))
                    cases.append(dict(part="uneven", kind=kind, n=n, L=L, noise=noise, speed=speed,
                                      angle=angle, param=param, seed=seed))
    return cases


def blur_cases(quick=False):
    """Single exposures for part 6: the texture moves by ``param`` px during one exposure, no flashes."""
    sizes, lengths, ls, seeds = ((128,), (2.0, 6.0), (1, 4), 1) if quick else (BLUR_SIZES, BLUR_LENGTHS, LS, 3)
    cases = []
    for n in sizes:
        for d in lengths:
            for L in ls:
                for _ in range(seeds):
                    seed = SEED_BLUR + len(cases)
                    angle = float(np.random.default_rng(seed).uniform(0, 180))
                    cases.append(dict(part="blur", kind="blur", n=n, L=L, param=d, noise=STATIC_NOISE,
                                      angle=angle, seed=seed))
    return cases


def raw_cases(recordings, quick=False):
    """Part 7: the camera-texture crops of parts 1 and 2 with their seeds, so at the same positions,
    but without the added noise."""
    return [dict(c, part="raw", kind="surface_raw")
            for c in calibration_cases(quick) + held_out_cases(recordings, quick) if c["kind"] == "surface"]


# Images

_recordings: dict = {}


def recording(path):
    """Frames of a recording (rows x cols x frames, uint8), loaded once per process."""
    if path not in _recordings:
        from scipy.io import loadmat
        _recordings[path] = loadmat(path)["recording"][:, :, 0, :]
    return _recordings[path]


def texture():
    """data/surface_texture.png as float; common.init_worker loads it in the workers."""
    if common._texture is None:
        common.init_worker()
    return common._texture


def power_law_texture(n, beta, rng):
    """Static texture with power spectrum ~ |f|**-beta, the central n x n of a larger periodic one."""
    m = int(2 ** np.ceil(np.log2(n + 2 * SIM_MARGIN)))
    f = np.fft.fftfreq(m)
    FX, FY = np.meshgrid(f, f)
    r = np.hypot(FX, FY)
    r[0, 0] = 1.0 / m                       # the mean is removed below
    t = np.real(np.fft.ifft2(np.fft.fft2(rng.standard_normal((m, m))) * r ** (-beta / 2)))
    r0 = (m - n) // 2
    t = t[r0:r0 + n, r0:r0 + n]
    return (t - t.mean()) / t.std(ddof=1)


def with_noise(x, rng, level=STATIC_NOISE):
    """``x`` plus white noise with ``level`` times the standard deviation of ``x``."""
    return x + level * x.std(ddof=1) * rng.standard_normal(x.shape)


def null_image(case):
    """Image without motion of a calibration or held-out case; returns (image, where it comes from)."""
    n, kind = case["n"], case["kind"]
    rng = np.random.default_rng(case["seed"])
    if kind == "white":
        return rng.standard_normal((n, n)), ""
    if kind == "random":
        return simulate_pulsed_image(n, (0.0, 0.0), case["L"], rng=rng), ""
    if kind in ("surface", "surface_raw"):     # the same crop for the same seed; raw: no added noise
        r, c = rng.integers(0, TEXTURE_SIDE - n + 1, size=2)
        crop = texture()[r:r + n, c:c + n]
        return (with_noise(crop, rng) if kind == "surface" else crop), f"r{r} c{c}"
    if kind == "corr_len":
        return simulate_pulsed_image(n, (0.0, 0.0), case["L"], corr_len=case["param"], rng=rng), ""
    if kind == "power_law":
        return with_noise(power_law_texture(n, case["param"], rng), rng), ""
    if kind == "vid1":
        rec = recording(case["file"])
        H, W, F = rec.shape
        f = int(rng.integers(0, F))
        r, c = int(rng.integers(0, H - n + 1)), int(rng.integers(0, W - n + 1))
        return rec[r:r + n, c:c + n, f].astype(float), f"frame {f + 1} r{r} c{c}"
    raise ValueError(f"unknown kind {kind}")


def real_image(case):
    """Crop of one camera frame at its cell of a centred grid of non-overlapping crops."""
    rec = recording(case["file"])
    H, W, _ = rec.shape
    n, (nr, nc), (i, j) = case["n"], case["grid"], case["cell"]
    r = (H - nr * n) // 2 + i * n
    c = (W - nc * n) // 2 + j * n
    return rec[r:r + n, c:c + n, case["frame"]].astype(float), f"frame {case['frame'] + 1} r{r} c{c}"


def _texture_crop(rng, n, reach):
    """Random square crop of the camera texture, large enough for shifts up to ``reach`` px each way."""
    side = int(np.ceil(n + 2 * reach + SIM_MARGIN)) + 2
    side = min(TEXTURE_SIDE, side + side % 2)
    r, c = rng.integers(0, TEXTURE_SIDE - side + 1, size=2)
    return texture()[r:r + side, c:c + side], f"r{r} c{c}"


def pulsed_at(tex, n, v, positions):
    """Mean of exposures at ``positions`` (in steps) of ``tex`` moving by ``v`` per step, central n x n.

    Exact Fourier shifts as in simulate_pulsed_image, but at any positions,
    so the steps between the exposures need not be equal.
    """
    m = tex.shape[0]
    f = np.fft.fftfreq(m)
    FX, FY = np.meshgrid(f, f)
    phase = -2j * np.pi * (FX * v[0] + FY * v[1])
    D = np.zeros((m, m), dtype=complex)
    for p in positions:
        D += np.exp(phase * p)
    big = np.real(np.fft.ifft2(np.fft.fft2(tex - tex.mean()) * D / len(positions)))
    r0 = (m - n) // 2
    return big[r0:r0 + n, r0:r0 + n]


def flash_positions(case):
    """Times of the flashes of an uneven case in steps, centred on 0."""
    K = 2 * case["L"] + 1
    if case["kind"] == "alternating":
        a = case["param"]
        p = np.concatenate([[0.0], np.cumsum([1 - a if i % 2 == 0 else 1 + a for i in range(K - 1)])])
    else:
        p = np.delete(np.arange(K + 1, dtype=float), int(case["param"]))
    return p - (p.max() + p.min()) / 2


def moving_image(case):
    """Image of a power, uneven or blur case: a random crop of the camera texture moves by v per step."""
    n, L = case["n"], case["L"]
    rng = np.random.default_rng((case["seed"], 1))      # not the stream that drew the case
    a = np.deg2rad(case["angle"])
    speed = case.get("speed", 1.0)
    v = (speed * np.cos(a), -speed * np.sin(a))
    vm = max(abs(v[0]), abs(v[1]))
    if case["part"] == "power":
        tex, where = _texture_crop(rng, n, L * vm)
        return simulate_pulsed_image(n, v, L, texture=tex, noise=case["noise"], rng=rng), where
    if case["part"] == "blur":            # v is the unit vector; positions cover the blur length evenly
        d = case["param"]
        pos = np.linspace(-d / 2, d / 2, int(np.ceil(d / BLUR_STEP)) + 1)
        tex, where = _texture_crop(rng, n, d / 2)
        return with_noise(pulsed_at(tex, n, v, pos), rng, case["noise"]), where
    pos = flash_positions(case)
    tex, where = _texture_crop(rng, n, np.abs(pos).max() * vm)
    bp = pulsed_at(tex, n, v, pos)
    return (with_noise(bp, rng, case["noise"]) if case["noise"] > 0 else bp), where


# Estimates

def independent_points(est):
    """m of the threshold: independent points of the speed profile, as estimate_velocity counts them."""
    npts = len(est.diag["profile_j"])
    return min(npts, (INDEPENDENT_PAD * est.n * npts) // est.npad)


def null_k(quality, m):
    """Statistic k = atanh(quality) * sqrt(m - 5); the quality passes tanh(c/sqrt(m - 5)) iff k >= c."""
    q = float(quality)
    if abs(q) >= 1:
        return math.copysign(math.inf, q)
    return math.atanh(q) * math.sqrt(m - NULL_OFFSET)


def quality_threshold(m, c, min_quality=DEFAULT_MIN_QUALITY):
    """Quality needed with the constant c, as estimate_velocity computes it for c = NULL_K."""
    return max(min_quality, math.inf if m <= NULL_OFFSET else math.tanh(c / math.sqrt(m - NULL_OFFSET)))


def passes(row, c):
    """Valid with the constant c in place of NULL_K (quality, range and isotropy)."""
    return bool(row["quality"] >= quality_threshold(int(row["m_eff"]), c)
                and row["in_range"] and not row["isotropic"])


def score_at(est, speed):
    """Fit score at ``speed`` on the speed grid (linear in log speed); NaN outside the grid."""
    v, s = est.diag["speed_grid"], est.diag["speed_score"]
    if not v[0] <= speed <= v[-1]:
        return math.nan
    return float(np.interp(np.log(speed), np.log(v), s))


def evaluate(case, img, where):
    """CSV row of one image: the case, the estimate, m, k and the scores at half and twice the speed."""
    est = estimate_velocity(img, L=case["L"], roi=case["n"], diagnostics=True)
    m = independent_points(est)
    row = dict(part=case["part"], kind=case["kind"], n=est.n, L=case["L"], seed=case.get("seed", ""),
               param=case.get("param", ""), where=where, est_angle=est.angle_deg, est_speed=est.speed,
               quality=est.quality, quality_perp=est.quality_perp, geba=est.geba, in_range=int(est.in_range),
               isotropic=int(est.isotropic), threshold=est.threshold, valid=int(est.valid), m_eff=m,
               k=null_k(est.quality, m), score_half=score_at(est, est.speed / 2),
               score_double=score_at(est, 2 * est.speed))
    # the threshold and the validity must be those estimate_velocity applied
    if not math.isclose(est.threshold, quality_threshold(m, NULL_K), rel_tol=1e-12) \
            or passes(row, NULL_K) != est.valid:
        raise RuntimeError(f"m = {m} does not reproduce the threshold {est.threshold} of estimate_velocity; "
                           f"make independent_points() match fftvel.core")
    if case["part"] in ("power", "uneven"):
        row.update(noise=case["noise"], speed=case["speed"], angle=case["angle"])
        row["angle_err"] = abs(angle_difference(est.angle_deg, case["angle"]))
        row["speed_err"] = abs(est.speed / case["speed"] - 1)
        row["correct"] = int(row["angle_err"] <= TOL_ANGLE and row["speed_err"] <= TOL_SPEED)
        row["failure"] = failure(row, est.speed / case["speed"])
    if case["part"] == "blur":            # one exposure blurred over d px; fooled: speed = d / K
        row.update(noise=case["noise"], angle=case["angle"], blur_px=case["param"],
                   k_speed=est.speed * est.exposures, angle_err=abs(angle_difference(est.angle_deg, case["angle"])))
    if case["part"] == "real":
        row["frame"] = case["frame"] + 1
    return row


def failure(row, ratio):
    """Kind of a wrong estimate: 'double' or 'half' the speed, 'angle' only, 'other'; '' if correct."""
    if row["correct"]:
        return ""
    if abs(math.log2(ratio) - 1) < DOUBLE_TOL:
        return "double"
    if abs(math.log2(ratio) + 1) < DOUBLE_TOL:
        return "half"
    return "angle" if row["speed_err"] <= TOL_SPEED else "other"


def run_case(case):
    """One case -> CSV row (top-level, for the worker pool)."""
    if case["part"] in ("power", "uneven", "blur"):
        img, where = moving_image(case)
    elif case["part"] == "real":
        img, where = real_image(case)
    else:
        img, where = null_image(case)
    row = evaluate(case, img, where)
    if "file" in case:                    # camera crop (8-bit values as float)
        row["saturated"] = float(np.mean(img >= SATURATED))
    return row


# Statistics

def upper_bound(k, n, confidence=CONFIDENCE):
    """One-sided Clopper-Pearson upper confidence bound of a rate with k events in n trials."""
    if n == 0:
        return math.nan
    if k >= n:
        return 1.0
    from scipy.stats import beta
    return float(beta.ppf(confidence, k + 1, n - k))


def auc(pos, neg):
    """Probability that a positive scores higher than a negative (ties count half); NaN if one is empty."""
    pos = np.asarray([p for p in pos if np.isfinite(p)], dtype=float)
    neg = np.asarray([x for x in neg if np.isfinite(x)], dtype=float)
    if pos.size == 0 or neg.size == 0:
        return math.nan
    d = pos[:, None] - neg[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / d.size)


def k_stats(rows):
    """Largest k, its row, the 99th percentile and median, and the counts at the constants."""
    k = np.array([r["k"] for r in rows])
    out = dict(n=len(rows), k_max=float(k.max()), k_p99=float(np.percentile(k, 99)),
               k_p50=float(np.median(k)), top=rows[int(np.argmax(k))], valid=int(sum(r["valid"] for r in rows)))
    for c in CONSTANTS:
        out[f"above_{c:g}"] = int((k >= c).sum())
        out[f"valid_{c:g}"] = int(sum(passes(r, c) for r in rows))
    return out


def gated(r):
    """In range and not isotropic: the other two conditions of valid."""
    return bool(r["in_range"]) and not bool(r["isotropic"])


def half_ratio(r):
    return r["score_half"] / r["quality"]


def roc_rows(power, null, real, uneven):
    """Accepted counts of correct and wrong moving estimates, null images and camera frames per rule.

    Rules (all need in range and not isotropic): ``k_floor``, the estimator's
    rule with the constant c, quality >= max(0.6, tanh(c/sqrt(m-5)));
    ``k``, the same without the floor 0.6; ``geba``, GEBA >= g alone;
    ``valid_and_geba``, the current rule (c = 9) and GEBA >= g. These use
    the moving images of part 4 (sample 'power'). ``valid_and_half_below``,
    the current rule and score_half/quality < t, uses part 4 and the uneven
    steps of part 5 (sample 'power+uneven').
    """
    rules = [("k_floor", K_THRESHOLDS, lambda r, t: passes(r, t)),
             ("k", K_THRESHOLDS, lambda r, t: gated(r) and r["k"] >= t),
             ("geba", GEBA_THRESHOLDS, lambda r, t: gated(r) and r["geba"] >= t),
             ("valid_and_geba", GEBA_THRESHOLDS, lambda r, t: passes(r, NULL_K) and r["geba"] >= t),
             ("valid_and_half_below", HALF_CUTS,
              lambda r, t: passes(r, NULL_K) and not half_ratio(r) >= t)]     # NaN (off the grid) passes
    out = []
    for rule, thresholds, accept in rules:
        moving = power + uneven if rule == "valid_and_half_below" else power
        groups = dict(correct=[r for r in moving if r["correct"]], wrong=[r for r in moving if not r["correct"]],
                      null=null, real=real)
        for t in thresholds:
            row = dict(rule=rule, sample="power+uneven" if rule == "valid_and_half_below" else "power",
                       threshold=float(t))
            for name, sel in groups.items():
                a = int(sum(accept(r, t) for r in sel))
                row[f"accepted_{name}"] = a
                row[f"n_{name}"] = len(sel)
                row[f"share_{name}"] = a / len(sel) if sel else math.nan
            out.append(row)
    return out


def half_groups(power, uneven):
    """Valid moving estimates by sample and outcome, for the half-speed check."""
    out = {}
    for sample, rows in (("equal steps", power), ("alternating", [r for r in uneven if r["kind"] == "alternating"]),
                         ("dropped", [r for r in uneven if r["kind"] == "dropped"]), ("all", power + uneven)):
        val = [r for r in rows if r["valid"]]
        out[sample] = {"correct": [r for r in val if r["correct"]],
                       "double": [r for r in val if r["failure"] == "double"],
                       "other wrong": [r for r in val if not r["correct"] and r["failure"] != "double"]}
    return out


def power_groups(rows):
    """(label, rows) of the size groups in the figure."""
    groups = [("32 and 48 px", (32, 48)), ("64 and 96 px", (64, 96)), ("128 and 192 px", (128, 192)),
              ("256 px", (256,))]
    return [(label, [r for r in rows if int(r["n"]) in sizes]) for label, sizes in groups]


def analyse(cal, held, real, power, uneven, blur, raw):
    """All numbers of the study from the rows of the seven per-image CSV files."""
    null = cal + held
    sim = [r for r in null if r["kind"] != "vid1"]          # simulated images without motion
    res = dict(cal=cal, held=held, real=real, power=power, uneven=uneven, blur=blur, raw=raw, null=null, sim=sim,
               cal_stats=k_stats(cal), held_stats=k_stats(held), real_stats=k_stats(real), sim_stats=k_stats(sim),
               raw_stats=k_stats(raw))
    res["roc"] = roc_rows(power, null, real, uneven)
    res["half"] = half_groups(power, uneven)
    correct = [r for r in power if r["correct"]]
    curve = {label: [share_valid(sel, c) for c in CURVE] for label, sel in power_groups(correct)}
    curve["null"] = [share_valid(null, c) for c in CURVE]
    res["curve"] = curve
    return res


# Run

def run(pool, quick, recordings=None):
    """All seven parts; writes the CSV files and the figure and prints the numbers."""
    recordings = Path(recordings if recordings is not None else common.RECORDINGS)
    parts = {"threshold_calibration.csv": calibration_cases(quick),
             "threshold_held_out.csv": held_out_cases(recordings, quick),
             "threshold_real.csv": real_cases(recordings, quick),
             "threshold_power.csv": power_cases(quick),
             "threshold_uneven.csv": uneven_cases(quick),
             "threshold_blur.csv": blur_cases(quick),
             "threshold_raw.csv": raw_cases(recordings, quick)}
    rows = {}
    for name, cases in parts.items():
        print(f"{name}: {len(cases)} images", flush=True)
        rows[name] = common.run_cases(cases, pool, worker=run_case)
        common.write_csv(name, rows[name])
    res = analyse(*(rows[n] for n in CSV))
    common.write_csv(ROC_CSV, res["roc"])
    plot(res)
    summary(res)
    return res


# Figure

def plot(res):
    fig, (ax1, ax2) = common.figure(2, 9.6)
    sets = [(np.array([r["k"] for r in res[key]]), key) for key in ("cal", "held", "real")]
    lo = float(np.floor(min(k.min() for k, _ in sets))) - 1
    hi = K_MAX_SHOWN
    bins = np.arange(lo, hi + 0.25, 0.25)
    clip = [np.clip(k, lo, hi - 0.01) for k, _ in sets]
    nvalid = [int(sum(r["valid"] for r in res[key])) for key in ("cal", "held", "real")]
    ax1.hist(clip[0], bins=bins, color=BLUE, edgecolor=SURFACE, linewidth=0.4,
             label=f"calibration ({sets[0][0].size} images, {nvalid[0]} valid)")
    ax1.hist(clip[1], bins=bins, histtype="step", color=ORANGE, lw=1.6,
             label=f"held out ({sets[1][0].size}, {nvalid[1]} valid)")
    ax1.hist(clip[2], bins=bins, histtype="step", color=AQUA, lw=1.6,
             label=f"single camera frames ({sets[2][0].size}, {nvalid[2]} valid)")
    ax1.axvline(NULL_K, color=INK, ls="--", lw=1.2)
    ax1.text(NULL_K + 0.15, 0.95, f"constant {NULL_K:g}", color=INK, fontsize=9, va="top",
             transform=ax1.get_xaxis_transform())
    st = res["cal_stats"]
    ax1.axvline(st["k_max"], color=INK2, ls=":", lw=1.2)
    ax1.text(st["k_max"] - 0.15, 0.72, f"largest in\ncalibration\n{st['k_max']:.2f}", color=INK2, fontsize=8.5,
             ha="right", va="top", transform=ax1.get_xaxis_transform())
    over = [r for r in res["real"] + res["held"] if r["k"] >= hi]
    common.style(ax1, "k = atanh(quality) * sqrt(m - 5)", "images (no repeated exposures)")
    ax1.set_xlim(lo, hi)
    ax1.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="upper left")
    if over:
        ax1.text(0.99, 0.42, f"{len(over)} images with k > {hi:g}\nin the last bin,\n"
                 f"{int(sum(r['valid'] for r in over))} of them valid", color=INK2, fontsize=8.5, ha="right",
                 transform=ax1.transAxes)

    correct = [r for r in res["power"] if r["correct"]]
    for (label, sel), col in zip(power_groups(correct), (BLUE, ORANGE, AQUA, INK2)):
        ax2.plot(CURVE, res["curve"][label], color=col, lw=2, label=f"correct, {label} ({len(sel)})")
    ax2.plot(CURVE, res["curve"]["null"], color=MUTED, lw=2, ls="--",
             label=f"no motion, calibration and held out ({len(res['null'])})")
    ax2.axvline(NULL_K, color=INK, ls="--", lw=1.2)
    common.style(ax2, "constant c in tanh(c / sqrt(m - 5))", "marked valid (%)")
    ax2.set_ylim(0, 104)
    ax2.set_xlim(CURVE[0], CURVE[-1])
    ax2.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="center left")
    p9 = share_valid(correct, NULL_K)
    common.title(fig, "Validity threshold: images without motion and correct estimates",
                 f"{st['n']} calibration images: largest k {st['k_max']:.2f}, 99th percentile {st['k_p99']:.2f}; "
                 f"with c = {NULL_K:g}, {p9:.1f} % of {len(correct)} correct moving estimates are valid")
    common.save(fig, "threshold.png")


# Numbers

def _bound(k, n):
    return f"{k} of {n} valid, {100 * CONFIDENCE:.0f} % upper bound {100 * upper_bound(k, n):.2f} %"


def _nanstat(f, x):
    x = np.asarray(x, dtype=float)
    return f(x[np.isfinite(x)]) if np.isfinite(x).any() else math.nan


def share_valid(rows, c):
    """Percentage of ``rows`` valid with the constant c; NaN for no rows."""
    return 100 * float(np.mean([passes(r, c) for r in rows])) if rows else math.nan


def _num(x):
    """Float of a cell (from a row or a CSV file); NaN if empty."""
    return math.nan if x == "" or x is None else float(x)


def _spread(x, unit=""):
    """'median (quartiles a to b, range c to d)' of the values ``x``."""
    q = np.percentile(x, [0, 25, 50, 75, 100])
    return f"{q[2]:.2f}{unit} (quartiles {q[1]:.2f} to {q[3]:.2f}, range {q[0]:.2f} to {q[4]:.2f})"


def _saturation_split(rows):
    """'valid of crops with >= SAT_SHARE saturated pixels; valid of the others'."""
    hi = [r for r in rows if _num(r.get("saturated", "")) >= SAT_SHARE]
    lo = [r for r in rows if not _num(r.get("saturated", "")) >= SAT_SHARE]
    return (f"crops with at least {100 * SAT_SHARE:.0f} % saturated pixels: {int(sum(r['valid'] for r in hi))} of "
            f"{len(hi)} valid; others: {int(sum(r['valid'] for r in lo))} of {len(lo)} valid; fully saturated "
            f"(constant) crops: {sum(_num(r.get('saturated', '')) == 1 for r in rows)}")


def summary(res):
    """Print every number the documentation quotes."""
    st = res["cal_stats"]
    top = st["top"]
    print(f"\n1. Calibration, {st['n']} images without motion:")
    print(f"   largest k = {st['k_max']:.3f} ({top['kind']}, n = {int(top['n'])}, L = {int(top['L'])}, "
          f"seed {int(top['seed'])}, quality {top['quality']:.4f}, m = {int(top['m_eff'])}); "
          f"99th percentile {st['k_p99']:.3f}; median {st['k_p50']:.3f}")
    for kind in NULL_KINDS:
        sel = [r["k"] for r in res["cal"] if r["kind"] == kind]
        print(f"   {kind:9s}: largest k {max(sel):.3f}, 99th percentile {np.percentile(sel, 99):.3f}")
    for L in sorted({int(r["L"]) for r in res["cal"]}):
        sel = [r["k"] for r in res["cal"] if int(r["L"]) == L]
        print(f"   L = {L:2d}   : largest k {max(sel):.3f}, 99th percentile {np.percentile(sel, 99):.3f}")
    print("   n    m  threshold  largest k  largest quality")
    for n in sorted({int(r["n"]) for r in res["cal"]}):
        sel = [r for r in res["cal"] if int(r["n"]) == n]
        print(f"   {n:3d} {int(sel[0]['m_eff']):4d}  {sel[0]['threshold']:9.3f}  {max(r['k'] for r in sel):9.3f}  "
              f"{max(r['quality'] for r in sel):15.3f}")
    print("   " + ", ".join(f"k >= {c:g}: {st[f'above_{c:g}']}" for c in CONSTANTS)
          + f"; {_bound(st['valid'], st['n'])}")

    hs = res["held_stats"]
    print(f"\n2. Held out, {hs['n']} images without motion:")
    print(f"   {_bound(hs['valid'], hs['n'])}; largest k {hs['k_max']:.3f} ({hs['top']['kind']}, "
          f"n = {int(hs['top']['n'])}, L = {int(hs['top']['L'])}, {hs['top']['where']}); 99th percentile "
          f"{hs['k_p99']:.3f}; " + ", ".join(f"k >= {c:g}: {hs[f'above_{c:g}']}" for c in CONSTANTS))
    for kind in NULL_KINDS + NEW_KINDS:
        sel = [r for r in res["held"] if r["kind"] == kind]
        if sel:
            print(f"   {kind:9s}: {_bound(int(sum(r['valid'] for r in sel)), len(sel))}, "
                  f"largest k {max(r['k'] for r in sel):.3f}, k >= {NULL_K:g}: "
                  f"{int(sum(r['k'] >= NULL_K for r in sel))}")
    sim = [r for r in res["held"] if r["kind"] != "vid1"]
    print(f"   simulated textures only: {_bound(int(sum(r['valid'] for r in sim)), len(sim))}")
    for r in res["held"]:
        if r["valid"]:
            print(f"   valid: {r['kind']} n = {int(r['n'])} L = {int(r['L'])} {r['where']} quality "
                  f"{r['quality']:.3f} (threshold {r['threshold']:.3f}), k {r['k']:.2f}, GEBA {r['geba']:.2f}, "
                  f"speed {r['est_speed']:.2f} at {r['est_angle']:.1f} deg"
                  + (f", {100 * _num(r['saturated']):.1f} % saturated"
                     if math.isfinite(_num(r.get("saturated", ""))) else ""))
    print(f"   vid1 {_saturation_split([r for r in res['held'] if r['kind'] == 'vid1'])}")
    null = res["null"]
    print(f"   calibration and held out together: {_bound(int(sum(r['valid'] for r in null)), len(null))}")
    ss = res["sim_stats"]
    print(f"   simulated images without motion (calibration and held out without vid1), {ss['n']} images: "
          f"{_bound(ss['valid'], ss['n'])}; largest k {ss['k_max']:.3f} ({ss['top']['part']}, {ss['top']['kind']}, "
          f"n = {int(ss['top']['n'])}, L = {int(ss['top']['L'])}, seed {int(ss['top']['seed'])}); "
          f"99th percentile {ss['k_p99']:.3f}; 99.9th percentile "
          f"{np.percentile([r['k'] for r in res['sim']], 99.9):.3f}; "
          + ", ".join(f"k >= {c:g}: {ss[f'above_{c:g}']}" for c in CONSTANTS))
    print("   by size: " + "; ".join(
        f"{n} px {len(sel)} images, k >= 8.6: {sum(r['k'] >= 8.6 for r in sel)}, "
        f"k >= {NULL_K:g}: {sum(r['k'] >= NULL_K for r in sel)}"
        for n in SIZES for sel in [[r for r in res["sim"] if int(r["n"]) == n]] if sel))

    rs = res["real_stats"]
    print(f"\n3. Single camera frames, {rs['n']} crops:")
    print(f"   {_bound(rs['valid'], rs['n'])}; largest k {rs['k_max']:.3f} ({rs['top']['kind']}, "
          f"{rs['top']['where']}, n = {int(rs['top']['n'])}, L = {int(rs['top']['L'])}, "
          f"quality {rs['top']['quality']:.3f}, in range {int(rs['top']['in_range'])})")
    print("   video  roi  valid  CP bound  k >= 9  quality >= threshold  of these in range  isotropic  largest k")
    for video in VIDEOS:
        for roi in sorted({int(r["n"]) for r in res["real"]}):
            sel = [r for r in res["real"] if r["kind"] == video and int(r["n"]) == roi]
            if sel:
                qp = [r for r in sel if r["quality"] >= r["threshold"]]
                print(f"   {video}  {roi:4d}  {int(sum(r['valid'] for r in sel)):2d}/{len(sel):3d}  "
                      f"{100 * upper_bound(int(sum(r['valid'] for r in sel)), len(sel)):6.2f} %  "
                      f"{int(sum(r['k'] >= NULL_K for r in sel)):6d}  {len(qp):20d}  "
                      f"{int(sum(r['in_range'] for r in qp)):17d}  {int(sum(r['isotropic'] for r in sel)):9d}  "
                      f"{max(r['k'] for r in sel):9.2f}")
    for video in VIDEOS:
        sel = [r for r in res["real"] if r["kind"] == video]
        print(f"   {video} all crops: {_bound(int(sum(r['valid'] for r in sel)), len(sel))}; {_saturation_split(sel)}")
    print("   valid crops (K = 2L + 1; K * speed is the length a blur of one exposure would have):")
    for video in VIDEOS:
        val = [r for r in res["real"] if r["kind"] == video and r["valid"]]
        if not val:
            continue
        print(f"   {video}: speed {_spread([r['est_speed'] for r in val], ' px')}; K * speed "
              f"{_spread([r['est_speed'] * (2 * int(r['L']) + 1) for r in val], ' px')}; direction "
              f"{_spread([r['est_angle'] for r in val], ' deg')}; L {sorted({int(r['L']) for r in val})}; "
              f"frames {sorted(int(r['frame']) for r in val)}")
        if video == "vid1":                # the static scene: every valid crop
            for r in val:
                print(f"      n = {int(r['n'])} L = {int(r['L']):2d} {r['where']}: speed {r['est_speed']:.2f} px, "
                      f"K * speed {r['est_speed'] * (2 * int(r['L']) + 1):.1f} px, {r['est_angle']:.1f} deg, "
                      f"k {r['k']:.2f}, GEBA {r['geba']:.2f}, {100 * _num(r['saturated']):.1f} % saturated")
    qp = [r for r in res["real"] if r["quality"] >= r["threshold"]]
    print(f"   quality at or above the threshold: {len(qp)} crops; their speeds "
          f"{min((r['est_speed'] for r in qp), default=math.nan):.3f} to "
          f"{max((r['est_speed'] for r in qp), default=math.nan):.3f} px, in range: "
          f"{int(sum(r['in_range'] for r in qp))}")

    power = res["power"]
    correct = [r for r in power if r["correct"]]
    wrong = [r for r in power if not r["correct"]]
    print(f"\n4. Power, {len(power)} moving images (equal steps): {len(correct)} correct, {len(wrong)} wrong")
    print("   share of correct estimates that are valid (%) with the constant c; wrong estimates, null images "
          "and camera frames valid (count):")
    print("   c       " + "".join(f"{f'n={n}':>8s}" for n in SIZES) + "     all   wrong    null    real")
    for c in CONSTANTS:
        cells = []
        for n in SIZES:
            sel = [r for r in correct if int(r["n"]) == n]
            cells.append(f"{share_valid(sel, c):8.1f}" if sel else f"{'':8s}")
        print(f"   {c:<6g}  " + "".join(cells) + f"{share_valid(correct, c):8.1f}"
              f"{sum(passes(r, c) for r in wrong):8d}{sum(passes(r, c) for r in null):8d}"
              f"{sum(passes(r, c) for r in res['real']):8d}")
    print("   correct estimates per size: " + ", ".join(
        f"n={n}: {sum(int(r['n']) == n for r in correct)}/{sum(int(r['n']) == n for r in power)}" for n in SIZES))
    for nl in NOISE_LEVELS:
        sel = [r for r in power if r["noise"] == nl]
        cc = [r for r in sel if r["correct"]]
        print(f"   noise {nl:4.2f}: {len(cc)} of {len(sel)} correct; valid among correct with c = "
              + ", ".join(f"{c:g}: {share_valid(cc, c):.1f} %" for c in CONSTANTS))
    kinds = {}
    for r in wrong:
        kinds[r["failure"]] = kinds.get(r["failure"], 0) + 1
    vw = [r for r in wrong if r["valid"]]
    listed = "; ".join(f"{r['failure']} n={int(r['n'])} L={int(r['L'])} noise {r['noise']:g} "
                       f"speed {r['speed']:.2f} -> {r['est_speed']:.2f}, angle error {r['angle_err']:.1f}"
                       for r in vw)
    print(f"   wrong estimates by kind: {dict(sorted(kinds.items()))}; valid: {len(vw)} of {len(power) - len(correct)}"
          f" ({listed})")
    nv = int(sum(r["valid"] for r in power))
    print(f"   valid estimates that are correct: {nv - len(vw)} of {nv} ({100 * (nv - len(vw)) / max(nv, 1):.1f} %)")

    print("\n   Quality rule against GEBA (every rule also needs in range and not isotropic;"
          " null = calibration + held out):")
    print(f"   {'rule':15s} {'threshold':>9s} {'correct %':>9s} {'wrong':>6s} {'null':>5s} {'real':>5s}")
    show = {"k_floor": (7.0, 8.0, 8.6, 9.0, 10.0, 11.0), "k": (7.0, 8.0, 9.0, 10.0),
            "geba": (1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0), "valid_and_geba": (1.2, 1.3, 1.4, 1.5, 1.6)}
    for r in res["roc"]:
        if r["threshold"] in show.get(r["rule"], ()):
            print(f"   {r['rule']:15s} {r['threshold']:9g} {100 * r['share_correct']:9.1f} "
                  f"{r['accepted_wrong']:6d} {r['accepted_null']:5d} {r['accepted_real']:5d}")
    pos = [r for r in correct if gated(r)]
    neg = [r for r in null if gated(r)]
    wg = [r for r in wrong if gated(r)]

    def aucs(rows, keys):
        return ", ".join(f"{key} {auc([r[key] for r in pos], [r[key] for r in rows]):.4f}" for key in keys)
    print(f"   AUC correct moving against null (all in range, not isotropic): {aucs(neg, ('k', 'quality', 'geba'))}")
    print(f"   AUC correct against wrong moving (all in range, not isotropic): {aucs(wg, ('k', 'quality', 'geba'))}")
    cur = next(r for r in res["roc"] if r["rule"] == "k_floor" and r["threshold"] == NULL_K)
    alone = [r for r in res["roc"] if r["rule"] == "geba" and r["share_correct"] >= cur["share_correct"]]
    if alone:
        g = alone[-1]
        print(f"   GEBA alone, highest threshold that keeps the share of correct estimates of the current rule "
              f"({100 * cur['share_correct']:.1f} %): GEBA >= {g['threshold']:g} gives "
              f"{100 * g['share_correct']:.1f} %, {g['accepted_wrong']} wrong, {g['accepted_null']} null, "
              f"{g['accepted_real']} real (current rule: {cur['accepted_wrong']} wrong, {cur['accepted_null']} null, "
              f"{cur['accepted_real']} real)")
    lowest = min((r["geba"] for r in power if r["valid"] and r["correct"]), default=math.nan)
    print(f"   lowest GEBA of a correct valid estimate: {lowest:.3f}; GEBA of the wrong valid ones: "
          + ", ".join(f"{r['geba']:.3f}" for r in vw))

    print("\n5. Half-speed failure: score at half (and twice) the estimated speed / quality, valid estimates")
    for sample, g in res["half"].items():
        parts = []
        for name in ("correct", "double", "other wrong"):
            sel = g[name]
            if not sel:
                parts.append(f"{name} 0")
                continue
            h = [half_ratio(r) for r in sel]
            d = [r["score_double"] / r["quality"] for r in sel]
            parts.append(f"{name} {len(sel)} (half: median {_nanstat(np.median, h):.3f}, "
                         f"range {_nanstat(np.min, h):.3f} to {_nanstat(np.max, h):.3f}; twice: median "
                         f"{_nanstat(np.median, d):.3f})")
        print(f"   {sample:11s}: " + "; ".join(parts))
        if g["correct"] and g["double"]:
            def twice(rows):
                return [r["score_double"] / r["quality"] for r in rows]
            print(f"   {'':11s}  AUC twice-speed against correct: score at half "
                  f"{auc([half_ratio(r) for r in g['double']], [half_ratio(r) for r in g['correct']]):.3f}, "
                  f"score at twice {auc(twice(g['double']), twice(g['correct'])):.3f}")
    print("   rejecting valid estimates with score_half/quality >= t (equal and uneven steps together):")
    print(f"   {'t':>5s} {'correct kept':>13s} {'twice kept':>11s} {'other wrong kept':>17s} {'null kept':>10s}")
    allg = res["half"]["all"]
    for r in res["roc"]:
        if r["rule"] == "valid_and_half_below":
            t = r["threshold"]
            dk = sum(not half_ratio(x) >= t for x in allg["double"])
            ok = sum(not half_ratio(x) >= t for x in allg["other wrong"])
            print(f"   {t:5.2f} {r['accepted_correct']:6d} of {len(allg['correct']):4d} "
                  f"{dk:4d} of {len(allg['double']):4d} {ok:8d} of {len(allg['other wrong']):5d} "
                  f"{r['accepted_null']:10d}")
    for kind in UNEVEN_KINDS:
        sel = [r for r in res["uneven"] if r["kind"] == kind]
        kinds = {}
        for r in sel:
            if not r["correct"]:
                kinds[r["failure"]] = kinds.get(r["failure"], 0) + 1
        vw = [r for r in sel if r["valid"] and not r["correct"]]
        vk = {}
        for r in vw:
            vk[r["failure"]] = vk.get(r["failure"], 0) + 1
        print(f"   {kind}: {len(sel)} images, {int(sum(r['correct'] for r in sel))} correct, "
              f"{int(sum(r['valid'] for r in sel))} valid, {len(vw)} valid but wrong; wrong by kind "
              f"{dict(sorted(kinds.items()))}; valid but wrong by kind {dict(sorted(vk.items()))}"
              + (f"; estimated / true speed of these {_spread([r['est_speed'] / r['speed'] for r in vw])}"
                 if vw else ""))

    blur = res["blur"]
    print(f"\n6. Single exposures blurred over d px, no flashes ({len(blur)} images): valid count, and for the valid "
          f"ones K * speed / d and the direction error")
    for d in sorted({r["blur_px"] for r in blur}):
        sel = [r for r in blur if r["blur_px"] == d]
        val = [r for r in sel if r["valid"]]
        ratio = [r["k_speed"] / d for r in val]
        print(f"   d = {d:4g} px: {_bound(len(val), len(sel))}"
              + (f"; K * speed / d median {np.median(ratio):.3f} (range {min(ratio):.3f} to {max(ratio):.3f}), "
                 f"direction error up to {max(r['angle_err'] for r in val):.1f} deg, within {TOL_ANGLE:g} deg: "
                 f"{sum(r['angle_err'] <= TOL_ANGLE for r in val)}" if val else ""))
    long = [r for r in blur if r["blur_px"] >= 4]
    if long:
        val = [r for r in long if r["valid"]]
        rest = [r for r in long if not r["valid"]]
        print(f"   d >= 4 px: {len(val)} of {len(long)} valid; direction error of the valid ones within "
              f"{TOL_ANGLE:g} deg: {sum(r['angle_err'] <= TOL_ANGLE for r in val)}, within 5 deg: "
              f"{sum(r['angle_err'] <= 5 for r in val)}; of the {len(rest)} not valid: quality below the threshold "
              f"{sum(r['quality'] < r['threshold'] for r in rest)}, out of range "
              f"{sum(not r['in_range'] for r in rest)}, isotropic {sum(bool(r['isotropic']) for r in rest)}")

    raw, rws = res["raw"], res["raw_stats"]
    if raw:
        noisy = {int(r["seed"]): r for r in res["null"] if r["kind"] == "surface"}
        top = rws["top"]
        twin = noisy.get(int(top["seed"]))
        print(f"\n7. Camera texture without the added noise ({len(raw)} crops at the positions of the surface crops "
              f"of parts 1 and 2):")
        print(f"   {_bound(rws['valid'], rws['n'])}; largest k {rws['k_max']:.3f} (n = {int(top['n'])}, "
              f"L = {int(top['L'])}, seed {int(top['seed'])}, {top['where']}, quality {top['quality']:.4f}, "
              f"speed {top['est_speed']:.2f} at {top['est_angle']:.1f} deg)"
              + (f", the same crop with the added noise: k {twin['k']:.3f}" if twin else "")
              + f"; 99th percentile {rws['k_p99']:.3f}; " + ", ".join(f"k >= {c:g}: {rws[f'above_{c:g}']}"
                                                                     for c in CONSTANTS))
        sel = [r for r in raw if int(r["seed"]) in noisy]
        print(f"   the same {len(sel)} crops with the added noise: largest k "
              f"{max((noisy[int(r['seed'])]['k'] for r in sel), default=math.nan):.3f}, "
              f"{sum(noisy[int(r['seed'])]['k'] >= NULL_K for r in sel)} with k >= {NULL_K:g}")
        print("   by size, largest k: " + ", ".join(
            f"{n} px {max(r['k'] for r in raw if int(r['n']) == n):.3f}" for n in SIZES
            if any(int(r["n"]) == n for r in raw)))


def replot():
    res = analyse(*(common.read_csv(n) for n in CSV))
    plot(res)
    summary(res)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--quick", action="store_true",
                   help="small sample; writes to experiments/results_quick unless --out is given")
    p.add_argument("--out", type=Path, metavar="DIR",
                   help="folder for the CSV files, the figure goes to DIR/figures (default: experiments/results "
                        "and docs/figures; with --quick experiments/results_quick)")
    p.add_argument("--jobs", type=int, default=os.cpu_count() or 1,
                   help="worker processes (default: number of CPUs)")
    p.add_argument("--replot", action="store_true",
                   help="redraw the figure and print the numbers from the CSV files in the output folder")
    p.add_argument("--recordings", type=Path, metavar="PATH",
                   help="folder with vid1.mat to vid4.mat (default: data/recordings)")
    a = p.parse_args(argv)
    if a.jobs < 1:
        p.error("--jobs must be at least 1")
    folder = a.recordings if a.recordings is not None else common.RECORDINGS
    out = (a.out if a.out is not None else common.QUICK if a.quick else common.COMMITTED).resolve()
    if a.quick and out == common.COMMITTED:
        p.error("--quick must not write to experiments/results, which holds the full-run results; "
                "leave out --out or choose another folder")
    if a.quick and out / "figures" == common.README_FIGURES:
        p.error("--quick must not write to docs/figures, which holds the figures of the full run; "
                "leave out --out or choose another folder")
    if not a.replot:
        missing = [f for f in FILES if not (folder / f).is_file()]
        if missing:
            p.error(f"--recordings {folder}: {', '.join(missing)} not found")
    common.set_output(out)
    import matplotlib
    matplotlib.use("Agg")
    common.FIGURES.mkdir(parents=True, exist_ok=True)
    print(f"CSV files: {common.RESULTS}\nfigures:   {common.FIGURES}")
    if a.replot:
        missing = common.missing_csv(CSV)
        if missing:
            sys.exit(f"--replot: {', '.join(missing)} not found in {common.RESULTS}")
        replot()
        return
    with Pool(a.jobs, initializer=common.init_worker) as pool:
        run(pool, a.quick, folder.resolve())


if __name__ == "__main__":
    main()
