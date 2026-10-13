"""Comparison with simple single-image baselines on the images of the noise
study: cepstrum, autocorrelation, corrected peak spacing, and phase
correlation between two single exposures as an upper reference.

    python experiments/baselines.py [--quick] [--out DIR] [--replot] [--jobs N]

The images are those of noise.py: random velocities (speed log-uniform 1 to
32 px per step, any direction), 9 exposures (L = 4), 256 x 256 px, random
and camera texture, 7 noise levels from 0 to 2 times the standard deviation
of the image, 60 images per texture and noise level (--quick: 4 levels, 12
images each). The cases are generated with the same seeds as in noise.py, so
the fftvel estimates here are those of noise.csv. Each method sees the same
pulsed image:

* fftvel: estimate_velocity with the defaults (the estimate of noise.csv).
* spacing: the direction of fftvel and speed = npad / spacing of the main
  lines in the fftvel profile (the 2025 student approach, done properly).
* cepstrum: strongest peak of the power cepstrum.
* autocorr: strongest side peak of the autocorrelation.
* phasecorr: phase correlation between two consecutive single exposures of
  the same motion, each with noise at the same noise/signal ratio as the
  pulsed image, with the peak refined by upsampling. It needs two frames,
  so it is an upper reference, not a competitor. The phasecorr rows also
  hold the result of fftvel.phase_correlation on the same two frames. It
  now refines its peak by upsampling as well; its earlier parabolic peak
  fit was off by about 0.1 px for shifts between pixels.

A result is correct when the direction is within 2 deg and the speed within
5 % of the true values. Only fftvel has a validity flag; the baselines
always return a velocity (or none when they find no peak, which counts as
not correct). The parameters of the baselines and the sub-pixel refinement
of each peak (a parabola for cepstrum and autocorrelation, upsampling for
the phase correlation) were chosen on separate images with other seeds
(rng 77, seeds from 200000), not on the images of the study.

Writes baselines.csv (one row per image and method) and baselines.png; the
summary prints the share correct, the median errors of the correct results
and, for fftvel, the shares valid and wrong but valid, the share correct by
speed and the precision of the phase correlation without noise. A full run
needs about 27 min of CPU time, so about 7 to 10 min on 4 cores. --quick
writes to experiments/results_quick; --replot redraws the figure and prints
the table from the CSV file in the output folder.

Columns of baselines.csv: the case (texture, speed, angle, L, noise, seed),
method, est_angle, est_speed, angle_err, speed_err, correct; on fftvel rows
also valid and quality; on phasecorr rows also vector_err (px, signed shift
against the true v), and vector_err_parabolic and correct_parabolic for
fftvel.phase_correlation (named after its earlier parabolic peak fit).
"""

from __future__ import annotations

import argparse
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy import signal

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import common  # noqa: E402
from common import AQUA, BLUE, INK2, ORANGE  # noqa: E402
from fftvel import angle_difference, estimate_velocity, phase_correlation, simulate_pulsed_image  # noqa: E402

CSV = ["baselines.csv"]
YELLOW = "#eda100"           # categorical slot 4 of the reference palette (common.py has slots 1-3)

METHODS = ("fftvel", "spacing", "cepstrum", "autocorr", "phasecorr")
LABELS = {"fftvel": "fftvel (template fit)", "spacing": "corrected peak spacing", "cepstrum": "cepstrum",
          "autocorr": "autocorrelation", "phasecorr": "phase correlation, 2 frames (reference)"}
STYLE = {"fftvel": (BLUE, "o", "-"), "spacing": (YELLOW, "D", "-"), "cepstrum": (ORANGE, "s", "-"),
         "autocorr": (AQUA, "^", "-"), "phasecorr": (INK2, "v", "--")}
TEXTURES = (("random", "random texture"), ("surface", "camera texture"))

# Parameters of the baselines (chosen on images with other seeds than the study)
PAD = 2                      # zero padding of the FFT, as in fftvel (linear autocorrelation needs 2)
MIN_LAG = 1.5                # px; cepstrum and autocorrelation peaks closer to 0 are not used
POWER_FLOOR = 1e-6           # cepstrum: log(P + POWER_FLOOR * mean(P)), only against log(0)
PROMINENCE_K = 2.0           # spacing: main lines stand out by this many robust standard deviations
MAD_SCALE = 1.4826           # robust standard deviation = MAD_SCALE * median absolute deviation
LINE_TOL = 0.2               # spacing: peaks further than this many spacings from a multiple are dropped
NOISE_STREAM = 1             # second seed of the noise added to the two single exposures
UPSAMPLE = 16                # phase correlation: peak located on a grid of 1/UPSAMPLE px around the maximum

ANGLE_TOL, SPEED_TOL = 2.0, 0.05     # correct: direction within 2 deg, speed within 5 %
SPEED_BANDS = ((1, 2), (2, 4), (4, 8), (8, 16), (16, 32))           # px per step, for the summary


# The images

def cases(quick=False):
    """The cases of noise.py, in the same order and with the same seeds."""
    levels = [0, 0.5, 1, 2] if quick else [0, 0.25, 0.5, 0.75, 1, 1.5, 2]
    per_level = 12 if quick else 60
    rng = np.random.default_rng(2025)
    out = []
    for texture in ("random", "surface"):
        for nl in levels:
            for _ in range(per_level):
                s = float(np.exp(rng.uniform(np.log(1), np.log(32))))
                a = float(rng.uniform(0, 180))
                out.append((texture, s, a, 4, nl, 5000 + len(out)))
    return out


def true_velocity(case):
    """Velocity ``(vx, vy)`` in px per step of a case (x along columns, y down), as in common.simulate_case."""
    texture, speed, angle, L, noise, seed = case
    a = np.deg2rad(angle)
    return speed * np.cos(a), -speed * np.sin(a)


def simulate(case):
    """Pulsed image and single exposures of one case; ``bp`` equals ``common.simulate_case(case)``.

    Needs the camera texture that ``common.init_worker`` loads.
    """
    texture, speed, angle, L, noise, seed = case
    tex = common._texture if texture == "surface" else None
    return simulate_pulsed_image(common.N, true_velocity(case), L, texture=tex, noise=noise, rng=seed,
                                 return_frames=True)


# Helpers of the baselines

def _vertex(ym, y0, yp):
    """Offset of the vertex of the parabola through three samples (0 if it is not a maximum)."""
    den = ym - 2 * y0 + yp
    return 0.5 * (ym - yp) / den if den < 0 else 0.0


def _angle_speed(dx, dy):
    """Direction in [0, 180) deg (counter-clockwise, y up) and length of a shift in image coordinates."""
    if not (np.isfinite(dx) and np.isfinite(dy)):
        return float("nan"), float("nan")
    return float(np.degrees(np.arctan2(-dy, dx)) % 180), float(np.hypot(dx, dy))


def _power(bp):
    """|FFT|^2 of the mean-free, Hann-windowed square image, zero padded to PAD times its size."""
    x = np.asarray(bp, dtype=float)
    n = x.shape[0]
    w = np.hanning(n)
    return np.abs(np.fft.fft2((x - x.mean()) * np.outer(w, w), (PAD * n, PAD * n))) ** 2


def _strongest_peak(R, rmax):
    """Shift ``(dx, dy)`` of the strongest local maximum of the lag map ``R`` (zero lag at index 0).

    Only 3 x 3 local maxima at lags from MIN_LAG to ``rmax`` px count. The
    position is refined by a parabola along each axis; for cepstrum and
    autocorrelation this gave more correct results on the tuning images
    than the upsampling of :func:`_upsampled_peak`. Returns NaNs when there
    is no such maximum.
    """
    R = np.fft.fftshift(R)
    m = R.shape[0]
    c = m // 2
    k = np.arange(m) - c
    lag = np.hypot(*np.meshgrid(k, k))
    centre = R[1:-1, 1:-1]
    peak = np.zeros(R.shape, dtype=bool)
    peak[1:-1, 1:-1] = np.all([centre > R[1 + dy:m - 1 + dy, 1 + dx:m - 1 + dx]
                               for dy in (-1, 0, 1) for dx in (-1, 0, 1) if dy or dx], axis=0)
    peak &= (lag >= MIN_LAG) & (lag <= rmax)
    if not peak.any():
        return float("nan"), float("nan")
    idx = np.flatnonzero(peak)
    r, q = divmod(int(idx[np.argmax(R.flat[idx])]), m)
    dy = r - c + _vertex(R[r - 1, q], R[r, q], R[r + 1, q])
    dx = q - c + _vertex(R[r, q - 1], R[r, q], R[r, q + 1])
    return dx, dy


def _upsampled_peak(F, dx0, dy0):
    """Sub-pixel maximum of the lag map ``real(ifft2(F))`` next to the integer lag ``(dx0, dy0)``.

    The lag map is evaluated by a matrix DFT of ``F`` on a grid of
    1/UPSAMPLE px within 1 px of ``(dx0, dy0)``, the band-limited
    interpolation of the samples (the upsampled cross-correlation of
    Guizar-Sicairos, Thurman and Fienup, 2008), and the maximum on that grid
    is refined by a parabola along each axis.
    """
    my, mx = F.shape
    u = np.arange(-UPSAMPLE, UPSAMPLE + 1) / UPSAMPLE
    ey = np.exp(2j * np.pi * np.outer(dy0 + u, np.fft.fftfreq(my) * my) / my)
    ex = np.exp(2j * np.pi * np.outer(np.fft.fftfreq(mx) * mx, dx0 + u) / mx)
    R = np.real(ey @ F @ ex)
    r, q = np.unravel_index(int(np.argmax(R)), R.shape)
    last = 2 * UPSAMPLE
    oy = _vertex(R[r - 1, q], R[r, q], R[r + 1, q]) if 0 < r < last else 0.0
    ox = _vertex(R[r, q - 1], R[r, q], R[r, q + 1]) if 0 < q < last else 0.0
    return float(dx0 + u[q] + ox / UPSAMPLE), float(dy0 + u[r] + oy / UPSAMPLE)


# The baselines

def cepstrum_velocity(bp):
    """Velocity ``(angle_deg, speed)`` from the power cepstrum of the pulsed image.

    The log power spectrum of K shifted copies is that of the texture plus
    log|sum exp(-2 pi i f.v k)|^2, whose inverse FFT has positive peaks at
    the multiples of v, the largest at v. The power cepstrum (inverse FFT of
    log|FFT|^2) is the classical tool for echoes and repeated copies in a
    signal; here it is taken of the Hann-windowed image, and the strongest
    local maximum at least MIN_LAG px from 0 and at most n/2 (the largest
    speed fftvel measures) gives v, up to its sign.
    """
    P = _power(bp)
    C = np.real(np.fft.ifft2(np.log(P + POWER_FLOOR * P.mean())))
    return _angle_speed(*_strongest_peak(C, np.shape(bp)[0] / 2))


def autocorrelation_velocity(bp):
    """Velocity ``(angle_deg, speed)`` from the autocorrelation of the pulsed image.

    The autocorrelation of K copies shifted by v is that of the texture
    repeated at the multiples m*v with weights K - |m|, so the first side
    peak at v is the strongest; looking for repeats in the autocorrelation
    is the usual first approach to a periodic or repeated pattern. The
    strongest local maximum at least MIN_LAG px from 0 and at most n/2 is
    taken. A local maximum cannot lie on the falling flank of the central
    lobe: when v is shorter than the width of the texture's own
    autocorrelation, the peak at v merges with the lobe and a later multiple
    (2v, 3v) or a random peak is found instead.
    """
    R = np.real(np.fft.ifft2(_power(bp)))
    return _angle_speed(*_strongest_peak(R, np.shape(bp)[0] / 2))


def spacing_speed(est):
    """Speed from the spacing of the main lines in the profile along the fftvel direction.

    ``est`` is the result of ``estimate_velocity(bp, diagnostics=True)``; its
    profile is the mean log spectrum along the direction found, at
    frequency bins ``j``. The main lines of the Dirichlet pattern lie at
    ``j = m * npad / speed``. This is the approach of the 2025 student
    project (speed = 2N / mean spacing of the peaks found with findpeaks)
    with its faults removed: the smooth texture spectrum is taken out by a
    cubic in ``j``; a peak counts as a main line only if both its
    prominence and its height above the median exceed PROMINENCE_K robust
    standard deviations (MAD_SCALE * median absolute deviation) of the
    detrended profile (the side lobes between main lines are about as
    prominent as the main lines in a clean image, but lower); each peak is
    refined by a parabola; the spacing is the least-squares fit of
    ``j = m * spacing`` to the peaks, with the line at ``j = 0`` that is
    always there, the median gap as start value, and peaks more than
    LINE_TOL spacings from a multiple left out. Returns NaN when no main
    line is found. The first main line must lie in the profile, so speeds
    below npad / max(j), about 2.2 px for 256 px images, are out of reach.
    """
    j = np.asarray(est.diag["profile_j"], dtype=float)
    g = np.asarray(est.diag["profile"], dtype=float)
    x = j / j.max()
    X = np.column_stack([x**3, x**2, x, np.ones_like(x)])
    res = g - X @ np.linalg.lstsq(X, g, rcond=None)[0]
    med = np.median(res)
    sigma = MAD_SCALE * np.median(np.abs(res - med))
    idx, _ = signal.find_peaks(res, prominence=PROMINENCE_K * sigma)
    idx = idx[res[idx] - med >= PROMINENCE_K * sigma]
    if idx.size == 0:
        return float("nan")
    pk = np.array([j[i] + _vertex(res[i - 1], res[i], res[i + 1]) for i in idx])
    s0 = np.median(np.diff(np.concatenate([[0.0], pk])))
    m = np.rint(pk / s0)
    keep = (m >= 1) & (np.abs(pk / s0 - m) <= LINE_TOL)
    if not keep.any():
        return float("nan")
    spacing = np.sum(m[keep] * pk[keep]) / np.sum(m[keep] ** 2)
    return float(est.npad / spacing)


def single_exposures(case, frames):
    """The two consecutive single exposures L and L+1 of a case, with noise, for the phase correlation.

    ``frames`` are those of :func:`simulate`. Each frame gets white noise
    with the noise/signal ratio of the case, relative to its own standard
    deviation, as the pulsed image does (rng ``[seed, NOISE_STREAM]``).
    """
    texture, speed, angle, L, noise, seed = case
    rng = np.random.default_rng([seed, NOISE_STREAM])
    pair = []
    for q in (L, L + 1):
        f = frames[:, :, q]
        if noise > 0:
            f = f + noise * f.std(ddof=1) * rng.standard_normal(f.shape)
        pair.append(f)
    return pair


def phase_correlation_shift(a, b):
    """Shift ``(dx, dy)`` that moves the content of ``a`` onto ``b`` (upper reference, two images).

    Phase correlation, the standard way to register two images by the shift
    theorem of the Fourier transform, computed as fftvel.phase_correlation
    does (mean-free, Hann-windowed images, normalised cross-power spectrum),
    with the peak refined by this study's own upsampling
    (:func:`_upsampled_peak`). It was written when fftvel.phase_correlation
    still fitted a parabola through three samples of the peak, which is
    about one pixel wide; that fit was off by up to about 0.17 px for shifts
    between pixels, more than 5 % of a speed below 3 px. For consecutive
    exposures the shift is v itself, with its sign.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    h, w = a.shape
    win = np.outer(np.hanning(h), np.hanning(w))
    X = np.fft.fft2((b - b.mean()) * win) * np.conj(np.fft.fft2((a - a.mean()) * win))
    X = X / np.maximum(np.abs(X), np.finfo(float).eps)
    i, j = np.unravel_index(int(np.argmax(np.real(np.fft.ifft2(X)))), X.shape)
    return _upsampled_peak(X, j - w if j > w / 2 else j, i - h if i > h / 2 else i)


# One image

def _errors(case, est_angle, est_speed):
    """Direction error (deg, modulo 180), relative speed error and correct (0 or 1, 0 for NaN)."""
    texture, speed, angle, L, noise, seed = case
    angle_err = abs(angle_difference(est_angle, angle))
    speed_err = abs(est_speed / speed - 1)
    return angle_err, speed_err, int(bool(angle_err <= ANGLE_TOL and speed_err <= SPEED_TOL))


def _row(case, method, est_angle, est_speed, **extra):
    texture, speed, angle, L, noise, seed = case
    angle_err, speed_err, correct = _errors(case, est_angle, est_speed)
    return dict(texture=texture, speed=speed, angle=angle, L=L, noise=noise, seed=seed, method=method,
                est_angle=est_angle, est_speed=est_speed, angle_err=angle_err, speed_err=speed_err,
                correct=correct, **extra)


def _vector_err(case, dx, dy):
    vx, vy = true_velocity(case)
    return float(np.hypot(dx - vx, dy - vy))


def evaluate(case):
    """Rows of all methods for one case, in the order of METHODS (top-level for Windows spawn)."""
    bp, frames = simulate(case)
    est = estimate_velocity(bp, L=case[3], diagnostics=True)
    a, b = single_exposures(case, frames)
    dx, dy = phase_correlation_shift(a, b)
    pdx, pdy, _ = phase_correlation(a, b)
    return [
        _row(case, "fftvel", est.angle_deg, est.speed, valid=int(est.valid), quality=est.quality),
        _row(case, "spacing", est.angle_deg, spacing_speed(est)),
        _row(case, "cepstrum", *cepstrum_velocity(bp)),
        _row(case, "autocorr", *autocorrelation_velocity(bp)),
        _row(case, "phasecorr", *_angle_speed(dx, dy), vector_err=_vector_err(case, dx, dy),
             vector_err_parabolic=_vector_err(case, pdx, pdy),
             correct_parabolic=_errors(case, *_angle_speed(pdx, pdy))[2]),
    ]


def run(pool, quick):
    rows = [r for rs in common.run_cases(cases(quick), pool, evaluate) for r in rs]
    common.write_csv("baselines.csv", rows)
    plot(rows)
    return rows


# Table and figure

def _median(values):
    return float(np.median(values)) if len(values) else float("nan")


def table(rows):
    """Metrics per method, texture and noise level (noise 'all': every level).

    Shares in %: correct, and for fftvel valid and wrong but valid (of all
    images); median direction error (deg) and speed error (%) of the correct
    results.
    """
    out = []
    for texture, _ in TEXTURES:
        levels = sorted({r["noise"] for r in rows if r["texture"] == texture})
        for nl in [*levels, "all"]:
            for method in METHODS:
                sel = [r for r in rows if r["texture"] == texture and r["method"] == method
                       and (nl == "all" or r["noise"] == nl)]
                if not sel:
                    continue
                ok = [r for r in sel if r["correct"]]
                t = dict(texture=texture, noise=nl, method=method, images=len(sel),
                         correct=100 * len(ok) / len(sel),
                         angle_err=_median([r["angle_err"] for r in ok]),
                         speed_err=100 * _median([r["speed_err"] for r in ok]))
                if method == "fftvel":
                    t["valid"] = 100 * np.mean([r["valid"] for r in sel])
                    t["wrong_valid"] = 100 * np.mean([r["valid"] * (1 - r["correct"]) for r in sel])
                out.append(t)
    return out


def summary(rows):
    """Print the table the docs quote."""
    t = table(rows)
    print("\nBaselines: share correct in %, [median direction error in deg / speed error in % of the correct]")
    head = "".join(f"{m:>21s}" for m in METHODS)
    for texture, _ in TEXTURES:
        print(f"{texture:8s} {'noise':>5s}{head}   fftvel valid  wrong+valid")
        for nl in dict.fromkeys(x["noise"] for x in t if x["texture"] == texture):
            cells = {x["method"]: x for x in t if x["texture"] == texture and x["noise"] == nl}
            line = "".join(f"{c['correct']:7.1f} [{c['angle_err']:5.2f}/{c['speed_err']:5.2f}]"
                           for c in (cells[m] for m in METHODS))
            f = cells["fftvel"]
            level = "all" if nl == "all" else f"{nl:4.2f}"
            print(f"{'':8s} {level:>5s}{line}   {f['valid']:12.1f} {f['wrong_valid']:12.1f}")
    print("By true speed, both textures and all noise levels: share correct in %")
    print(f"{'speed px':>14s}{''.join(f'{m:>11s}' for m in METHODS)}")
    for lo, hi in SPEED_BANDS:
        sel = [r for r in rows if lo <= r["speed"] < hi]
        if sel:
            shares = [100 * np.mean([r["correct"] for r in sel if r["method"] == m]) for m in METHODS]
            print(f"{lo:6g} to {hi:<4g}{''.join(f'{s:11.1f}' for s in shares)}")
    pc = [r for r in rows if r["method"] == "phasecorr"]
    clean = [r for r in pc if r["noise"] == 0]
    if clean:
        err = np.array([r["vector_err"] for r in clean])
        par = np.array([r["vector_err_parabolic"] for r in clean])
        print(f"Phase correlation without noise, median / largest error of the shift in px: "
              f"{np.median(err):.4f} / {err.max():.4f} with the upsampled peak, "
              f"{np.median(par):.4f} / {par.max():.4f} with fftvel.phase_correlation")
    if pc:
        print(f"Phase correlation correct, all images: {100 * np.mean([r['correct'] for r in pc]):.1f} % "
              f"with the upsampled peak, {100 * np.mean([r['correct_parabolic'] for r in pc]):.1f} % "
              f"with fftvel.phase_correlation")
    n = len({(r["texture"], r["noise"], r["seed"]) for r in rows})
    print(f"({n} images; the phase correlation uses two single exposures, the others one pulsed image)")


def plot(rows):
    fig, axes = common.figure(2, 9.6)
    for ax, (texture, label) in zip(axes, TEXTURES):
        levels = sorted({r["noise"] for r in rows if r["texture"] == texture})
        for method in METHODS:
            colour, marker, ls = STYLE[method]
            share = [100 * np.mean([r["correct"] for r in rows if r["texture"] == texture
                                    and r["method"] == method and r["noise"] == nl]) for nl in levels]
            ax.plot(levels, share, ls=ls, marker=marker, color=colour, lw=2, ms=5, label=LABELS[method],
                    zorder=3 if method == "fftvel" else 2)
        common.style(ax, "noise / signal (standard deviation)", "correct estimates (%)")
        ax.set_ylim(0, 104)
        ax.text(0.98, 0.97, label, transform=ax.transAxes, ha="right", va="top", color=INK2, fontsize=10)
    axes[0].legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="lower left")
    per = len([r for r in rows if r["method"] == "fftvel" and r["texture"] == "random" and r["noise"] == 0])
    common.title(fig, "fftvel against simple single-image baselines",
                 f"images of the noise study, {per} per texture and noise level, 1 to 32 px per flash; "
                 f"correct = within 2 deg and 5 %; validity not applied")
    common.save(fig, "baselines.png")


def replot():
    rows = common.read_csv("baselines.csv")
    plot(rows)
    summary(rows)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--quick", action="store_true",
                   help="small sample; writes to experiments/results_quick unless --out is given")
    p.add_argument("--out", type=Path, metavar="DIR",
                   help="folder for the CSV file, the figure goes to DIR/figures (default: experiments/results "
                        "and docs/figures; with --quick experiments/results_quick)")
    p.add_argument("--jobs", type=int, default=os.cpu_count(),
                   help="worker processes (default: number of CPUs)")
    p.add_argument("--replot", action="store_true",
                   help="redraw the figure and print the table from the CSV file in the output folder")
    a = p.parse_args(argv)
    if a.jobs < 1:
        p.error("--jobs must be at least 1")
    out = (a.out if a.out is not None else common.QUICK if a.quick else common.COMMITTED).resolve()
    if a.quick and out == common.COMMITTED:
        p.error("--quick must not write to experiments/results, which holds the full-run results; "
                "leave out --out or choose another folder")
    if a.quick and out / "figures" == common.README_FIGURES:
        p.error("--quick must not write to docs/figures, which holds the figures of the full run; "
                "leave out --out or choose another folder")
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
        rows = run(pool, a.quick)
    summary(rows)


if __name__ == "__main__":
    main()
