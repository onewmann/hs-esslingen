"""Model violations: an independent forward model of a multi-flash camera
image and how the estimator fails when the image breaks its assumptions.

    python experiments/model_violations.py [--quick] [--out DIR] [--replot] [--jobs N]

The estimator fits the ideal model of simulate_pulsed_image: K equal,
instantaneous flashes at equal steps, a moving texture and nothing else. Good
results on that simulator only verify the implementation. This study renders
the images with its own forward model, which does not call
simulate_pulsed_image, and breaks one assumption at a time.

Forward model (one image, ``n`` = 256 px, ``L`` = 4, true K = 9):

* Scene: a square crop of ``n + 448`` px of data/surface_texture.png (real
  camera texture; random position, rotation by a multiple of 90 deg and
  mirroring) for even image numbers, a texture of fftvel.make_texture
  (correlation length 1.5 px, relative contrast 0.07 as the camera texture)
  for odd ones. The radiance is scaled to a mean of 1.
* Flashes: flash k shows the scene shifted by ``t_k * v`` with
  ``t_k = k - (K-1)/2``. The shift is a Fourier shift (phase ramp on the
  FFT of the crop, written here, not taken from fftvel); the window is the
  central ``n`` x ``n`` px, far enough from the crop border that no wrapped
  content reaches it. A cubic spline shift (scipy.ndimage.shift, order 3)
  was the first choice, but in the ideal case it moves the estimate by up
  to 1.36 % and 0.56 deg at speeds below 2 px per step (above: 0.21 % and
  0.17 deg), because the spline blurs each flash differently with the
  sub-pixel part of its shift. A camera samples the shifted scene without such a blur, so the
  Fourier shift is the closer model; the check below reports both.
* Flash duration: each flash is the mean of sub-steps at most 0.25 px apart
  that cover ``duty * v`` around ``t_k * v`` (the simulator uses a closed
  form instead).
* Camera: the image is scaled so that its mean is 40 % of full scale and
  clipped at full scale (saturation); shot noise is Poisson for the given
  full well, plus Gaussian read noise of 7 e- rms; quantisation rounds to
  ``2**bits - 1`` levels. The reference full well, 32000 e-, is that of the
  Sony IMX174 and IMX249 sensors (1936 x 1216 px, as in the Basler
  recordings in data/recordings). Twice the exposure has the same effect
  as twice the full well.

Factors, each swept with all else ideal (the first level is the ideal
image, except for the two K factors and the two with sensor noise):

* pulse_width: flash duration / flash step (duty) 0 to 0.9.
* timing_jitter: flash k fires at ``t_k + s * e_k`` with ``e_k`` standard
  normal and independent for each flash (the jitter of a clocked flash
  trigger), s = 0 to 20 % of the flash period. Its position is off by
  ``s * e_k * |v|`` px along v.
* energy_jitter: flash energy ``1 + s * e_k``, s = 0 to 40 %.
* shot_noise: full well 100000 down to 1000 e- at 40 % exposure (a smaller
  full well is the same as a darker exposure), read noise 7 e-.
* quantisation: 12, 10, 8, 6 and 4 bits (8 bits at a quarter of the
  exposure are about 6 bits at full exposure).
* offset: a constant offset (stray light, flare or a raised black level)
  that makes up the given share of the mean signal, with the total held at
  40 % of the 32000 e- full well and its shot and read noise. Without noise
  a constant offset changes nothing, because the estimator removes the
  mean. Ambient light on the moving surface is not constant: it shows the
  smeared scene (the factor continuous).
* continuous: share of the signal from continuous light, which shows the
  scene smeared over the whole exposure, here 10 flash steps centred on
  the 9 flashes.
* continuous_k: the same with an exposure of exactly 9 steps. The smear
  over K steps has its spectral zeros where the flash pattern has them.
* defocus: disc point spread function of radius 1 to 8 px, without noise.
* defocus_noisy: radius 1 to 5 px with the noise of the 32000 e- sensor.
* edge: a static, sharp, straight brightness step (random direction,
  through a point near the centre) of 0.25 to 2 times the mean signal;
  from 1.5 on, the bright side reaches full scale and saturates.
* periodic: the scene is a rectangular tile of the texture (width w,
  height 3/4 w, made seamless by the periodic plus smooth decomposition of
  Moisan 2011), repeated; w = 256 down to 64 px.
* glints: two static Gaussian spots (sigma 1.5 px, random positions in the
  central half) with peak 0.3 to 100 times the mean signal; the sensor
  clips them at full scale, so bright glints become saturated discs.
* assumed_K: the image has 9 flashes, the estimator is told 7 to 11
  (``exposures``).
* true_K: the image has 7 to 11 flashes, the estimator assumes 9.

Each level has 32 images (16 per texture; --quick: 4), with speed
log-uniform from 0.5 to 30 px per step and direction uniform. Image i has
the same velocity, scene and random draws (jitter, edge, glint positions)
at every level of a factor; only the strength changes. Velocity, scene,
factor draws and sensor noise come from separate seed streams (see
:func:`stream`). Each image is
estimated with the default options and, where a matching option exists,
with it set to the true value (pulse_width: ``duty``, glints:
``clip_outliers=True``, true_K: ``exposures``; for assumed_K the true value
is the level 9). A result is correct when the direction is within 2 deg
and the speed within 5 %.

Check (run first): the forward model in the ideal case against
simulate_pulsed_image on the same scene and velocity, both estimated with
the defaults; also with the cubic spline shift, with duty 0.5 (sub-steps
against the closed form of the simulator) and with K = 8 (even K). The
simulator uses the same Fourier shift, so the tests also render images
directly, without an FFT: integer shifts by array slicing, and sub-pixel
shifts, flash duration and continuous light on a scene made of sinusoids,
whose shifted and smeared values are known in closed form.

Writes violations.csv (one row per image and estimator variant; the check
rows have factor ``check``) and violations.png (one panel per factor: share
correct and valid, and share wrong but valid, against the level). The
summary prints the check, the table per factor and level, the levels from
which a factor gives wrong but valid results, and per factor how many valid
results are blur-like (the columns quality_blur, blur_length and blur_like
of the estimate) and how many images have more than 1 % of their pixels
saturated (the column saturated). A full run takes about
20 min on 4 cores. --quick writes to experiments/results_quick; --replot
redraws the figure and prints the tables from the CSV file in the output
folder. Every image has fixed seeds, so a run gives the same CSV file on the
same machine.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from multiprocessing import Pool
from pathlib import Path
from typing import NamedTuple

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import common  # noqa: E402
from common import BLUE, GRID, INK, INK2, MUTED, ORANGE, SURFACE  # noqa: E402
from fftvel import angle_difference, estimate_velocity, make_texture, simulate_pulsed_image  # noqa: E402
from fftvel.core import SAT_WARN  # noqa: E402

CSV = ["violations.csv"]

L = 4                        # the estimator's default; true K = 2L+1 = 9
TRUE_K = 2 * L + 1
REACH = 224                  # px of scene beyond the window on each side (largest shift about 200 px)
SPEEDS = (0.5, 30.0)         # px per flash step, log-uniform
PER_LEVEL, PER_LEVEL_QUICK = 32, 4
SEED = 7100                  # first entry of every seed of this study
SCENE_KEY = 100              # seed key of the scene, distinct from the factor codes

RANDOM_CONTRAST = 0.07       # std / mean of the random texture (camera texture: 0.0695)
SUBSTEP = 0.25               # px; largest distance between the sub-steps of a flash or of continuous light
EXPOSURE = 0.4               # mean signal as a fraction of full scale
READ_NOISE = 7.0             # e- rms
REF_FULL_WELL = 32000        # e-, Sony IMX174 / IMX249
GLINT_SIGMA = 1.5            # px
TILE_ASPECT = 0.75           # height / width of the periodic tile
DISC_SUPERSAMPLE = 8         # samples per pixel and axis for the disc point spread function

ANGLE_TOL, SPEED_TOL = 2.0, 0.05             # correct: direction within 2 deg, speed within 5 %
AGREE_SPEED, AGREE_ANGLE = 1e-3, 0.1         # check: forward model against simulator, 0.1 % and 0.1 deg
CHECKS = ("fourier", "spline3", "duty_0.5", "K_8")


def _num(x):
    return f"{x:g}"


def _none(x):
    return "none" if x == 0 else f"{x:g}"


def _kilo(x):
    return "none" if x == 0 else (f"{x / 1000:g}k" if x >= 1000 else f"{x:g}")


class Factor(NamedTuple):
    name: str
    code: int                # seed stream of the factor (fixed, so the order of FACTORS does not matter)
    title: str
    xlabel: str
    levels: tuple
    fmt: object
    option: str | None       # estimator option set to the true value in the variant 'matched'


FACTORS = (
    Factor("pulse_width", 1, "Flash duration", "duration / flash step (duty)",
           (0, 0.1, 0.25, 0.5, 0.75, 0.9), _num, "duty"),
    Factor("timing_jitter", 2, "Flash timing jitter", "standard deviation (% of the flash period)",
           (0, 1, 2, 5, 10, 20), _num, None),
    Factor("energy_jitter", 3, "Flash energy jitter", "standard deviation (%)",
           (0, 2, 5, 10, 20, 40), _num, None),
    Factor("shot_noise", 4, "Shot and read noise", "full well (e-), mean at 40 %",
           (0, 100000, 32000, 10000, 3200, 1000), _kilo, None),
    Factor("quantisation", 5, "Quantisation", "bits, mean at 40 % of full scale",
           (0, 12, 10, 8, 6, 4), _none, None),
    Factor("offset", 6, "Constant offset (stray light)", "share of the signal (32k e- noise)",
           (0, 0.25, 0.5, 0.75, 0.9), _num, None),
    Factor("continuous", 7, "Continuous light, exposure 10 steps", "share of the signal",
           (0, 0.1, 0.25, 0.5, 0.75, 0.9), _num, None),
    Factor("continuous_k", 14, "Continuous light, exposure 9 steps", "share of the signal",
           (0, 0.25, 0.5, 0.75, 0.9), _num, None),
    Factor("defocus", 8, "Defocus", "disc radius (px)",
           (0, 1, 2, 3, 5, 8), _num, None),
    Factor("defocus_noisy", 15, "Defocus with sensor noise", "disc radius (px), 32k e- noise",
           (0, 1, 2, 3, 5), _num, None),
    Factor("edge", 9, "Static sharp edge", "step / mean signal",
           (0, 0.25, 0.5, 1, 1.5, 2), _num, None),
    Factor("periodic", 10, "Periodic texture", "tile width (px), height 3/4 of it",
           (0, 256, 192, 160, 128, 64), _none, None),
    Factor("glints", 11, "Two static glints", "peak / mean signal",
           (0, 0.3, 1, 3, 10, 100), _num, "clip_outliers"),
    Factor("assumed_K", 12, "Wrong K assumed", "K given to the estimator (true 9)",
           (7, 8, 9, 10, 11), _num, None),
    Factor("true_K", 13, "Wrong number of flashes", "flashes in the image (assumed 9)",
           (7, 8, 9, 10, 11), _num, "exposures"),
)
BY_NAME = {f.name: f for f in FACTORS}


def stream(i, *key):
    """Random generator of image ``i`` for one purpose, given by ``key``.

    numpy's SeedSequence pads the seed with zeros, so ``[a, b]`` and
    ``[a, b, 0]`` give the same stream. No key used here ends in 0: velocity
    ``()``, scene ``(SCENE_KEY,)``, draws of a factor ``(code,)`` and sensor
    noise of a level ``(code, level_index + 1)``.
    """
    return np.random.default_rng([SEED, i, *key])


# The scene

class Image(NamedTuple):
    index: int
    texture: str             # 'surface' (camera texture) or 'random' (make_texture)
    speed: float
    angle: float             # deg, counter-clockwise from +x, y up
    v: tuple                 # (vx, vy) in px per step, y down


def image_params(i):
    """Velocity and texture of image ``i``; the same at every level and factor."""
    rng = stream(i)
    speed = float(np.exp(rng.uniform(np.log(SPEEDS[0]), np.log(SPEEDS[1]))))
    angle = float(rng.uniform(0, 180))
    a = np.deg2rad(angle)
    return Image(i, "surface" if i % 2 == 0 else "random", speed, angle,
                 (speed * math.cos(a), -speed * math.sin(a)))


def scene(img, n):
    """Radiance of the scene, a ``(n + 2*REACH)`` square with mean 1.

    The camera texture needs ``common.init_worker`` (or ``common._texture``)
    to be loaded.
    """
    m = n + 2 * REACH
    rng = stream(img.index, SCENE_KEY)
    if img.texture == "surface":
        tex = common._texture
        if tex.ndim == 3:
            tex = tex[:, :, :3].mean(axis=2)
        oy, ox = rng.integers(0, np.array(tex.shape[:2]) - m + 1)
        crop = np.rot90(tex[oy:oy + m, ox:ox + m], int(rng.integers(4)))
        if rng.integers(2):
            crop = crop[:, ::-1]
        crop = np.array(crop, dtype=float)
        return crop / crop.mean()
    return 1 + RANDOM_CONTRAST * make_texture(m, rng=rng)


def periodic_component(u):
    """Periodic part of the periodic plus smooth decomposition (Moisan 2011).

    Removes the jumps at the borders of ``u`` so that tiling it leaves no
    seams; the mean is kept.
    """
    u = np.asarray(u, dtype=float)
    h, w = u.shape
    b = np.zeros_like(u)
    b[0, :] += u[-1, :] - u[0, :]
    b[-1, :] += u[0, :] - u[-1, :]
    b[:, 0] += u[:, -1] - u[:, 0]
    b[:, -1] += u[:, 0] - u[:, -1]
    den = (2 * np.cos(2 * np.pi * np.arange(h) / h)[:, None]
           + 2 * np.cos(2 * np.pi * np.arange(w) / w)[None, :] - 4)
    den[0, 0] = 1.0
    s = np.fft.fft2(b) / den
    s[0, 0] = 0.0
    return u - np.real(np.fft.ifft2(s))


def tiled(R, width):
    """``R`` replaced by a seamless rectangular tile from its centre, repeated."""
    m = R.shape[0]
    h = int(round(TILE_ASPECT * width))
    c = m // 2
    tile = periodic_component(R[c - h // 2:c - h // 2 + h, c - width // 2:c - width // 2 + width])
    return np.tile(tile, (m // h + 1, m // width + 1))[:m, :m]


def disc_otf(radius, m):
    """Transfer function on the ``m`` x ``m`` FFT grid of a disc of ``radius`` px (area-sampled)."""
    r = int(np.ceil(radius)) + 1
    c = np.arange(-r, r + 1)
    o = (np.arange(DISC_SUPERSAMPLE) + 0.5) / DISC_SUPERSAMPLE - 0.5
    x = (c[:, None] + o[None, :]).ravel()
    inside = x[:, None] ** 2 + x[None, :] ** 2 <= radius ** 2
    k = inside.reshape(c.size, DISC_SUPERSAMPLE, c.size, DISC_SUPERSAMPLE).mean(axis=(1, 3))
    psf = np.zeros((m, m))
    psf[np.ix_(c % m, c % m)] = k / k.sum()
    return np.fft.fft2(psf)


# Flashes and rendering

def flash_times(K, jitter=0.0, rng=None):
    """Times of the K flashes in flash periods, centred on 0.

    With ``jitter`` > 0 each flash is off its nominal time by ``jitter * e_k``
    (``e_k`` standard normal, independent for each flash), as for a trigger
    driven by a clock; the result is centred again.
    """
    t = np.arange(K) - (K - 1) / 2
    if jitter == 0:
        return t
    t = t + jitter * rng.standard_normal(K)
    return t - t.mean()


def substeps(length):
    """Sub-steps along a path of ``length`` px, at most SUBSTEP px apart: their positions as fractions
    of the path (centred, from -0.5 to 0.5) and their count."""
    count = max(1, int(np.ceil(length / SUBSTEP)))
    return ((np.arange(count) + 0.5) / count - 0.5), count


def phase_sum(shifts, weights, m):
    """``sum_j w_j exp(-2 pi i f . s_j)`` on the ``m`` x ``m`` FFT grid; ``shifts`` are ``(sy, sx)`` in px.

    Multiplying the FFT of the scene by it gives the weighted sum of the
    scene shifted by each ``s_j`` (Fourier shift).
    """
    shifts = np.atleast_2d(np.asarray(shifts, dtype=float))
    w = np.asarray(weights, dtype=float)
    f = np.fft.fftfreq(m)
    ey = np.exp(-2j * np.pi * shifts[:, :1] * f[None, :])
    ex = np.exp(-2j * np.pi * shifts[:, 1:] * f[None, :])
    return (ey * w[:, None]).T @ ex


def render(R, v, times, energies, n, duty=0.0, continuous=0.0, exposure_steps=None, otf=None):
    """Window of the mean of the flash images (Fourier shift), relative to the mean radiance 1.

    Flash k shows the scene shifted by ``times[k] * v``, weighted by
    ``energies[k] / K``; with ``duty > 0`` each flash is the mean of its
    sub-steps over ``duty * v``. A share ``continuous`` of the result is
    continuous light over the whole exposure, ``exposure_steps`` (an
    integer, default K) steps centred on the flashes, as sub-steps. ``otf``
    (transfer function of a blur) multiplies the FFT of the scene.
    """
    m = R.shape[0]
    vy, vx = v[1], v[0]
    K = len(times)
    T = np.fft.fft2(R)
    if otf is not None:
        T = T * otf
    t = np.asarray(times, dtype=float)
    D = phase_sum(np.column_stack([t * vy, t * vx]), np.asarray(energies) / K, m)
    speed = math.hypot(vx, vy)
    if duty > 0:
        tau, count = substeps(duty * speed)
        D = D * phase_sum(np.column_stack([duty * tau * vy, duty * tau * vx]), np.full(count, 1 / count), m)
    if continuous > 0:
        # the exposure as E intervals of one step, each covered by the same sub-steps: the sum over
        # all sub-steps is the product of the sum over the interval centres and the sum over one interval
        E = K if exposure_steps is None else int(exposure_steps)
        centres = np.arange(E) - (E - 1) / 2
        tau, count = substeps(speed)
        C = (phase_sum(np.column_stack([centres * vy, centres * vx]), np.full(E, 1 / E), m)
             * phase_sum(np.column_stack([tau * vy, tau * vx]), np.full(count, 1 / count), m))
        D = (1 - continuous) * D + continuous * C
    r0 = (m - n) // 2
    return np.real(np.fft.ifft2(T * D))[r0:r0 + n, r0:r0 + n]


def render_spline(R, v, times, n):
    """The same as :func:`render` without duty and continuous light, with cubic spline shifts.

    Uses scipy.ndimage.shift (order 3) on the spline coefficients of the
    scene, cut to the window and a border of 3 px: the same values as
    shifting the whole scene. Only used in the check.
    """
    from scipy import ndimage
    pad = 3
    m = R.shape[0]
    r0 = (m - n) // 2
    coeff = ndimage.spline_filter(R, order=3, mode="mirror")
    out = np.zeros((n, n))
    for t in times:
        sy, sx = t * v[1], t * v[0]
        iy, ix = int(np.floor(sy)), int(np.floor(sx))
        patch = coeff[r0 - iy - pad:r0 - iy + n + pad, r0 - ix - pad:r0 - ix + n + pad]
        moved = ndimage.shift(patch, (sy - iy, sx - ix), order=3, prefilter=False, mode="mirror")
        out += moved[pad:pad + n, pad:pad + n]
    return out / len(times)


def static_edge(n, rng):
    """Sharp straight step from 0 to 1 (pixel-area coverage), random direction, near the centre."""
    phi = rng.uniform(0, np.pi)
    cy, cx = n / 2 + rng.uniform(-n / 8, n / 8, 2)
    y, x = np.mgrid[0:n, 0:n] + 0.5
    return np.clip((x - cx) * np.cos(phi) + (y - cy) * np.sin(phi) + 0.5, 0, 1)


def glint_map(n, rng):
    """Two Gaussian spots of peak 1 at random positions in the central half of the window."""
    pos = rng.uniform(n / 4, 3 * n / 4, (2, 2))
    y, x = np.mgrid[0:n, 0:n] + 0.5
    return sum(np.exp(-((x - px) ** 2 + (y - py) ** 2) / (2 * GLINT_SIGMA ** 2)) for py, px in pos)


def sensor(rel, full_well=0, bits=0, rng=None):
    """Camera output for the relative image ``rel`` (mean about 1) and its noise / signal ratio.

    Mean at EXPOSURE of full scale; with ``full_well > 0`` Poisson shot
    noise and READ_NOISE; clipped to full scale; with ``bits > 0`` rounded to
    ``2**bits - 1`` levels (the result is then in digital numbers).
    """
    E = EXPOSURE * rel
    clean = np.clip(E, 0, 1)
    if full_well > 0:
        E = (rng.poisson(np.maximum(E, 0) * full_well) + READ_NOISE * rng.standard_normal(E.shape)) / full_well
    E = np.clip(E, 0, 1)
    if bits > 0:
        top = 2 ** bits - 1
        out = np.rint(E * top)
        E = out / top
    else:
        out = E
    nsr = 0.0
    if full_well > 0 or bits > 0:
        spread = np.std(clean)
        nsr = float(np.std(E - clean) / spread) if spread > 0 else float("nan")     # nan: flat image
    return out, nsr


def make_image(factor, level_index, i, n=None):
    """Pulsed image of image ``i`` at one level of a factor; returns ``(image, params, nsr)``."""
    n = common.N if n is None else n
    f = BY_NAME[factor]
    level = f.levels[level_index]
    img = image_params(i)
    R = scene(img, n)
    frng = stream(i, f.code)                        # the same draws at every level
    lrng = stream(i, f.code, level_index + 1)       # sensor noise of this level
    K = int(level) if factor == "true_K" else TRUE_K
    times = flash_times(K, level / 100 if factor == "timing_jitter" else 0.0, frng)
    energies = np.ones(K)
    if factor == "energy_jitter":
        energies = np.maximum(1 + level / 100 * frng.standard_normal(K), 0)
    if factor == "periodic" and level > 0:
        R = tiled(R, int(level))
    defocus = factor in ("defocus", "defocus_noisy") and level > 0
    otf = disc_otf(level, R.shape[0]) if defocus else None
    continuous = factor in ("continuous", "continuous_k")
    rel = render(R, img.v, times, energies, n, duty=level if factor == "pulse_width" else 0.0,
               continuous=level if continuous else 0.0,
               exposure_steps=K + 1 if factor == "continuous" else K, otf=otf)
    if factor == "edge":
        rel = rel + level * static_edge(n, frng)
    if factor == "glints":
        rel = rel + level * glint_map(n, frng)
    full_well = level if factor == "shot_noise" else 0
    if factor == "offset":
        rel = (1 - level) * rel + level
    if factor in ("offset", "defocus_noisy"):
        full_well = REF_FULL_WELL
    bits = int(level) if factor == "quantisation" else 0
    out, nsr = sensor(rel, full_well, bits, lrng)
    return out, img, nsr


def variants(factor, level):
    """Estimator options per variant: the defaults, and the matching option set to the true value."""
    f = BY_NAME[factor]
    base = {"L": L}
    if factor == "assumed_K":
        base["exposures"] = int(level)
    out = {"default": base}
    if f.option == "duty":
        out["matched"] = {**base, "duty": float(level)}
    elif f.option == "clip_outliers":
        out["matched"] = {**base, "clip_outliers": True}
    elif f.option == "exposures":
        out["matched"] = {**base, "exposures": int(level)}
    return out


# One image

def _scores(img, e):
    angle_err = abs(angle_difference(e.angle_deg, img.angle))
    speed_err = abs(e.speed / img.speed - 1)
    return dict(est_speed=e.speed, est_angle=e.angle_deg, quality=e.quality, quality_perp=e.quality_perp,
                isotropic=int(e.isotropic), in_range=int(e.in_range), valid=int(e.valid),
                quality_blur=e.quality_blur, blur_length=e.blur_length, blur_like=int(e.blur_like),
                saturated=e.saturated, angle_err=angle_err, speed_err=speed_err,
                correct=int(bool(angle_err <= ANGLE_TOL and speed_err <= SPEED_TOL)))


def evaluate(case):
    """Rows of one image at one level, one per estimator variant (top-level for Windows spawn)."""
    factor, level_index, i, n = case
    level = BY_NAME[factor].levels[level_index]
    bp, img, nsr = make_image(factor, level_index, i, n)
    rows = []
    for name, opts in variants(factor, level).items():
        e = estimate_velocity(bp, **opts)
        rows.append(dict(factor=factor, level=level, level_index=level_index, variant=name, image=i,
                         texture=img.texture, speed=img.speed, angle=img.angle, nsr=nsr, **_scores(img, e)))
    return rows


def check(case):
    """Row comparing the forward model with simulate_pulsed_image on one image (top-level)."""
    name, i, n = case
    n = common.N if n is None else n
    img = image_params(i)
    R = scene(img, n)
    K, duty = TRUE_K, 0.0
    if name == "duty_0.5":
        duty = 0.5
    if name == "K_8":
        K = 8
    times = flash_times(K)
    if name == "spline3":
        bp = render_spline(R, img.v, times, n)
    else:
        bp = render(R, img.v, times, np.ones(K), n, duty=duty)
    ref = simulate_pulsed_image(n, img.v, L=L, exposures=K, duty=duty, texture=R)
    opts = dict(L=L, exposures=K, duty=duty)
    e = estimate_velocity(sensor(bp)[0], **opts)
    s = estimate_velocity(ref, **opts)
    return [dict(factor="check", level=0, level_index=0, variant=name, image=i, texture=img.texture,
                 speed=img.speed, angle=img.angle, nsr=0.0, **_scores(img, e),
                 sim_speed=s.speed, sim_angle=s.angle_deg, sim_valid=int(s.valid),
                 dev_speed=e.speed / s.speed - 1, dev_angle=angle_difference(e.angle_deg, s.angle_deg))]


def cases(quick=False, n=None):
    """Check cases and factor cases, in the order of the CSV file."""
    per = PER_LEVEL_QUICK if quick else PER_LEVEL
    checks = [(name, i, n) for name in CHECKS for i in range(per)]
    sweeps = [(f.name, li, i, n) for f in FACTORS for li in range(len(f.levels)) for i in range(per)]
    return checks, sweeps


def run(pool, quick):
    checks, sweeps = cases(quick)
    rows = [r for rs in common.run_cases(checks, pool, check) for r in rs]
    print_check(rows)
    rows += [r for rs in common.run_cases(sweeps, pool, evaluate) for r in rs]
    common.write_csv("violations.csv", rows)
    plot(rows)
    return rows


# Tables and figure

def _median(values):
    return float(np.median(values)) if len(values) else float("nan")


def check_table(rows):
    """Largest deviation of the forward model from the simulator per check."""
    out = []
    for name in CHECKS:
        sel = [r for r in rows if r["factor"] == "check" and r["variant"] == name]
        if not sel:
            continue
        ds = [abs(r["dev_speed"]) for r in sel]
        da = [abs(r["dev_angle"]) for r in sel]
        out.append(dict(check=name, images=len(sel), max_speed=100 * max(ds), median_speed=100 * _median(ds),
                        max_angle=max(da), median_angle=_median(da),
                        valid_agree=int(sum(r["valid"] == r["sim_valid"] for r in sel)),
                        agree=bool(max(ds) <= AGREE_SPEED and max(da) <= AGREE_ANGLE)))
    return out


def print_check(rows):
    t = check_table(rows)
    print("\nCheck: forward model against simulate_pulsed_image, same scene and velocity "
          f"(agreement: within {100 * AGREE_SPEED:g} % and {AGREE_ANGLE:g} deg)")
    print(f"{'check':>10s} {'images':>6s} {'max speed %':>12s} {'median %':>9s} {'max angle':>10s} "
          f"{'median':>8s} {'same valid':>11s}")
    for c in t:
        print(f"{c['check']:>10s} {c['images']:6d} {c['max_speed']:12.2e} {c['median_speed']:9.2e} "
              f"{c['max_angle']:10.2e} {c['median_angle']:8.2e} {c['valid_agree']:5d} of {c['images']:<3d}"
              f"{'  agrees' if c['agree'] else ''}")
    return t


def table(rows):
    """Shares in % (correct and valid, wrong but valid, not valid) and median errors of the valid results."""
    out = []
    for f in FACTORS:
        sel_f = [r for r in rows if r["factor"] == f.name]
        for level in sorted({r["level"] for r in sel_f}, key=lambda x: _order(f, x)):
            for variant in ("default", "matched"):
                sel = [r for r in sel_f if r["level"] == level and r["variant"] == variant]
                if not sel:
                    continue
                val = [r for r in sel if r["valid"]]
                cv = sum(1 for r in val if r["correct"])
                out.append(dict(factor=f.name, level=level, variant=variant, images=len(sel),
                                correct_valid=100 * cv / len(sel), wrong_valid=100 * (len(val) - cv) / len(sel),
                                not_valid=100 * (len(sel) - len(val)) / len(sel),
                                speed_err=100 * _median([r["speed_err"] for r in val]),
                                angle_err=_median([r["angle_err"] for r in val]),
                                nsr=float(np.mean([r["nsr"] for r in sel]))))
    return out


def _order(f, level):
    """Position of a level in the sweep (levels of an older run that are not in FACTORS go last)."""
    return f.levels.index(level) if level in f.levels else len(f.levels) + level


def reference_level(f):
    """The level the other levels of a factor are compared with: K = 9, otherwise the first level."""
    return TRUE_K if f.name in ("assumed_K", "true_K") else f.levels[0]


def onsets(rows, t=None):
    """Wrong but valid results per factor and variant.

    The share at the reference level, and for each other level with a
    wrong but valid result its share, the median ratio of estimated to true
    speed and the median direction error of those results; the first level
    with 10 % or more.
    """
    t = table(rows) if t is None else t
    out = []
    for f in FACTORS:
        ref = reference_level(f)
        for variant in ("default", "matched"):
            sel = [x for x in t if x["factor"] == f.name and x["variant"] == variant]
            if not sel:
                continue
            base = [x["wrong_valid"] for x in sel if x["level"] == ref]
            levels = []
            for x in sel:
                if x["level"] == ref or x["wrong_valid"] == 0:
                    continue
                bad = [r for r in rows if r["factor"] == f.name and r["variant"] == variant
                       and r["level"] == x["level"] and r["valid"] and not r["correct"]]
                levels.append(dict(level=x["level"], share=x["wrong_valid"], count=len(bad),
                                   speed_ratio=_median([r["est_speed"] / r["speed"] for r in bad]),
                                   angle_err=_median([r["angle_err"] for r in bad])))
            ten = [lv["level"] for lv in levels if lv["share"] >= 10]
            out.append(dict(factor=f.name, variant=variant, reference=ref, ref_share=base[0] if base else 0.0,
                            levels=levels, first10=ten[0] if ten else None))
    return out


def summary(rows, check=True):
    """Print the check (unless ``check`` is False), the table and the onset of wrong but valid results."""
    if check:
        print_check(rows)
    t = table(rows)
    print("\nModel violations: shares in % of the images; median errors of the valid results")
    print(f"{'factor':>14s} {'level':>7s} {'variant':>8s} {'images':>6s} {'correct+valid':>14s} "
          f"{'wrong+valid':>12s} {'not valid':>10s} {'speed err %':>12s} {'angle err':>10s} {'noise/sig':>10s}")
    for x in t:
        f = BY_NAME[x["factor"]]
        print(f"{x['factor']:>14s} {f.fmt(x['level']):>7s} {x['variant']:>8s} {x['images']:6d} "
              f"{x['correct_valid']:14.1f} {x['wrong_valid']:12.1f} {x['not_valid']:10.1f} "
              f"{x['speed_err']:12.3f} {x['angle_err']:10.3f} {x['nsr']:10.3f}")
    print("\nWrong but valid results: share at the reference level (K = 9 or the first level); for each other "
          "level with such results\nits share, and the median speed ratio (estimated / true) and direction error "
          "of those results")
    for o in onsets(rows, t):
        f = BY_NAME[o["factor"]]
        where = "; ".join(f"{f.fmt(lv['level'])}: {lv['share']:.1f} %, x{lv['speed_ratio']:.3g}, "
                          f"{lv['angle_err']:.2g} deg" for lv in o["levels"]) or "none"
        print(f"{o['factor']:>14s} {o['variant']:>8s}: at {f.fmt(o['reference'])} {o['ref_share']:4.1f} %; {where}")
    n = len({(r["factor"], r["level"], r["image"]) for r in rows if r["factor"] != "check"})
    print(f"({n} images in {len(t)} rows of the table)")
    sweep = [r for r in rows if r["factor"] != "check"]
    if sweep and all("blur_like" in r and "saturated" in r for r in sweep):    # not in runs before these fields
        print_blur(sweep)


def blur_table(rows):
    """Per factor and variant, all levels: valid results, correct and wrong, and how many of them are
    blur-like (``reject_blur`` would make them not valid); images with more than SAT_WARN saturated pixels."""
    out = []
    for f in FACTORS:
        for variant in ("default", "matched"):
            sel = [r for r in rows if r["factor"] == f.name and r["variant"] == variant]
            if not sel:
                continue
            ok = [r for r in sel if r["valid"] and r["correct"]]
            bad = [r for r in sel if r["valid"] and not r["correct"]]
            out.append(dict(factor=f.name, variant=variant, images=len(sel), correct_valid=len(ok),
                            correct_blur=int(sum(r["blur_like"] for r in ok)), wrong_valid=len(bad),
                            wrong_blur=int(sum(r["blur_like"] for r in bad)),
                            saturated=int(sum(r["saturated"] > SAT_WARN for r in sel))))
    return out


def print_blur(rows):
    print(f"\nBlur check and saturation, all levels: valid results that are blur-like, and images with more than "
          f"{100 * SAT_WARN:g} % of their pixels saturated")
    print(f"{'factor':>14s} {'variant':>8s} {'images':>6s} {'correct+valid blur-like':>24s} "
          f"{'wrong+valid blur-like':>22s} {'saturated':>10s}")
    for x in blur_table(rows):
        print(f"{x['factor']:>14s} {x['variant']:>8s} {x['images']:6d} "
              f"{x['correct_blur']:13d} of {x['correct_valid']:<4d}  {x['wrong_blur']:11d} of {x['wrong_valid']:<4d}  "
              f"{x['saturated']:9d}")


def grid_figure(nrows, ncols, width, height):
    """Figure with a grid of panels and room for title and subtitle, as common.figure."""
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(nrows, ncols, figsize=(width, height), layout="constrained")
    fig.get_layout_engine().set(rect=(0, 0, 1, 1 - 0.62 / height))
    fig.patch.set_facecolor(SURFACE)
    return fig, axes


def grid_title(fig, text, sub=None):
    """common.title for a figure of any height (same distances from the top in inches)."""
    h = fig.get_figheight()
    fig.text(0.012 * 6.4 / fig.get_figwidth(), 1 - 0.147 / h, text, ha="left", va="center", color=INK,
             fontsize=12, fontweight="bold")
    if sub:
        fig.text(0.012 * 6.4 / fig.get_figwidth(), 1 - 0.40 / h, sub, ha="left", va="center", color=INK2,
                 fontsize=9)


def plot(rows):
    from matplotlib.lines import Line2D
    t = table(rows)
    ncols = 4
    nrows = -(-(len(FACTORS) + 1) // ncols)
    fig, axes = grid_figure(nrows, ncols, 12.8, 2.75 * nrows + 0.62)
    axes = axes.ravel()
    for ax, f in zip(axes, FACTORS):
        sel = [x for x in t if x["factor"] == f.name]
        levels = sorted({x["level"] for x in sel}, key=lambda lv: _order(f, lv))
        pos = np.arange(len(levels))
        for variant, ls, mk in (("default", "-", "o"), ("matched", "--", "s")):
            cells = {x["level"]: x for x in sel if x["variant"] == variant}
            if not cells:
                continue
            ax.plot(pos, [cells[lv]["correct_valid"] for lv in levels], ls=ls, marker=mk, color=BLUE, lw=2, ms=4)
            ax.plot(pos, [cells[lv]["wrong_valid"] for lv in levels], ls=ls, marker=mk, color=ORANGE, lw=2, ms=4)
        common.style(ax, f.xlabel, "images (%)" if list(axes).index(ax) % ncols == 0 else "")
        ax.set_xticks(pos, [f.fmt(lv) for lv in levels])
        ax.set_ylim(-3, 104)
        ax.set_title(f.title, loc="left", color=INK, fontsize=10)
        ax.xaxis.label.set_fontsize(9)
    key = axes[len(FACTORS)]
    key.axis("off")
    handles = [Line2D([], [], color=BLUE, lw=2, marker="o", ms=4, label="correct and valid"),
               Line2D([], [], color=ORANGE, lw=2, marker="o", ms=4, label="wrong but valid"),
               Line2D([], [], color=MUTED, lw=2, ls="-", label="default options"),
               Line2D([], [], color=MUTED, lw=2, ls="--", marker="s", ms=4,
                      label="option set to the true value\n(duty, clip_outliers, exposures)")]
    key.legend(handles=handles, loc="upper left", frameon=False, fontsize=9, labelcolor=INK2)
    c = {x["check"]: x for x in check_table(rows)}
    if "fourier" in c:
        key.text(0.0, 0.18, f"ideal case against simulate_pulsed_image:\n"
                            f"largest difference {c['fourier']['max_speed']:.1e} % and "
                            f"{c['fourier']['max_angle']:.1e} deg", transform=key.transAxes, color=INK2,
                 fontsize=8.5, va="top")
    for ax in axes[len(FACTORS) + 1:]:
        ax.axis("off")
    per = len([r for r in rows if r["factor"] == "edge" and r["level"] == 0 and r["variant"] == "default"])
    grid_title(fig, "When the image breaks the model: correct and wrong valid estimates",
               f"independent forward model, one violation at a time, {per} images per level (camera and random "
               f"texture), 0.5 to 30 px per step, 9 flashes, 256 px; correct = within 2 deg and 5 %")
    for ax in axes[:len(FACTORS)]:
        ax.axhline(0, color=GRID, lw=0.8, zorder=0)
    common.save(fig, "violations.png")


def replot():
    rows = common.read_csv("violations.csv")
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
                   help="redraw the figure and print the tables from the CSV file in the output folder")
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
    summary(rows, check=False)          # run printed the check before the sweeps


if __name__ == "__main__":
    main()
