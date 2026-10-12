"""Velocity estimation from a single multi-flash (pulsed-exposure) image.

The model, the ripple profile WR and the quality measure GEBA follow Chapter 3
of J. Horn's doctoral thesis (Universitaet Karlsruhe (TH), 2006; full reference
in README.md). The speed fit, the validity threshold, the isotropy check and the
blur check were added in this implementation.

This module mirrors the MATLAB/Octave package ``matlab/+fftvel`` line by
line, so that both give the same numbers for the same input (see
``tests/parity``). Only NumPy is needed.

Conventions
-----------
* Images are 2-D arrays indexed ``[row, column]``.
* Image coordinates: x along columns (to the right), y along rows (down).
* Angles are in degrees, counter-clockwise from +x as seen on screen (y up),
  in ``[0, 180)``. A single pulsed image cannot tell ``v`` from ``-v``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

# Public names; the helpers below (log_spectrum, mean_projection, ripple_profile,
# fit_speed, fit_blur, null_quality, dirichlet_template, box_template,
# parabolic_offset, hann_window, crop_square, dilate, saturated_share) are
# imported from fftvel.core by the tests and experiments.
__all__ = [
    "Estimate",
    "PhysicalVelocity",
    "angle_difference",
    "estimate_velocity",
    "pulsed_from_stack",
    "to_physical",
]

# Numeric constants of the algorithm. tests/constants.json holds the same names
# and values with the meaning of each, and matlab/+fftvel/+internal/constants.m
# returns them in MATLAB; tests in both languages check all three against each
# other. Every upper-case number or tuple of numbers in this module must be
# listed there.

# Defaults of estimate_velocity
DEFAULT_L = 4                # half-count of the thesis, K = 2L+1 = 9 exposures
DEFAULT_PAD_FACTOR = 2       # zero padding factor of the FFT
DEFAULT_MIN_QUALITY = 0.6    # lowest fit quality accepted for valid

# Input checks
MIN_SIZE = 16                # smallest analysed square (and roi) in px
MIN_EXPOSURES = 2            # smallest number of exposures K

# Log spectrum
LOG_FLOOR = 0.05             # S = log(|B| + LOG_FLOOR*median|B|), mean|B| if the median is 0
BLANK_FLOOR = 1.0            # floor instead when |B| is 0 everywhere (blank image)
CLIP_SIGMA = 6               # clip_outliers: robust standard deviations from the median ...
MAD_SCALE = 1.4826           # ... with robust standard deviation = MAD_SCALE * MAD
CLIP_RADIUS = 3              # ... and the pixels within this distance (px) of such a pixel

# Direction
ANGLE_STEP = 1.0             # angle grid 0, ANGLE_STEP, ..., 180 - ANGLE_STEP degrees
RFRAC = 0.9                  # fraction of the spectrum radius used for WR and the fit
WEIGHT_FLOOR = 1e-12         # mean_projection: lower bound of the bin weights (empty bins)

# Speed
PROFILE_START = 2            # first bin of the half profile used for the speed fit
SPEED_STEP = 0.005           # step of the logarithmic speed grid
VMIN_FACTOR = 1.5            # vmin = VMIN_FACTOR*npad/(K*max(jpos)): first zero well inside the profile
MIN_LINE_SPACING = 4         # vmax = npad/MIN_LINE_SPACING: motion lines at least this many bins apart
GUARD_FACTOR = 4             # the speed grid extends down to vmin/GUARD_FACTOR (guard band)
TREND_DEGREE = 3             # degree of the polynomial in jpos that absorbs the texture spectrum
TEMPLATE_FLOOR = 1e-3        # dirichlet_template and box_template return 0.5*log(P + TEMPLATE_FLOOR)
FLAT_MODEL = 1e-6            # fit_speed, fit_blur: score 0 where the detrended model is below this share of it

# Validity
NULL_K = 9                   # null_quality(m) = tanh(NULL_K/sqrt(m - NULL_OFFSET)) ...
NULL_OFFSET = 5              # ... with TREND_DEGREE + 1 points taken by the polynomial, 1 by the fit
INDEPENDENT_PAD = 2          # m = min(points, floor(INDEPENDENT_PAD*n*points/npad)) independent points
ISO_GAP = 0.15               # isotropic if quality - quality_perp is below this ...
ISO_LOG_SPEED = 0.12         # ... and |log(speed_perp/speed)| is below this
BLUR_STEP = 0.5              # px between the samples of the continuous blur in box_template
BLUR_GAP = 0.0               # blur_like if quality - quality_blur is below this (and the speed is low)

# Saturation
SAT_SHARE = 0.005            # saturated: least share of pixels at the maximum (unless it is the integer maximum)
SAT_WARN = 0.01              # command line, stream and experiments warn above this share

# Simulation (simulate.py) and grey conversion of colour input
SIM_MARGIN = 32              # extra texture size in px beyond n + (K-1+duty)*max|v|
DEFAULT_CORR_LEN = 1.5       # correlation length of the random texture in px
GREY_WEIGHTS = (0.2989, 0.5870, 0.1140)    # R, G, B weights of MATLAB rgb2gray


@dataclass
class Estimate:
    """Result of :func:`estimate_velocity`.

    Attributes
    ----------
    angle_deg : direction of motion in ``[0, 180)`` degrees, counter-clockwise
        from +x as seen on screen.
    speed : ``|v|`` in pixels per exposure step.
    vx, vy : ``speed * (cos(angle), -sin(angle))`` in image coordinates
        (y pointing down). ``-(vx, vy)`` is equally possible.
    geba : ``max(WR) / mean(WR)`` of the ripple profile.
    quality : correlation between the measured profile and the model.
    quality_perp : quality of the same fit perpendicular to the direction.
    in_range : False if the best speed lies in the guard band below the
        measurable range (below ``vmin``, see :func:`fit_speed`) or on the
        last point of the speed grid, its top end at about ``npad/4``.
    isotropic : True if the perpendicular fit is about as good
        (``quality - quality_perp < 0.15``) and gives about the same speed
        (``|log(speed_perp/speed)| < 0.12``). Motion makes the spectrum
        anisotropic; defocus blur and square repeating textures do not.
    threshold : quality needed for this image size, the larger of
        ``min_quality`` and ``tanh(9/sqrt(m-5))``. ``m`` counts independent
        profile points: ``min(points, floor(2*n*points/npad))``, so zero
        padding beyond ``pad_factor = 2`` does not lower the threshold.
    valid : ``quality >= threshold`` and ``in_range`` and not ``isotropic``
        (and not ``blur_like`` with ``reject_blur``).
    L : ``L`` as passed.
    exposures : number of exposures ``K`` used in the model (``2L+1``
        unless ``exposures`` was given).
    n : side of the analysed central square in pixels.
    npad : FFT length.
    quality_blur : quality of the same fit along ``angle_deg`` with a single
        continuous exposure (motion blur of a free length, no flashes) in
        place of the ``K`` flashes.
    blur_length : blur length in pixels of that fit; a blur of ``d`` px
        looks like ``K`` flashes about ``d/K`` apart.
    blur_like : True if the blur model fits better (``quality_blur >
        quality``) and the speed is below ``npad/max(profile_j)``, about
        2.2 px per step, so that the repeated lines of the flashes
        (``npad/speed`` bins apart), which alone tell flashes from a blur,
        lie beyond the fitted profile. Most single exposures of a moving
        surface are blur-like, but below that speed so are most multi-flash
        images with noise: there the flag says that the image cannot tell
        flashes from a blur, not that it is a blur. Only ``reject_blur``
        makes it count for ``valid``.
    saturated : share of the analysed pixels at the maximum of the square,
        before any processing (:func:`saturated_share`): for an integer
        image whose maximum is the largest value of its type, the share at
        that value; otherwise the share at the maximum if it is at least
        0.5 % and the square is not constant, else 0. Saturated areas and
        glints can give wrong results that are marked valid.
    diag : None, or with ``diagnostics=True`` a dict of intermediate
        results: ``spectrum`` (log spectrum, ``npad``-by-``npad``),
        ``theta`` and ``wr`` (angle grid in degrees and ripple profile WR),
        ``profile_j``, ``profile`` and ``model`` (bins of the half profile
        along ``angle_deg``, the profile and the fitted model),
        ``speed_grid`` and ``speed_score`` (speed grid and fit quality per
        speed), ``speed_perp`` (speed of the perpendicular fit), and
        ``blur_grid`` and ``blur_score`` (blur lengths and quality of the
        blur fit).
    """

    angle_deg: float
    speed: float
    vx: float
    vy: float
    geba: float
    quality: float
    quality_perp: float
    in_range: bool
    isotropic: bool
    threshold: float
    valid: bool
    L: int
    exposures: int
    n: int
    npad: int
    quality_blur: float = float("nan")
    blur_length: float = float("nan")
    blur_like: bool = False
    saturated: float = 0.0
    diag: Optional[dict] = field(default=None, repr=False)


_KINDS = {
    "positive": "a positive integer",
    "exposures": f"an integer of at least {MIN_EXPOSURES}",
    "roi": f"an integer of at least {MIN_SIZE} (px)",
    "unit": "a finite real number in [0, 1]",
    "duty": "a finite real number in [0, 1)",
    "flag": "true or false",
    "positive_real": "a finite real number > 0",
}
_LOWEST = {"positive": 1, "exposures": MIN_EXPOSURES, "roi": MIN_SIZE}


def _check_option(x, name: str, kind: str):
    """Return the option value ``x`` as int, float or bool, or raise ValueError.

    ``x`` must be a real, finite number, not text and not a bool (except
    for ``kind='flag'``, which takes True/False or 0/1). ``kind`` is one of
    the keys of ``_KINDS``. Same rules and messages as
    ``matlab/+fftvel/+internal/check_option.m``.
    """
    what = _KINDS[kind]
    if isinstance(x, str):
        got = f"text '{x}'"
    elif x is None:
        got = "None"
    else:
        a = np.asarray(x)
        if a.ndim != 0:
            got = f"shape {a.shape}"
        elif a.dtype.kind == "b":
            if kind == "flag":
                return bool(a)
            got = "a logical value"
        elif a.dtype.kind == "c":
            got = "a complex value"
        elif a.dtype.kind not in "iuf":
            got = f"a value of type {type(x).__name__}"
        else:
            v = float(a)
            if kind == "flag":
                ok = v in (0, 1)
            elif kind == "unit":
                ok = 0 <= v <= 1
            elif kind == "duty":
                ok = 0 <= v < 1
            elif kind == "positive_real":
                ok = bool(np.isfinite(v)) and v > 0
            else:
                ok = bool(np.isfinite(v)) and v == np.floor(v) and v >= _LOWEST[kind]
            if ok:
                if kind == "flag":
                    return bool(v)
                return int(v) if kind in _LOWEST else v
            got = f"{v:g}"
    raise ValueError(f"{name} must be {what}, got {got}")


def angle_difference(a, b):
    """Smallest difference of two directions modulo 180 degrees.

    Returns ``a - b`` in degrees, wrapped to ``[-90, 90)``:
    ``mod(a - b + 90, 180) - 90``. A direction and its opposite count as the
    same, as for ``Estimate.angle_deg``, so ``abs(d)`` is the angle between
    the two directions, at most 90 degrees; ``d > 0`` when ``a`` lies
    counter-clockwise of ``b``. Element-wise for arrays (NumPy
    broadcasting); a float for two scalars. Same as
    ``fftvel.angle_difference`` in MATLAB.

    >>> angle_difference(179, 1), angle_difference(1, 179)
    (-2.0, 2.0)
    """
    d = np.mod(np.asarray(a, dtype=float) - np.asarray(b, dtype=float) + 90, 180) - 90
    return float(d) if d.ndim == 0 else d


def hann_window(n: int) -> np.ndarray:
    """Symmetric Hann window, equal to ``numpy.hanning`` and MATLAB ``hann``."""
    if n == 1:
        return np.ones(1)
    k = np.arange(n)
    return 0.5 - 0.5 * np.cos(2 * np.pi * k / (n - 1))


def crop_square(img: np.ndarray, n: Optional[int] = None) -> np.ndarray:
    """Central ``n``-by-``n`` crop (default: the largest central square)."""
    img = np.asarray(img)
    if img.ndim != 2:
        raise ValueError(f"expected a 2-D grey-scale image, got shape {img.shape}. Convert a colour "
                         f"image to grey. For a frame stack use bp, K = pulsed_from_stack(frames, "
                         f"exposures=K), which averages the first K frames, or all of them if K is not "
                         f"given and their number is odd; then pass exposures=K to estimate_velocity")
    h, w = img.shape
    if n is None:
        n = min(h, w)
    if n > min(h, w):
        raise ValueError(f"crop of {n} px does not fit a {h}x{w} image")
    r0 = (h - n) // 2
    c0 = (w - n) // 2
    return img[r0:r0 + n, c0:c0 + n]


def pulsed_from_stack(frames, *, L: Optional[int] = None, exposures: Optional[int] = None,
                      axis: int = -1):
    """Pulsed-exposure image from a stack of single exposures (video frames).

    Takes the frames along ``axis`` (default: the last axis) and drops all
    other axes of length 1, so a MATLAB ``getdata`` stack (rows x cols x 1 x
    frames) works as it is. Returns ``(bp, K)``: the mean of the first ``K``
    frames as a 2-D float image, and ``K``. ``K`` is ``exposures`` if given,
    else ``2L+1`` if ``L`` (the half-count of Horn's thesis) is given, else
    the number of frames, which must then be odd. Pass ``K`` on as
    ``estimate_velocity(bp, exposures=K)``. All arguments after ``frames``
    are keyword-only.

    Example: ``bp, K = pulsed_from_stack(frames, exposures=9)`` averages the
    first 9 frames.
    """
    a = np.asarray(frames)
    shape = a.shape
    if a.ndim < 3:
        raise ValueError(f"expected a stack of 2-D frames (rows, cols, frames), got shape {shape}")
    a = np.moveaxis(a, axis, -1)
    a = np.squeeze(a, axis=tuple(i for i in range(a.ndim - 1) if a.shape[i] == 1))
    if a.ndim != 3:
        raise ValueError(f"expected a stack of 2-D frames (rows, cols, frames), got shape {shape}")
    nf = a.shape[2]
    if nf < 2:
        raise ValueError(f"a frame stack needs at least 2 frames, got shape {shape}")
    if L is not None:
        L = _check_option(L, "L", "positive")
    if exposures is not None:
        K = _check_option(exposures, "exposures", "exposures")
    elif L is not None:
        K = 2 * L + 1
    elif nf % 2 == 0:
        raise ValueError(f"the stack has {nf} frames, an even number; pass L or exposures "
                         f"to choose how many to average")
    else:
        K = nf
    if K > nf:
        raise ValueError(f"{K} exposures need {K} frames, the stack has {nf}")
    return a[:, :, :K].astype(float).mean(axis=2), K


def dilate(mask: np.ndarray, radius: int) -> np.ndarray:
    """``mask`` grown by a disc: True where a True pixel lies within ``radius`` px.

    Pixels count as within the disc when ``dy**2 + dx**2 <= radius**2``;
    nothing is wrapped around the borders. clip_outliers uses it so that the
    faint rim of a glint goes with its clipped core.
    """
    mask = np.asarray(mask, dtype=bool)
    out = mask.copy()
    h, w = mask.shape
    for dy in range(-radius, radius + 1):
        for dx in range(-radius, radius + 1):
            if 0 < dy * dy + dx * dx <= radius * radius:
                out[max(dy, 0):h + min(dy, 0), max(dx, 0):w + min(dx, 0)] |= \
                    mask[max(-dy, 0):h + min(-dy, 0), max(-dx, 0):w + min(-dx, 0)]
    return out


def saturated_share(img) -> float:
    """Share of the pixels of ``img`` that sit at its maximum, if that looks like saturation, else 0.

    For an integer image whose maximum is the largest value of its type
    (255 for uint8, 65535 for uint16): the share of pixels at that value,
    however small. Otherwise (floating point or logical images, and integer
    images whose maximum lies below that of their type, such as 12-bit data
    in uint16): the share of pixels at the maximum when it is at least 0.5 %
    and the image is not constant, else 0. ``img`` is used as it is, before
    any conversion; estimate_velocity passes the analysed square.
    """
    a = np.asarray(img)
    top = a.max()
    if a.dtype.kind in "iu" and top == np.iinfo(a.dtype).max:
        return np.count_nonzero(a == top) / a.size
    if top > a.min():
        share = np.count_nonzero(a == top) / a.size
        if share >= SAT_SHARE:
            return share
    return 0.0


def log_spectrum(bp: np.ndarray, pad_factor: int = DEFAULT_PAD_FACTOR, roi: Optional[int] = None,
                 clip_outliers: bool = False):
    """Centred log-magnitude spectrum of the windowed, mean-free image.

    Returns ``(S, n, npad)``. ``S = log(|B| + 0.05 * median|B|)`` of the
    ``npad``-point FFT (``npad = pad_factor * n``) with zero frequency at index
    ``npad // 2``. The image is cropped to its central ``roi``-by-``roi``
    square (default and for ``roi`` larger than the image: the largest
    central square). With ``clip_outliers``, pixels more than 6 robust
    standard deviations (``1.4826 * MAD``) from the median, and all pixels
    within 3 px of them (:func:`dilate`), are set to the median first.
    """
    pad_factor = _check_option(pad_factor, "pad_factor", "positive")
    img = np.asarray(bp)
    side = None
    if roi is not None:
        roi = _check_option(roi, "roi", "roi")
        if img.ndim == 2:
            side = min(roi, *img.shape)
    x = np.asarray(crop_square(img, side), dtype=float)     # converts only the crop; no copy for float
    n = x.shape[0]
    if n < MIN_SIZE:
        raise ValueError(f"image must be at least {MIN_SIZE}x{MIN_SIZE} pixels")
    if not np.all(np.isfinite(x)):
        raise ValueError("image contains NaN or Inf values")
    if _check_option(clip_outliers, "clip_outliers", "flag"):
        med = np.median(x)
        s = MAD_SCALE * np.median(np.abs(x - med))
        if s > 0:
            x = np.where(dilate(np.abs(x - med) > CLIP_SIGMA * s, CLIP_RADIUS), med, x)
    npad = pad_factor * n
    x = x - x.mean()
    w = hann_window(n)
    B = np.fft.fftshift(np.fft.fft2(x * np.outer(w, w), (npad, npad)))
    A = np.abs(B)
    fl = LOG_FLOOR * np.median(A)
    if fl <= 0:
        fl = LOG_FLOOR * A.mean()
    if fl <= 0:
        fl = BLANK_FLOOR     # blank image: flat spectrum, the estimate comes out not valid
    S = np.log(A + fl)
    return S, n, npad


_count_cache: dict = {}


def mean_projection(S: np.ndarray, theta_deg, rmax: int):
    """Mean of ``S`` along parallel lines at each angle (a normalised Radon).

    Returns ``(M, j, C)``: ``M[:, k]`` is the mean of ``S`` along the lines at
    signed distance ``j = -rmax..rmax`` from the centre for angle
    ``theta_deg[k]``; ``C`` holds the bin weights (those of the last call with
    more than one angle are cached). Only samples inside the disc of radius
    ``rmax`` are used, each split linearly between its two neighbouring bins.
    """
    npad = S.shape[0]
    c = npad // 2
    k = np.arange(npad) - c
    U, V = np.meshgrid(k, -k)            # U: x to the right, V: y up
    inside = U**2 + V**2 <= rmax**2
    u = U[inside].astype(float)
    v = V[inside].astype(float)
    s = S[inside]
    nb = 2 * rmax + 3                    # bins -rmax-1 .. rmax+1
    theta = np.atleast_1d(np.asarray(theta_deg, dtype=float))

    key = (npad, rmax, theta.tobytes())
    C = _count_cache.get(key)
    have_C = C is not None
    if not have_C:
        C = np.zeros((2 * rmax + 1, theta.size))
    M = np.zeros((2 * rmax + 1, theta.size))
    for a, th in enumerate(theta):
        t = th * np.pi / 180
        r = u * np.cos(t) + v * np.sin(t)
        j0 = np.floor(r)
        f = r - j0
        b = (j0 + rmax + 1).astype(np.int64)
        w = s * f
        acc = np.bincount(b, s - w, nb) + np.bincount(b + 1, w, nb)
        if not have_C:
            cnt = np.bincount(b, 1 - f, nb) + np.bincount(b + 1, f, nb)
            C[:, a] = cnt[1:-1]
        M[:, a] = acc[1:-1] / np.maximum(C[:, a], WEIGHT_FLOOR)
    if not have_C and theta.size > 1:     # keep the full angle grid, not single-angle calls
        _count_cache.clear()
        _count_cache[key] = C
    j = np.arange(-rmax, rmax + 1)
    return M, j, C


def ripple_profile(M: np.ndarray, j: np.ndarray, rfrac: float = RFRAC) -> np.ndarray:
    """Total variation ``sum(|diff(M)|)`` per column over ``|j| <= rfrac*max(j)``."""
    sel = np.abs(j) <= rfrac * j.max()
    return np.abs(np.diff(M[sel], axis=0)).sum(axis=0)


def parabolic_offset(y: np.ndarray, i: int, periodic: bool) -> float:
    """Offset of the vertex of the parabola through ``y[i-1], y[i], y[i+1]``."""
    n = y.size
    if periodic:
        ym = y[(i - 1) % n]
        yp = y[(i + 1) % n]
    else:
        if i <= 0 or i >= n - 1:
            return 0.0
        ym = y[i - 1]
        yp = y[i + 1]
    den = ym - 2 * y[i] + yp
    if den < 0:
        return float(0.5 * (ym - yp) / den)
    return 0.0


def _window_autocorr(tau: np.ndarray, n: int) -> np.ndarray:
    """Normalised autocorrelation of the ``n``-point Hann window (closed form).

    Matches the discrete autocorrelation of ``hann_window(n)`` at integer lags
    to better than 1e-7 for ``n >= 64`` (3e-5 at ``n = 16``).
    """
    x = np.abs(tau) / (n - 1)
    r = np.zeros_like(x)
    inside = x < 1
    x = x[inside]
    r[inside] = (1 - x) * (2 + np.cos(2 * np.pi * x)) / 3 + np.sin(2 * np.pi * x) / (2 * np.pi)
    return r


def dirichlet_template(jpos: np.ndarray, v: np.ndarray, L: int, npad: int, n: int,
                       angle_deg: float, exposures: Optional[int] = None,
                       duty: float = 0.0) -> np.ndarray:
    """Expected log spectrum along ``v`` for ``K`` exposures, one column per speed.

    ``K = 2L+1`` unless ``exposures`` is given. For a random texture the
    expected power of the windowed pulsed image at bin ``j`` of the
    projection along ``angle_deg`` is::

        P(j) = [K + 2*sum_{m=1}^{K-1} (K-m) * rho(m*v) * cos(2*pi*j*m*|v|/npad)] / K**2

    with ``rho(m*v) = rho_w(m*vx) * rho_w(m*vy)`` the autocorrelation of the
    separable Hann window at the offset between exposures ``m`` steps apart.
    The closed form includes the blur by the window exactly and cannot alias.
    ``duty > 0`` (exposure time / time step) adds the motion blur within each
    exposure: ``P`` is multiplied by ``sinc(j*|v|*duty/npad)**2``.
    Returns ``0.5*log(P + 1e-3)``.
    """
    jpos = np.asarray(jpos, dtype=float).reshape(-1, 1)
    v = np.asarray(v, dtype=float).reshape(1, -1)
    c = abs(np.cos(np.deg2rad(angle_deg)))
    s = abs(np.sin(np.deg2rad(angle_deg)))
    K = 2 * L + 1 if exposures is None else int(exposures)
    P = K * np.ones((jpos.size, v.size))
    for m in range(1, K):
        w = (K - m) * _window_autocorr(m * v * c, n) * _window_autocorr(m * v * s, n)
        P += 2 * np.cos(2 * np.pi * jpos * (m * v) / npad) * w
    P = np.maximum(P / K**2, 0)
    if duty > 0:
        P = P * np.sinc(jpos * v * duty / npad) ** 2
    return 0.5 * np.log(P + TEMPLATE_FLOOR)


def box_template(jpos: np.ndarray, d: np.ndarray, npad: int, n: int, angle_deg: float) -> np.ndarray:
    """Expected log spectrum along the motion for one continuous exposure, one column per blur length.

    A single exposure during which the texture moves by ``d`` px along
    ``angle_deg`` (motion blur, no flashes) multiplies the spectrum by
    ``sinc(j*d/npad)**2``. As in :func:`dirichlet_template` the window is
    included exactly: the blur is sampled every ``h = 0.5`` px, so that::

        P(j) = (h/d) * [1 + 2*sum_{q>=1} max(1 - q*h/d, 0) * rho(q*h) * cos(2*pi*j*q*h/npad)]

    with ``rho`` the autocorrelation of the window. For ``d = K*h`` this is
    the template of ``K`` exposures ``h`` apart. Returns ``0.5*log(P + 1e-3)``.
    """
    jpos = np.asarray(jpos, dtype=float).reshape(-1, 1)
    d = np.asarray(d, dtype=float).reshape(1, -1)
    c = abs(np.cos(np.deg2rad(angle_deg)))
    s = abs(np.sin(np.deg2rad(angle_deg)))
    u = BLUR_STEP * np.arange(1, int(np.ceil(d.max() / BLUR_STEP)) + 1)      # offsets 0 < u <= max(d)
    rho = _window_autocorr(u * c, n) * _window_autocorr(u * s, n)
    C = np.cos(2 * np.pi * jpos * u / npad) * rho
    W = np.maximum(1 - u.reshape(-1, 1) / d, 0)
    P = np.maximum(BLUR_STEP / d * (1 + 2 * (C @ W)), 0)
    return 0.5 * np.log(P + TEMPLATE_FLOOR)


def _template_scores(g: np.ndarray, jpos: np.ndarray, T: np.ndarray):
    """Correlation of the detrended profile ``g`` with each detrended column of ``T``.

    A cubic polynomial in ``jpos`` is removed from both; columns whose
    detrended norm is below ``1e-6`` of their norm (a flat model) score 0.
    Returns ``(score, gp, Tp)``: the scores and the detrended profile and
    templates. Shared by :func:`fit_speed` and :func:`fit_blur`.
    """
    x = jpos / jpos.max()
    X = np.column_stack([x**p for p in range(TREND_DEGREE, -1, -1)])     # cubic: x**3, x**2, x, 1
    Q, _ = np.linalg.qr(X)
    gp = g - Q @ (Q.T @ g)
    Tp = T - Q @ (Q.T @ T)
    tn = np.sqrt((Tp**2).sum(axis=0))
    score = (gp @ Tp) / (np.linalg.norm(gp) * tn + np.finfo(float).tiny)
    score[tn <= FLAT_MODEL * np.sqrt((T**2).sum(axis=0))] = 0     # flat model: no information
    return score, gp, Tp


def fit_speed(g: np.ndarray, jpos: np.ndarray, L: int, npad: int, n: int, angle_deg: float,
              exposures: Optional[int] = None, duty: float = 0.0):
    """Speed from the projection along the direction of motion ``angle_deg``.

    Correlates the detrended half profile ``g(jpos)`` with the detrended
    Dirichlet model (:func:`dirichlet_template`, ``K = 2L+1`` exposures
    unless ``exposures`` is given, motion blur ``duty``) for speeds on a
    logarithmic grid. A cubic polynomial in ``jpos`` absorbs the smooth
    texture spectrum. The measurable range runs from
    ``vmin = 1.5*npad/(K*max(jpos))`` to ``npad/4``; the grid extends down
    to ``vmin/4`` as a guard band, and ``fit['in_range']`` is False when the
    best speed lands there or on the last grid point (the top end). With
    ``pad_factor > 2`` the grid reaches speeds at which exposures one step
    apart no longer overlap within the window; the model is then flat
    (``duty = 0``), and its score is set to 0 instead of a ratio of
    rounding errors.
    Returns ``(speed, quality, fit)`` with
    ``fit = {'v', 'score', 'vmin', 'in_range', 'model'}``.
    """
    jpos = np.asarray(jpos, dtype=float).ravel()
    g = np.asarray(g, dtype=float).ravel()
    K = 2 * L + 1 if exposures is None else int(exposures)
    vmin = VMIN_FACTOR * npad / (K * jpos.max())
    vmax = npad / MIN_LINE_SPACING
    nv = int(np.floor((np.log(vmax) - np.log(vmin)) / SPEED_STEP)) + 1
    kext = int(np.floor(np.log(GUARD_FACTOR) / SPEED_STEP))
    v = np.exp(np.log(vmin) + np.arange(-kext, nv) * SPEED_STEP)

    T = dirichlet_template(jpos, v, L, npad, n, angle_deg, K, duty)
    score, gp, Tp = _template_scores(g, jpos, T)

    k = int(np.argmax(score))
    quality = float(score[k])
    dk = parabolic_offset(score, k, periodic=False)
    speed = float(np.exp(np.log(v[k]) + dk * SPEED_STEP))
    tk = Tp[:, k]
    model = (g - gp) + (gp @ tk) / (tk @ tk) * tk
    in_range = bool(kext <= k < v.size - 1)
    return speed, quality, {"v": v, "score": score, "vmin": vmin, "in_range": in_range, "model": model}


def fit_blur(g: np.ndarray, jpos: np.ndarray, npad: int, n: int, angle_deg: float):
    """Length of a continuous motion blur that explains the profile along ``angle_deg``.

    The same fit as :func:`fit_speed` with :func:`box_template` (one long
    exposure, no flashes) in place of the Dirichlet model. A blur of ``d``
    px has the main lobe of ``K`` flashes ``d/K`` apart; only the periodic
    side lines of the flashes tell them apart. The grid of blur lengths is
    ``K`` times the speed grid of :func:`fit_speed` (it does not depend on
    ``K``), from ``1.5*npad/max(jpos)/4`` up to ``npad/4``, where the zeros
    of the sinc are four bins apart. Returns ``(length, quality, fit)``
    with ``fit = {'d', 'score'}``.
    """
    jpos = np.asarray(jpos, dtype=float).ravel()
    g = np.asarray(g, dtype=float).ravel()
    dref = VMIN_FACTOR * npad / jpos.max()                  # K * vmin of fit_speed
    dmax = npad / MIN_LINE_SPACING
    nd = int(np.floor((np.log(dmax) - np.log(dref)) / SPEED_STEP)) + 1
    kext = int(np.floor(np.log(GUARD_FACTOR) / SPEED_STEP))
    d = np.exp(np.log(dref) + np.arange(-kext, nd) * SPEED_STEP)
    score, _, _ = _template_scores(g, jpos, box_template(jpos, d, npad, n, angle_deg))
    k = int(np.argmax(score))
    dk = parabolic_offset(score, k, periodic=False)
    length = float(np.exp(np.log(d[k]) + dk * SPEED_STEP))
    return length, float(score[k]), {"d": d, "score": score}


def null_quality(m: int) -> float:
    """Fit quality that images without motion stay below, for ``m`` independent profile points.

    Under the null hypothesis the Fisher transform ``atanh(q)`` of the best
    correlation scales with ``1/sqrt(m-5)`` (four points are taken by the
    cubic, one by the fit). On the 1008 calibration images without motion
    of experiments/calibrate_threshold.py (white noise, random and camera
    textures, 32 to 256 px, L = 1 to 10, ``pad_factor`` 2) the largest value
    was ``tanh(8.41/sqrt(m-5))``; on all 2268 simulated images without
    motion, held-out textures included, it was ``tanh(8.64/sqrt(m-5))``, and
    none of them was valid with 9.
    """
    return float("inf") if m <= NULL_OFFSET else float(np.tanh(NULL_K / np.sqrt(m - NULL_OFFSET)))


def estimate_velocity(bp, *, L: int = DEFAULT_L, exposures: Optional[int] = None,
                      duty: float = 0.0, roi: Optional[int] = None,
                      clip_outliers: bool = False, pad_factor: int = DEFAULT_PAD_FACTOR,
                      min_quality: float = DEFAULT_MIN_QUALITY, reject_blur: bool = False,
                      diagnostics: bool = False) -> Estimate:
    """Estimate direction and speed of the motion in a pulsed-exposure image.

    Parameters
    ----------
    bp : 2-D array
        Sum or mean of ``K`` exposures taken at equal time steps while the
        texture moves by ``v`` per step: ``K`` strobe flashes within one
        camera frame (a multi-flash image), or the mean of ``K`` consecutive
        video frames (see :func:`pulsed_from_stack`). Non-square images are
        cropped to the central square, never resized.
    L : int
        Half-count of Horn's thesis, ``K = 2L+1``, used when ``exposures``
        is not given; default 4 (9 exposures).
    exposures : int, optional
        Number of exposures ``K >= 2``, even or odd: the flashes in one
        multi-flash image, or the averaged frames. Overrides ``L``; default
        ``2L+1``. For slow motion a wrong count scales the speed by
        ``K_true/K`` without lowering the quality.
    duty : float
        Duration of one exposure (flash pulse, or camera exposure time)
        divided by the time step, in ``[0, 1)``; default 0 (short flashes).
        The motion blur within each exposure, ``duty*|v|`` pixels, must
        either be given here or stay well below about 3 px.
    roi : int, optional
        Side in pixels of the central square that is analysed (at least 16);
        default and values larger than the image: the largest central square.
        The run time grows with ``roi**2``.
    clip_outliers : bool
        Set pixels more than 6 robust standard deviations (``1.4826*MAD``)
        from the median, and all pixels within 3 px of them, to the median
        before the FFT, against static bright spots such as glints or hot
        pixels. Default False.
    pad_factor : int
        Zero padding factor of the FFT; default 2.
    min_quality : float
        Lowest fit quality accepted for ``valid``, in ``[0, 1]``; default
        0.6. Small images need more (see ``Estimate.threshold``).
    reject_blur : bool
        Also require ``Estimate.blur_like`` to be False for ``valid``,
        against single exposures of a moving surface (motion blur, no
        flashes). Most multi-flash images of slow motion (below about
        2.2 px per step) with noise are blur-like too, so this rejects most
        valid results there. Default False.
    diagnostics : bool
        Keep intermediate results in ``Estimate.diag`` (needed for plotting).
        Default False.

    All arguments after ``bp`` are keyword-only, e.g.
    ``estimate_velocity(bp, exposures=9)``. All options are checked; invalid
    values raise ValueError. :func:`to_physical` converts the result to
    length units per second.

    Examples
    --------
    >>> from fftvel import simulate_pulsed_image
    >>> bp = simulate_pulsed_image(256, (5, -3), 4, rng=0)    # L = 4, so K = 2L+1 = 9 flashes
    >>> est = estimate_velocity(bp, exposures=9)
    >>> print(f"{est.angle_deg:.1f} deg, {est.speed:.2f} px per flash, valid {est.valid}")
    31.0 deg, 5.83 px per flash, valid True
    """
    L = _check_option(L, "L", "positive")
    K = 2 * L + 1 if exposures is None else _check_option(exposures, "exposures", "exposures")
    duty = _check_option(duty, "duty", "duty")
    if roi is not None:
        roi = _check_option(roi, "roi", "roi")
    pad_factor = _check_option(pad_factor, "pad_factor", "positive")
    min_quality = _check_option(min_quality, "min_quality", "unit")
    clip_outliers = _check_option(clip_outliers, "clip_outliers", "flag")
    reject_blur = _check_option(reject_blur, "reject_blur", "flag")
    diagnostics = _check_option(diagnostics, "diagnostics", "flag")
    S, n, npad = log_spectrum(bp, pad_factor, roi, clip_outliers)
    saturated = saturated_share(crop_square(np.asarray(bp), n))     # the analysed square as given
    rmax = npad // 2 - 1

    theta = np.arange(0, 180, ANGLE_STEP)
    M, j, _ = mean_projection(S, theta, rmax)
    wr = ripple_profile(M, j, RFRAC)
    i = int(np.argmax(wr))
    angle = float((theta[i] + ANGLE_STEP * parabolic_offset(wr, i, periodic=True)) % 180)
    geba = float(wr[i] / wr.mean()) if wr.mean() > 0 else 1.0   # 1: no ripple at all

    m, _, _ = mean_projection(S, [angle], rmax)
    m = m[:, 0]
    pos = (j >= PROFILE_START) & (j <= RFRAC * rmax)
    g = 0.5 * (m[pos] + m[::-1][pos])
    jpos = j[pos]
    speed, quality, fit = fit_speed(g, jpos, L, npad, n, angle, K, duty)

    # Isotropy check: the same fit perpendicular to the direction found
    ap = (angle + 90.0) % 180.0
    mp, _, _ = mean_projection(S, [ap], rmax)
    mp = mp[:, 0]
    speed_perp, quality_perp, _ = fit_speed(0.5 * (mp[pos] + mp[::-1][pos]), jpos, L, npad, n,
                                            ap, K, duty)
    isotropic = bool(quality - quality_perp < ISO_GAP
                     and abs(np.log(speed_perp / speed)) < ISO_LOG_SPEED)

    # Blur check: one continuous exposure (motion blur, no flashes) along the same direction. Only the
    # repeated lines of the flashes, npad/speed bins apart, tell the two apart; below npad/max(jpos)
    # they lie beyond the profile
    blur_length, quality_blur, fit_b = fit_blur(g, jpos, npad, n, angle)
    blur_like = bool(quality - quality_blur < BLUR_GAP and speed * jpos.max() < npad)

    a = np.deg2rad(angle)
    npts = jpos.size
    m_eff = min(npts, (INDEPENDENT_PAD * n * npts) // npad)    # independent points: padding interpolates
    threshold = max(min_quality, null_quality(m_eff))
    valid = bool(quality >= threshold and fit["in_range"] and not isotropic
                 and not (reject_blur and blur_like))
    est = Estimate(
        angle_deg=angle,
        speed=speed,
        vx=float(speed * np.cos(a)),
        vy=float(-speed * np.sin(a)),
        geba=geba,
        quality=quality,
        quality_perp=float(quality_perp),
        in_range=fit["in_range"],
        isotropic=isotropic,
        threshold=threshold,
        valid=valid,
        L=L,
        exposures=K,
        n=n,
        npad=npad,
        quality_blur=quality_blur,
        blur_length=blur_length,
        blur_like=blur_like,
        saturated=float(saturated),
    )
    if diagnostics:
        est.diag = {
            "spectrum": S,
            "theta": theta,
            "wr": wr,
            "profile_j": jpos,
            "profile": g,
            "model": fit["model"],
            "speed_grid": fit["v"],
            "speed_score": fit["score"],
            "speed_perp": speed_perp,
            "blur_grid": fit_b["d"],
            "blur_score": fit_b["score"],
        }
    return est


@dataclass(frozen=True)
class PhysicalVelocity:
    """Result of :func:`to_physical`: the velocity in length units per second.

    Attributes
    ----------
    speed : ``|v|`` in the length unit of ``gsd`` per second (mm/s for
        ``gsd`` in mm per px).
    vx, vy : velocity in image coordinates (x along columns, y along rows
        pointing down) in the same unit; ``-(vx, vy)`` is equally possible.
    angle_deg : direction as in ``Estimate.angle_deg``, unchanged.
    """

    speed: float
    vx: float
    vy: float
    angle_deg: float


def to_physical(est, gsd: float, period: float) -> PhysicalVelocity:
    """Velocity in length units per second from an estimate in pixels per step.

    ``gsd`` is the ground sample distance, the length on the moving surface
    that one pixel covers (for example mm per px, measured from an image of
    a ruler or a checkerboard in the plane of the surface). ``period`` is the
    time step in seconds: the time between two flashes for a multi-flash
    image, or between two frames for averaged frames. Then::

        speed = est.speed * gsd / period

    and ``vx``, ``vy`` likewise, in the length unit of ``gsd`` per second.
    ``est`` is an :class:`Estimate` or any object with the attributes
    ``speed``, ``vx``, ``vy`` and ``angle_deg``. It is converted whether or
    not it is valid, so check ``est.valid`` first. ``gsd`` and ``period``
    must be finite and greater than 0 (ValueError otherwise). Same as
    ``fftvel.to_physical`` in MATLAB.

    >>> e = Estimate(90.0, 2.5, 0.0, -2.5, 3.0, 0.99, 0.5, True, False, 0.7, True, 4, 9, 256, 512)
    >>> to_physical(e, gsd=0.04, period=0.001).speed      # 2.5 px * 0.04 mm/px / 1 ms
    100.0
    """
    gsd = _check_option(gsd, "gsd", "positive_real")
    period = _check_option(period, "period", "positive_real")
    missing = [k for k in ("speed", "vx", "vy", "angle_deg") if not hasattr(est, k)]
    if missing:
        raise TypeError(f"est must be a result of estimate_velocity (an Estimate); "
                        f"{type(est).__name__} has no {', '.join(missing)}")
    return PhysicalVelocity(
        speed=float(est.speed * gsd / period),
        vx=float(est.vx * gsd / period),
        vy=float(est.vy * gsd / period),
        angle_deg=float(est.angle_deg),
    )
