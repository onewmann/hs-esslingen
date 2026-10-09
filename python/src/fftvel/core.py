"""Velocity estimation from a single pulsed-exposure image.

The model, the ripple profile WR and the quality measure GEBA follow Chapter 3
of J. Horn's doctoral thesis (Universitaet Karlsruhe, 2007; full reference in
README.md). The speed fit and the validity threshold were added in this
implementation.

This module mirrors the MATLAB/Octave functions in ``matlab/`` line by line,
so that both give the same numbers for the same input (see ``tests/parity``).
Only NumPy is needed.

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

__all__ = [
    "Estimate",
    "estimate_velocity",
    "log_spectrum",
    "mean_projection",
    "ripple_profile",
    "fit_speed",
    "null_quality",
    "dirichlet_template",
    "parabolic_offset",
    "hann_window",
    "crop_square",
]

RFRAC = 0.9          # fraction of the spectrum radius used for WR and the fit
SPEED_STEP = 0.005   # step of the logarithmic speed grid


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
    in_range : False if the best speed lies below the measurable range.
    threshold : quality needed for this image size, the larger of
        ``min_quality`` and ``tanh(9/sqrt(m-5))`` for ``m`` profile points.
    valid : ``quality >= threshold`` and ``in_range``.
    """

    angle_deg: float
    speed: float
    vx: float
    vy: float
    geba: float
    quality: float
    in_range: bool
    threshold: float
    valid: bool
    L: int
    n: int
    npad: int
    diag: Optional[dict] = field(default=None, repr=False)


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
        raise ValueError("expected a 2-D grey-scale image")
    h, w = img.shape
    if n is None:
        n = min(h, w)
    if n > min(h, w):
        raise ValueError(f"crop of {n} px does not fit a {h}x{w} image")
    r0 = (h - n) // 2
    c0 = (w - n) // 2
    return img[r0:r0 + n, c0:c0 + n]


def log_spectrum(bp: np.ndarray, pad_factor: int = 2):
    """Centred log-magnitude spectrum of the windowed, mean-free image.

    Returns ``(S, n, npad)``. ``S = log(|B| + 0.05 * median|B|)`` of the
    ``npad``-point FFT (``npad = pad_factor * n``) with zero frequency at index
    ``npad // 2``.
    """
    if pad_factor < 1 or int(pad_factor) != pad_factor:
        raise ValueError("pad_factor must be a positive integer")
    pad_factor = int(pad_factor)
    x = crop_square(np.asarray(bp, dtype=float))
    n = x.shape[0]
    if n < 16:
        raise ValueError("image must be at least 16x16 pixels")
    if not np.all(np.isfinite(x)):
        raise ValueError("image contains NaN or Inf values")
    npad = pad_factor * n
    x = x - x.mean()
    w = hann_window(n)
    B = np.fft.fftshift(np.fft.fft2(x * np.outer(w, w), (npad, npad)))
    A = np.abs(B)
    fl = 0.05 * np.median(A)
    if fl <= 0:
        fl = 0.05 * A.mean()
    if fl <= 0:
        fl = 1.0         # blank image: flat spectrum, the estimate comes out not valid
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
        M[:, a] = acc[1:-1] / np.maximum(C[:, a], 1e-12)
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
                       angle_deg: float) -> np.ndarray:
    """Expected log spectrum along ``v`` for ``2L+1`` exposures, one column per speed.

    For a random texture the expected power of the windowed pulsed image at
    bin ``j`` of the projection along ``angle_deg`` is::

        P(j) = [(2L+1) + 2*sum_{m=1}^{2L} (2L+1-m) * rho(m*v) * cos(2*pi*j*m*|v|/npad)] / (2L+1)**2

    with ``rho(m*v) = rho_w(m*vx) * rho_w(m*vy)`` the autocorrelation of the
    separable Hann window at the offset between exposures ``m`` steps apart.
    The closed form includes the blur by the window exactly and cannot alias.
    Returns ``0.5*log(P + 1e-3)``.
    """
    jpos = np.asarray(jpos, dtype=float).reshape(-1, 1)
    v = np.asarray(v, dtype=float).reshape(1, -1)
    c = abs(np.cos(np.deg2rad(angle_deg)))
    s = abs(np.sin(np.deg2rad(angle_deg)))
    K = 2 * L + 1
    P = K * np.ones((jpos.size, v.size))
    for m in range(1, 2 * L + 1):
        w = (K - m) * _window_autocorr(m * v * c, n) * _window_autocorr(m * v * s, n)
        P += 2 * np.cos(2 * np.pi * jpos * (m * v) / npad) * w
    P = np.maximum(P / K**2, 0)
    return 0.5 * np.log(P + 1e-3)


def fit_speed(g: np.ndarray, jpos: np.ndarray, L: int, npad: int, n: int, angle_deg: float):
    """Speed from the projection along the direction of motion ``angle_deg``.

    Correlates the detrended half profile ``g(jpos)`` with the detrended
    Dirichlet model (:func:`dirichlet_template`) for speeds on a logarithmic
    grid. A cubic polynomial in ``jpos`` absorbs the smooth texture spectrum.
    The measurable range runs from ``vmin = 1.5*npad/((2L+1)*max(jpos))`` to
    ``npad/4``; the grid extends down to ``vmin/4`` as a guard band, and
    ``fit['in_range']`` is False when the best speed lands there.
    Returns ``(speed, quality, fit)`` with
    ``fit = {'v', 'score', 'vmin', 'in_range', 'model'}``.
    """
    jpos = np.asarray(jpos, dtype=float).ravel()
    g = np.asarray(g, dtype=float).ravel()
    vmin = 1.5 * npad / ((2 * L + 1) * jpos.max())
    vmax = npad / 4
    nv = int(np.floor((np.log(vmax) - np.log(vmin)) / SPEED_STEP)) + 1
    kext = int(np.floor(np.log(4) / SPEED_STEP))
    v = np.exp(np.log(vmin) + np.arange(-kext, nv) * SPEED_STEP)

    T = dirichlet_template(jpos, v, L, npad, n, angle_deg)
    x = jpos / jpos.max()
    X = np.column_stack([x**3, x**2, x, np.ones_like(x)])
    Q, _ = np.linalg.qr(X)
    gp = g - Q @ (Q.T @ g)
    Tp = T - Q @ (Q.T @ T)
    tn = np.sqrt((Tp**2).sum(axis=0))
    score = (gp @ Tp) / (np.linalg.norm(gp) * tn + np.finfo(float).tiny)

    k = int(np.argmax(score))
    quality = float(score[k])
    dk = parabolic_offset(score, k, periodic=False)
    speed = float(np.exp(np.log(v[k]) + dk * SPEED_STEP))
    tk = Tp[:, k]
    model = (g - gp) + (gp @ tk) / (tk @ tk) * tk
    in_range = bool(kext <= k < v.size - 1)
    return speed, quality, {"v": v, "score": score, "vmin": vmin, "in_range": in_range, "model": model}


def null_quality(m: int) -> float:
    """Fit quality that images without motion stay below, for a profile of ``m`` points.

    Under the null hypothesis the Fisher transform ``atanh(q)`` of the best
    correlation scales with ``1/sqrt(m-5)`` (four points are taken by the
    cubic, one by the fit). On 1008 images without motion (white noise,
    random and camera textures, 32 to 256 px, L = 1 to 10) the largest value
    was ``tanh(8.6/sqrt(m-5))``; 9 adds a margin.
    """
    return float("inf") if m <= 5 else float(np.tanh(9 / np.sqrt(m - 5)))


def estimate_velocity(bp, L: int = 4, pad_factor: int = 2, min_quality: float = 0.6,
                      diagnostics: bool = False) -> Estimate:
    """Estimate direction and speed of the motion in a pulsed-exposure image.

    Parameters
    ----------
    bp : 2-D array
        Sum or mean of ``2L+1`` exposures taken at equal time steps while the
        texture moves by ``v`` per step. Non-square images are cropped to the
        central square, never resized.
    L : int
        Half the number of exposures.
    pad_factor : int
        Zero padding factor of the FFT.
    min_quality : float
        Lowest fit quality accepted for ``valid``; small images need more
        (see ``Estimate.threshold``).
    diagnostics : bool
        Keep intermediate results in ``Estimate.diag`` (needed for plotting).
    """
    if np.ndim(L) != 0 or L < 1 or int(L) != L:
        raise ValueError("L must be a positive integer")
    L = int(L)
    S, n, npad = log_spectrum(bp, pad_factor)
    rmax = npad // 2 - 1

    theta = np.arange(180, dtype=float)
    M, j, _ = mean_projection(S, theta, rmax)
    wr = ripple_profile(M, j, RFRAC)
    i = int(np.argmax(wr))
    angle = float((theta[i] + parabolic_offset(wr, i, periodic=True)) % 180)
    geba = float(wr[i] / wr.mean()) if wr.mean() > 0 else 1.0   # 1: no ripple at all

    m, _, _ = mean_projection(S, [angle], rmax)
    m = m[:, 0]
    pos = (j >= 2) & (j <= RFRAC * rmax)
    g = 0.5 * (m[pos] + m[::-1][pos])
    jpos = j[pos]
    speed, quality, fit = fit_speed(g, jpos, L, npad, n, angle)

    a = np.deg2rad(angle)
    threshold = max(min_quality, null_quality(jpos.size))
    est = Estimate(
        angle_deg=angle,
        speed=speed,
        vx=float(speed * np.cos(a)),
        vy=float(-speed * np.sin(a)),
        geba=geba,
        quality=quality,
        in_range=fit["in_range"],
        threshold=threshold,
        valid=bool(quality >= threshold and fit["in_range"]),
        L=L,
        n=n,
        npad=npad,
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
        }
    return est
