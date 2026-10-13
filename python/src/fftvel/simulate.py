"""Synthetic pulsed-exposure images with known velocity."""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from .core import DEFAULT_CORR_LEN, SIM_MARGIN, _check_option, crop_square
from .io import to_grey

__all__ = ["make_texture", "simulate_pulsed_image"]


def make_texture(m: int, corr_len: float = DEFAULT_CORR_LEN, rng=None) -> np.ndarray:
    """Random periodic ``m``-by-``m`` texture with autocorrelation close to ``exp(-r/corr_len)``.

    White noise is shaped with the matching power spectrum
    ``(1 + (2*pi*corr_len*f)**2)**-1.5``; the result has zero mean and unit
    standard deviation. ``rng`` is a seed or a ``numpy.random.Generator``
    (default: a new generator with a random seed).
    """
    rng = np.random.default_rng(rng)
    f = np.fft.fftfreq(m)
    FX, FY = np.meshgrid(f, f)
    psd = (1 + (2 * np.pi * corr_len) ** 2 * (FX**2 + FY**2)) ** -1.5
    t = np.real(np.fft.ifft2(np.fft.fft2(rng.standard_normal((m, m))) * np.sqrt(psd)))
    return (t - t.mean()) / t.std(ddof=1)


def _velocity(v) -> tuple:
    """``(vx, vy)`` as floats; ValueError unless ``v`` holds exactly two finite real numbers.

    Same rule and message as ``matlab/+fftvel/simulate_pulsed_image.m``.
    """
    try:
        a = None if isinstance(v, (str, bytes)) else np.asarray(v)
    except ValueError:                          # e.g. a ragged list
        a = None
    if a is None or a.size != 2 or a.dtype.kind not in "iuf" or not np.all(np.isfinite(a)):
        raise ValueError("v must hold two finite real numbers (vx, vy)")
    return float(a.flat[0]), float(a.flat[1])


def simulate_pulsed_image(n: int, v: Sequence[float], L: int, *,
                          exposures: Optional[int] = None, duty: float = 0.0,
                          texture: Optional[np.ndarray] = None, corr_len: float = DEFAULT_CORR_LEN,
                          noise: float = 0.0, rng=None, return_frames: bool = False):
    """Mean of ``2L+1`` exposures of a texture moving by ``v`` per step under a fixed window.

    ``L`` is the half-count of Horn's thesis: ``L = 4`` gives ``K = 2L+1 = 9``
    flashes. ``v = (vx, vy)`` in pixels per step, x along columns, y along
    rows (pointing down). Shifts are exact band-limited (Fourier) shifts, so
    ``v`` may be fractional; the texture is large enough that no wrapped
    content reaches the ``n``-by-``n`` window and there are no zero-filled
    borders.

    ``texture`` replaces the random texture (e.g. a real camera frame); it must
    be at least ``n + (K-1+duty)*max|v| + 32`` pixels in both directions. An
    RGB(A) texture is converted to grey with the weights 0.2989, 0.5870,
    0.1140. ``noise`` is the standard deviation of white noise added to the
    result, relative to the standard deviation of the noise-free image.

    ``exposures`` sets the number of exposures ``K`` and overrides ``L``
    (default ``2L+1``; an even ``K`` puts the exposures at half-integer steps
    around the centre). ``duty`` is the duration of each exposure divided by the time step,
    in ``[0, 1)``: the texture also moves during each exposure (box blur of
    ``duty*|v|`` pixels along ``v``). Default 0: instantaneous exposures.

    ``corr_len`` is the correlation length of the random texture in px
    (default 1.5, see :func:`make_texture`); ``rng`` seeds it and the noise.
    All arguments after ``L`` are keyword-only.

    Returns ``bp`` or ``(bp, frames)`` with ``frames`` of shape ``(n, n, K)``.

    Example::

        bp, frames = simulate_pulsed_image(256, (5, -3), 4, return_frames=True)   # L = 4, so K = 9
        est = estimate_velocity(bp, exposures=9, diagnostics=True)
        plot_estimate(bp, est, frames[:, :, 4])     # the middle exposure, k = 0
    """
    rng = np.random.default_rng(rng)
    L = _check_option(L, "L", "positive")
    K = 2 * L + 1 if exposures is None else _check_option(exposures, "exposures", "exposures")
    duty = _check_option(duty, "duty", "duty")
    vx, vy = _velocity(v)
    margin = ((K - 1) + duty) * max(abs(vx), abs(vy)) + SIM_MARGIN
    if texture is None:
        m = int(2 ** np.ceil(np.log2(n + margin)))
        tex = make_texture(m, corr_len, rng)
    else:
        try:
            tex = to_grey(texture)                        # RGB(A) -> grey, MATLAB rgb2gray weights
        except ValueError:
            raise ValueError(f"texture must be a 2-D grey-scale or an RGB(A) image, got shape "
                             f"{np.shape(texture)}") from None
        if min(tex.shape) < n + margin:
            raise ValueError(f"texture too small for n={n} and |v|={max(abs(vx), abs(vy))} "
                             f"with {K} exposures")
        side = min(tex.shape)
        tex = crop_square(tex, side - side % 2)
        tex = tex - tex.mean()
    m = tex.shape[0]
    f = np.fft.fftfreq(m)
    FX, FY = np.meshgrid(f, f)
    T = np.fft.fft2(tex)
    if duty > 0:                         # motion during each exposure: box blur along v
        T = T * np.sinc(duty * (FX * vx + FY * vy))
    r0 = (m - n) // 2
    win = slice(r0, r0 + n)

    frames = np.zeros((n, n, K)) if return_frames else None
    D = np.zeros((m, m), dtype=complex)
    for q in range(K):
        k = q - (K - 1) / 2              # -L..L for K = 2L+1
        e = np.exp(-2j * np.pi * (FX * vx + FY * vy) * k)
        D += e
        if return_frames:
            frames[:, :, q] = np.real(np.fft.ifft2(T * e))[win, win]
    big = np.real(np.fft.ifft2(T * D / K))
    bp = big[win, win]
    if noise > 0:
        bp = bp + noise * bp.std(ddof=1) * rng.standard_normal((n, n))
    return (bp, frames) if return_frames else bp
