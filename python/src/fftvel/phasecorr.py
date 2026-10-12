"""Shift between two frames by phase correlation, a check that is independent of the estimator."""

from __future__ import annotations

import numpy as np

from .core import hann_window, parabolic_offset

__all__ = ["phase_correlation"]

# Upsampling factor of the sub-pixel refinement: the lag map is evaluated on a
# grid of 1/UPSAMPLE px within 1 px of its largest sample. The same number is
# used in matlab/+fftvel/phase_correlation.m.
UPSAMPLE = 20


def phase_correlation(a: np.ndarray, b: np.ndarray):
    """Shift ``(dx, dy, peak)`` that moves the content of ``a`` onto ``b``.

    ``dx`` along columns, ``dy`` along rows (pointing down), in pixels, each
    in ``(-w/2, w/2]`` and ``(-h/2, h/2]``. Both images are made mean-free and
    Hann windowed; the lag map is the inverse FFT of their normalised
    cross-power spectrum, and its largest sample gives the shift to the
    nearest pixel.

    Sub-pixel refinement: the lag map is a trigonometric polynomial in the
    lag, so it can be evaluated between the samples directly from the
    cross-power spectrum (band-limited interpolation). It is evaluated on a
    grid of 1/UPSAMPLE px within 1 px of the largest sample, as two small
    matrix products with complex exponentials instead of a large zero-padded
    FFT, and the largest value on that grid is refined by a parabola along
    each axis. A parabola through three samples of the full-resolution lag
    map instead, whose peak is only about one pixel wide, is off by up to
    about 0.17 px for shifts between pixels. On Fourier-shifted 256 px crops
    of data/surface_texture.png without noise the refined shift is off by
    0.002 px (median; at most 0.005 px); the remaining error comes from
    content moving under the fixed window.

    ``peak`` is the height of the interpolated peak, the largest value on
    the fine grid (1 for two identical images, lower when content moves in
    or out of the window, changes or is noisy). Blank images give
    ``(0, 0, 0)``; NaN or Inf in either image gives ``(0, 0, nan)``.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    h, w = a.shape
    win = np.outer(hann_window(h), hann_window(w))
    A = np.fft.fft2((a - a.mean()) * win)
    B = np.fft.fft2((b - b.mean()) * win)
    X = B * np.conj(A)
    X = X / np.maximum(np.abs(X), np.finfo(float).eps)
    r = np.real(np.fft.ifft2(X))
    i, j = np.unravel_index(int(np.argmax(r)), r.shape)
    oy, ox, peak = _upsampled_peak(X, i, j)
    dy = i + oy
    dx = j + ox
    if dy > h / 2:
        dy -= h
    if dx > w / 2:
        dx -= w
    return float(dx), float(dy), float(peak)


def _upsampled_peak(X, i, j):
    """Offsets ``(oy, ox)`` of the maximum of the lag map ``real(ifft2(X))`` from the sample ``(i, j)``.

    Also returns the largest value on the grid. Row ``p`` and column ``q`` of
    ``R`` hold the lag map at ``(i + uy[p], j + ux[q])``.
    """
    h, w = X.shape
    u = np.arange(-UPSAMPLE, UPSAMPLE + 1) / UPSAMPLE
    uy = u if h > 1 else np.zeros(1)            # an axis of length 1 has no sub-pixel lag
    ux = u if w > 1 else np.zeros(1)
    ky = np.fft.fftfreq(h) * h                  # signed frequency indices, as in the FFT
    kx = np.fft.fftfreq(w) * w
    ey = np.exp(2j * np.pi * np.outer(i + uy, ky) / h)
    ex = np.exp(2j * np.pi * np.outer(kx, j + ux) / w)
    R = np.real(ey @ X @ ex) / (h * w)
    p, q = np.unravel_index(int(np.argmax(R)), R.shape)
    cy, cx = uy.size // 2, ux.size // 2
    if not R[p, q] > R[cy, cx]:                 # nothing higher than the sample itself: blank images, NaN
        p, q = cy, cx
    oy = uy[p] + parabolic_offset(R[:, q], p, periodic=False) / UPSAMPLE
    ox = ux[q] + parabolic_offset(R[p, :], q, periodic=False) / UPSAMPLE
    return oy, ox, R[p, q]
