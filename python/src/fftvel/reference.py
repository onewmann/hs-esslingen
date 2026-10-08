"""Independent reference: frame-to-frame shift by phase correlation."""

from __future__ import annotations

import numpy as np

from .core import hann_window, parabolic_offset

__all__ = ["phase_correlation"]


def phase_correlation(a: np.ndarray, b: np.ndarray):
    """Shift ``(dx, dy, peak)`` that moves the content of ``a`` onto ``b``.

    ``dx`` along columns, ``dy`` along rows (pointing down), with parabolic
    sub-pixel refinement. Both images are Hann windowed.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    h, w = a.shape
    win = np.outer(hann_window(h), hann_window(w))
    A = np.fft.fft2((a - a.mean()) * win)
    B = np.fft.fft2((b - b.mean()) * win)
    X = B * np.conj(A)
    r = np.real(np.fft.ifft2(X / np.maximum(np.abs(X), np.finfo(float).eps)))
    i, j = np.unravel_index(int(np.argmax(r)), r.shape)
    dy = i + parabolic_offset(r[:, j], i, periodic=True)
    dx = j + parabolic_offset(r[i, :], j, periodic=True)
    if dy > h / 2:
        dy -= h
    if dx > w / 2:
        dx -= w
    return float(dx), float(dy), float(r[i, j])
