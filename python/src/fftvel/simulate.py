"""Synthetic pulsed-exposure images with known velocity."""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from .core import crop_square

__all__ = ["make_texture", "simulate_pulsed_image"]


def make_texture(m: int, corr_len: float = 1.5, rng=None) -> np.ndarray:
    """Random periodic ``m``-by-``m`` texture with autocorrelation close to ``exp(-r/corr_len)``.

    White noise is shaped with the matching power spectrum
    ``(1 + (2*pi*corr_len*f)**2)**-1.5``; the result has zero mean and unit
    standard deviation.
    """
    rng = np.random.default_rng(rng)
    f = np.fft.fftfreq(m)
    FX, FY = np.meshgrid(f, f)
    psd = (1 + (2 * np.pi * corr_len) ** 2 * (FX**2 + FY**2)) ** -1.5
    t = np.real(np.fft.ifft2(np.fft.fft2(rng.standard_normal((m, m))) * np.sqrt(psd)))
    return (t - t.mean()) / t.std(ddof=1)


def simulate_pulsed_image(n: int, v: Sequence[float], L: int, corr_len: float = 1.5,
                          texture: Optional[np.ndarray] = None, noise: float = 0.0,
                          rng=None, return_frames: bool = False):
    """Mean of ``2L+1`` exposures of a texture moving by ``v`` per step under a fixed window.

    ``v = (vx, vy)`` in pixels per step, x along columns, y along rows
    (pointing down). Shifts are exact band-limited (Fourier) shifts, so ``v``
    may be fractional; the texture is large enough that no wrapped content
    reaches the ``n``-by-``n`` window and there are no zero-filled borders.

    ``texture`` replaces the random texture (e.g. a real camera frame); it must
    be at least ``n + 2*L*max|v| + 32`` pixels in both directions. ``noise`` is
    the standard deviation of white noise added to the result, relative to
    the standard deviation of the noise-free image.

    Returns ``bp`` or ``(bp, frames)`` with ``frames`` of shape ``(n, n, 2L+1)``.
    """
    rng = np.random.default_rng(rng)
    vx, vy = float(v[0]), float(v[1])
    margin = 2 * L * max(abs(vx), abs(vy)) + 32
    if texture is None:
        m = int(2 ** np.ceil(np.log2(n + margin)))
        tex = make_texture(m, corr_len, rng)
    else:
        tex = np.asarray(texture, dtype=float)
        if min(tex.shape) < n + margin:
            raise ValueError(f"texture too small for n={n} and |v|={max(abs(vx), abs(vy))} with L={L}")
        side = min(tex.shape)
        tex = crop_square(tex, side - side % 2)
        tex = tex - tex.mean()
    m = tex.shape[0]
    f = np.fft.fftfreq(m)
    FX, FY = np.meshgrid(f, f)
    T = np.fft.fft2(tex)
    r0 = (m - n) // 2
    win = slice(r0, r0 + n)

    frames = np.zeros((n, n, 2 * L + 1)) if return_frames else None
    D = np.zeros((m, m), dtype=complex)
    for k in range(-L, L + 1):
        e = np.exp(-2j * np.pi * (FX * vx + FY * vy) * k)
        D += e
        if return_frames:
            frames[:, :, k + L] = np.real(np.fft.ifft2(T * e))[win, win]
    big = np.real(np.fft.ifft2(T * D / (2 * L + 1)))
    bp = big[win, win]
    if noise > 0:
        bp = bp + noise * bp.std(ddof=1) * rng.standard_normal((n, n))
    return (bp, frames) if return_frames else bp
