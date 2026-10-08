"""The estimator of the original study project, ported to Python for comparison.

Processing chain of ``finaltest1.m`` (live Basler GUI, lines 102-139) as it
was before the rewrite: top-hat background removal, Wiener filter, Hann
window, zero padding to 2N, ``log(1+|FFT|)``, MATLAB-style Radon transform,
``beta = argmax WR`` and ``|v| = 2N / mean(spacing of findpeaks maxima)``
with an absolute ``MinPeakProminence`` of 0.05.

``radon_matlab`` and ``findpeaks_matlab`` reproduce MATLAB's ``radon``
(4 sub-pixels, linear binning, y up) and ``findpeaks`` semantics; the port was
checked against the original code run in GNU Octave.

Only used by ``run_benchmark.py`` to show what changed. Do not use it for
measurements.
"""

import numpy as np
from scipy import ndimage, signal

THETA = np.arange(180)


def _disk(r):
    y, x = np.mgrid[-r:r + 1, -r:r + 1]
    return (x * x + y * y) <= r * r


def _imopen(img, r):
    fp = _disk(r)
    e = ndimage.grey_erosion(img, footprint=fp, mode="constant", cval=np.inf)
    return ndimage.grey_dilation(e, footprint=fp, mode="constant", cval=-np.inf)


def _wiener2(g, n=5):
    lm = ndimage.uniform_filter(g, n, mode="constant")
    lv = ndimage.uniform_filter(g * g, n, mode="constant") - lm ** 2
    noise = lv.mean()
    gg = np.maximum(lv - noise, 0)
    lv = np.maximum(lv, noise)
    return (g - lm) / lv * gg + lm


def radon_matlab(img, theta):
    """MATLAB ``radon``: 4 sub-pixels per pixel, linear binning, y axis up."""
    M, N = img.shape
    xo = max(0, (N - 1) // 2)
    yo = max(0, (M - 1) // 2)
    rlast = int(np.ceil(np.hypot(M - 1 - yo, N - 1 - xo))) + 1
    rfirst = -rlast
    rsize = rlast - rfirst + 1
    x = np.arange(N) - xo
    y = yo - np.arange(M)
    R = np.zeros((rsize, len(theta)))
    nz = img != 0
    vals = img[nz] * 0.25
    yy, xx = np.nonzero(nz)
    X, Y = x[xx], y[yy]
    for k, th in enumerate(np.deg2rad(theta)):
        c, s = np.cos(th), np.sin(th)
        acc = np.zeros(rsize + 2)
        for dx in (-0.25, 0.25):
            for dy in (-0.25, 0.25):
                r = (X + dx) * c + (Y + dy) * s - rfirst
                ri = np.floor(r).astype(int)
                d = r - ri
                acc += np.bincount(ri, vals * (1 - d), minlength=rsize + 2)[:rsize + 2]
                acc += np.bincount(ri + 1, vals * d, minlength=rsize + 2)[:rsize + 2]
        R[:, k] = acc[:rsize]
    return R, np.arange(rfirst, rlast + 1)


def findpeaks_matlab(y, x, minprom=0.05, mindist=3):
    """MATLAB ``findpeaks``: local maxima, then prominence, then distance (tallest first)."""
    idx, _ = signal.find_peaks(y)
    if idx.size == 0:
        return np.array([])
    prom = signal.peak_prominences(y, idx)[0]
    idx = idx[prom >= minprom]
    if idx.size == 0:
        return np.array([])
    order = np.argsort(-y[idx], kind="stable")
    keep = np.ones(idx.size, bool)
    xs = x[idx]
    for i in order:
        if keep[i]:
            close = np.abs(xs - xs[i]) < mindist
            close[i] = False
            keep[close & keep] = False
    return x[np.sort(idx[keep])]


def original_estimate(bp):
    """Return ``(beta_deg, v_mag)`` exactly as the original code computed them.

    ``bp`` must be square and scaled to [0, 1] like the contrast-stretched
    frames of the original GUI.
    """
    n = bp.shape[0]
    pad = 2 * n
    bg = _imopen(bp, round(n / 20))
    bp_dn = _wiener2(np.maximum(bp - bg, 0), 5)
    h = np.hanning(n)
    p = (pad - n) // 2
    B = np.fft.fftshift(np.fft.fft2(np.pad(bp_dn * np.outer(h, h), p)))
    R, xp = radon_matlab(np.log1p(np.abs(B)), THETA)
    wr = np.abs(np.diff(R, axis=0)).sum(axis=0)
    idx = int(np.argmax(wr))
    proj = R[:, idx]
    mask = xp >= 0
    xp_p = xp[mask].astype(float)
    proj_p = proj[mask]
    locs = findpeaks_matlab(proj_p, xp_p).astype(float)
    sub = locs.copy()
    for k, loc in enumerate(locs):
        ix = int(np.argmin(np.abs(xp_p - loc)))
        if 0 < ix < len(proj_p) - 1:
            y1, y0, y2 = proj_p[ix - 1], proj_p[ix], proj_p[ix + 1]
            sub[k] = loc + (y1 - y2) / (2 * (y1 - 2 * y0 + y2))
    d = np.diff(np.sort(sub))
    return float(THETA[idx]), float(pad / d.mean()) if d.size else float("nan")
