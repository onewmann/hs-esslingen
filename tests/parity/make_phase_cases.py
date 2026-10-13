"""Write the parity cases of the phase correlation.

    python tests/parity/make_phase_cases.py

creates ``phase_cases.mat`` with image pairs ``a01, b01, ...`` and the
results of the Python implementation: ``expected`` holds one row
``dx, dy, peak`` per pair, ``shift`` the true shift ``dx, dy`` of the
content from ``a`` to ``b`` (NaN for the blank pair, which has none).
``tests/matlab/test_phase_correlation.m`` and
``tests/python/test_phase_correlation.py`` check that MATLAB/Octave and
Python both reproduce ``expected`` from the pairs to 1e-9. The pairs are
crops of data/surface_texture.png moved by exact band-limited (Fourier)
shifts of the whole texture, so ``b`` has no wrapped or zero-filled content.
This file is separate from make_cases.py and cases.mat, which hold the
cases of the estimator.
"""

from pathlib import Path

import matplotlib.image as mpimg
import numpy as np
from scipy.io import savemat

from fftvel import phase_correlation
from fftvel.simulate import make_texture

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

# (rows, cols, dx, dy, noise, kind); noise relative to the standard deviation of each crop
CASES = [
    (64, 64, 2.3, -1.6, 0.0, "float"),
    (96, 80, -4.55, 3.25, 0.0, "float"),       # non-square: dx along the 80 columns
    (65, 63, 0.5, -0.5, 0.0, "float"),         # odd sizes, half-pixel shift
    (96, 96, 7.8, 0.1, 0.25, "float"),         # noise
    (64, 64, 0.0, 0.0, 0.0, "float"),          # identical images: (0, 0) and peak 1
    (64, 64, 3.7, 2.2, 0.0, "uint8"),          # 8-bit frames, as from a camera
    (48, 48, 11.4, -10.6, 0.0, "random"),      # random texture, shift near a quarter of the square
    (64, 64, -0.25, 0.75, 1.0, "float"),       # strong noise: a broad peak
    (1, 96, -3.35, 0.0, 0.0, "float"),         # a single row
    (16, 16, np.nan, np.nan, 0.0, "blank"),    # constant images: (0, 0) and peak 0
]


def shifted(tex, dx, dy):
    """``tex`` moved by ``(dx, dy)`` px (x along columns, y down), periodic band-limited shift."""
    fy = np.fft.fftfreq(tex.shape[0])
    fx = np.fft.fftfreq(tex.shape[1])
    return np.real(np.fft.ifft2(np.fft.fft2(tex) * np.exp(-2j * np.pi * np.add.outer(fy * dy, fx * dx))))


def main():
    rng = np.random.default_rng(20261011)
    surface = mpimg.imread(REPO / "data" / "surface_texture.png").astype(float)[:512, :512]
    mats = {}
    expected = []
    for k, (h, w, dx, dy, noise, kind) in enumerate(CASES, start=1):
        tex = make_texture(256, 1.5, rng) if kind == "random" else surface
        r0, c0 = rng.integers(32, tex.shape[0] - max(h, w) - 32, 2)
        a = tex[r0:r0 + h, c0:c0 + w]
        if kind == "blank":
            a, b = np.ones((h, w)), np.full((h, w), 7.0)
        elif dx or dy:
            b = shifted(tex, dx, dy)[r0:r0 + h, c0:c0 + w]
        else:
            b = a.copy()
        if noise:
            a = a + noise * a.std() * rng.standard_normal(a.shape)
            b = b + noise * b.std() * rng.standard_normal(b.shape)
        if kind == "uint8":
            lo, hi = tex.min(), tex.max()
            a, b = (np.round(255 * (x - lo) / (hi - lo)).astype(np.uint8) for x in (a, b))
        else:
            a, b = a.astype(np.float32), b.astype(np.float32)
        mats[f"a{k:02d}"] = a
        mats[f"b{k:02d}"] = b
        expected.append(phase_correlation(a, b))
        gx, gy, peak = expected[-1]
        print(f"case {k}: {h}x{w} {kind}, true ({dx}, {dy}), got ({gx:.4f}, {gy:.4f}), peak {peak:.3f}")
    mats["expected"] = np.array(expected)
    mats["shift"] = np.array([[c[2], c[3]] for c in CASES], dtype=float)
    savemat(HERE / "phase_cases.mat", mats, do_compression=True)


if __name__ == "__main__":
    import argparse

    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()
    main()
