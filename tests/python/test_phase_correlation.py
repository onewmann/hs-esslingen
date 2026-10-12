"""fftvel.phase_correlation: conventions, sub-pixel accuracy and parity.

tests/parity/phase_cases.mat holds image pairs and the Python results
(written by tests/parity/make_phase_cases.py); tests/matlab/test_phase_correlation.m
checks MATLAB/Octave against the same numbers.
"""

import matplotlib.image as mpimg
import numpy as np
import pytest
from scipy.io import loadmat

from fftvel import phase_correlation
from fftvel.reference import UPSAMPLE, _upsampled_peak
from fftvel.simulate import make_texture


def shifted(tex, dx, dy):
    """``tex`` moved by ``(dx, dy)`` px (x along columns, y down), periodic band-limited shift."""
    fy = np.fft.fftfreq(tex.shape[0])
    fx = np.fft.fftfreq(tex.shape[1])
    return np.real(np.fft.ifft2(np.fft.fft2(tex) * np.exp(-2j * np.pi * np.add.outer(fy * dy, fx * dx))))


@pytest.fixture(scope="module")
def texture():
    # periodic, so shifts wrap cleanly; tests/parity/make_cases.py stores the same texture for the
    # MATLAB tests (texture256 in fixed_inputs.mat)
    return make_texture(256, 1.5, np.random.default_rng(11))


def test_parity_cases(repo):
    d = loadmat(repo / "tests" / "parity" / "phase_cases.mat")
    expected, shift = d["expected"], d["shift"]
    assert expected.shape == (10, 3) and shift.shape == (10, 2)
    for k in range(expected.shape[0]):
        a, b = d[f"a{k + 1:02d}"], d[f"b{k + 1:02d}"]
        got = phase_correlation(a, b)
        assert got == pytest.approx(expected[k], abs=1e-9), k + 1
        if k + 1 != 8 and not np.isnan(shift[k]).any():    # not the case with strong noise or the blank one
            assert got[:2] == pytest.approx(shift[k], abs=0.06), k + 1
    assert tuple(expected[9]) == (0, 0, 0)                  # blank images


def test_sign_and_axis_order(texture):
    # content of a moves by +5.3 columns (right) and -2.6 rows (up) in b
    a = texture[40:104, 20:116]                                     # 64 rows x 96 columns
    b = shifted(texture, 5.3, -2.6)[40:104, 20:116]
    dx, dy, peak = phase_correlation(a, b)
    assert (dx, dy) == pytest.approx((5.3, -2.6), abs=0.02)
    assert 0.5 < peak <= 1
    assert phase_correlation(b, a)[:2] == pytest.approx((-5.3, 2.6), abs=0.02)        # swapped: opposite sign
    assert phase_correlation(a.T, b.T)[:2] == pytest.approx((-2.6, 5.3), abs=0.02)    # transposed: axes swap


def test_result_is_wrapped_to_half_the_size():
    # dx in (-w/2, w/2], dy in (-h/2, h/2]: a periodic shift by 60 of 64 columns is -4
    a = make_texture(64, 1.5, np.random.default_rng(12))
    dx, dy, _ = phase_correlation(a, np.roll(a, (-5, 60), axis=(0, 1)))
    assert (dx, dy) == pytest.approx((-4, -5), abs=0.02)


def test_shifts_near_half_the_size(texture):
    # sub-pixel shifts just inside +-w/2 and +-h/2 keep their sign through the refinement and the
    # wrapping; whole periodic frames, so only the window limits the overlap
    a = texture[::4, ::2]                                           # 64 rows x 128 columns, periodic
    for d in ((63.6, 0.3), (-63.6, -0.3), (0.3, 31.6), (-0.3, -31.6)):
        assert phase_correlation(a, shifted(a, *d))[:2] == pytest.approx(d, abs=0.05), d


def test_identical_images():
    a = np.random.default_rng(3).standard_normal((48, 40))
    dx, dy, peak = phase_correlation(a, a)
    assert abs(dx) < 1e-12 and abs(dy) < 1e-12 and peak == pytest.approx(1, abs=1e-12)
    # integer input is converted as a whole: uint8 frames give the same result as their float copies
    u = np.round(100 + 20 * a).clip(0, 255).astype(np.uint8)
    v = np.roll(u, 3, axis=1)
    assert phase_correlation(u, v) == phase_correlation(u.astype(float), v.astype(float))


def test_blank_images_and_single_rows(texture):
    # nothing to correlate: no shift and peak 0, as with the parabolic fit used before
    assert phase_correlation(np.ones((32, 32)), np.full((32, 32), 7.0)) == (0.0, 0.0, 0.0)
    # NaN in an image: no shift and a NaN peak, not a shift to the corner of the fine grid
    a = texture[:32, :32].copy()
    a[5, 7] = np.nan
    with np.errstate(invalid="ignore"):
        dx, dy, peak = phase_correlation(a, texture[:32, :32])
    assert (dx, dy) == (0, 0) and np.isnan(peak)
    # a single row or column: the shift along it, exactly 0 across it
    a = texture[:1, 40:168]
    b = shifted(texture[:1, :], -5.7, 0)[:, 40:168]
    dx, dy, _ = phase_correlation(a, b)
    assert dx == pytest.approx(-5.7, abs=0.03) and dy == 0
    dx2, dy2, _ = phase_correlation(a.T, b.T)
    assert dx2 == 0 and dy2 == pytest.approx(dx, abs=1e-12)


def test_sub_pixel_accuracy(repo):
    # Fourier-shifted crops of the camera texture: the parabolic fit used before was off by up
    # to 0.17 px (0.12 px median) between pixels; the interpolated peak by about 0.01 px at most
    tex = mpimg.imread(repo / "data" / "surface_texture.png").astype(float)[:384, :384]
    rng = np.random.default_rng(5)
    err = []
    for frac in np.linspace(0, 1, 11)[:-1]:
        d = rng.integers(-6, 7, 2) + np.array([frac, 1 - frac])
        a = tex[128:256, 128:256]
        b = shifted(tex, *d)[128:256, 128:256]
        err.append(np.array(phase_correlation(a, b)[:2]) - d)
    err = np.array(err)
    assert np.abs(err).max() < 0.02
    assert np.abs(err.mean(axis=0)).max() < 0.005


def test_upsampled_peak_of_a_pure_shift():
    # the lag map of exp(-2 pi i k.t / m) is a Dirichlet kernel centred on t (off the 1/20 px grid)
    assert UPSAMPLE == 20
    for h, w in ((63, 65), (64, 64), (96, 128)):
        ky = np.fft.fftfreq(h) * h
        kx = np.fft.fftfreq(w) * w
        for ty, tx in ((3.317, -7.4523), (-0.4871, 0.2611), (0.0731, -0.9612)):
            X = np.exp(-2j * np.pi * np.add.outer(ky * ty / h, kx * tx / w))
            i, j = round(ty), round(tx)
            oy, ox, peak = _upsampled_peak(X, i, j)
            assert (i + oy, j + ox) == pytest.approx((ty, tx), abs=1e-4)
            assert 0.998 < peak <= 1 + 1e-12
