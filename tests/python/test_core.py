import numpy as np
import pytest
from scipy.io import loadmat

from fftvel import (crop_square, estimate_velocity, hann_window, log_spectrum, mean_projection,
                    simulate_pulsed_image)
from fftvel.core import parabolic_offset

from conftest import angle_error


def fixed_inputs(repo):
    """Images without motion saved by tests/parity/make_cases.py, shared with the MATLAB tests."""
    return loadmat(repo / "tests" / "parity" / "fixed_inputs.mat")


@pytest.mark.parametrize("v, expected", [
    ((2, 0), 0.0),
    ((0, -2), 90.0),
    ((2, -2), 45.0),
    ((2, 2), 135.0),
    ((-3, 1), 18.434948822922),       # -v = (3, -1): up and to the right
])
def test_angle_convention(v, expected):
    bp = simulate_pulsed_image(128, v, 4, rng=11)
    e = estimate_velocity(bp)
    assert abs(angle_error(e.angle_deg, expected)) < 0.5
    assert 0 <= e.angle_deg < 180
    vv = np.array([e.vx, e.vy])
    assert min(np.linalg.norm(vv - v), np.linalg.norm(vv + v)) < 0.03 * np.hypot(*v) + 0.02


@pytest.mark.parametrize("speed", [0.7, 3, 12, 30])
@pytest.mark.parametrize("angle", [15, 80, 125])
def test_known_velocity(speed, angle):
    a = np.deg2rad(angle)
    bp = simulate_pulsed_image(256, (speed * np.cos(a), -speed * np.sin(a)), 4, rng=42)
    e = estimate_velocity(bp, L=4)
    assert abs(angle_error(e.angle_deg, angle)) < 0.5
    assert abs(e.speed / speed - 1) < 0.01
    assert e.valid


def test_more_exposures_resolve_slow_motion():
    bp = simulate_pulsed_image(256, (0.3, 0), 20, rng=1)
    e = estimate_velocity(bp, L=20)
    assert abs(e.speed / 0.3 - 1) < 0.02


def test_real_texture(repo):
    import matplotlib.image as mpimg

    tex = mpimg.imread(repo / "data" / "surface_texture.png").astype(float)
    bp = simulate_pulsed_image(256, (5, -3), 4, texture=tex, rng=0)
    e = estimate_velocity(bp)
    assert abs(angle_error(e.angle_deg, np.rad2deg(np.arctan2(3, 5)))) < 0.5
    assert abs(e.speed / np.hypot(5, 3) - 1) < 0.01


@pytest.mark.parametrize("L, speed", [(1, 0.9), (1, 1.0), (3, 0.45), (4, 0.3)])
def test_motion_below_the_range_is_not_valid(L, speed):
    a = np.deg2rad(37)
    bp = simulate_pulsed_image(256, (speed * np.cos(a), -speed * np.sin(a)), L, rng=5)
    e = estimate_velocity(bp, L=L)
    assert not e.in_range and not e.valid


@pytest.mark.parametrize("seed", [1, 2, 3])
@pytest.mark.parametrize("n", [64, 128, 256])
def test_static_scene_is_not_valid(seed, n):
    e = estimate_velocity(simulate_pulsed_image(n, (0, 0), 4, rng=seed))
    assert not e.valid


def test_no_motion_is_not_valid_at_32_px(repo):
    # about 1 in 1000 images of 32 px without motion passes, so fixed images are used
    # (the same ones as in tests/matlab/test_known_velocity.m)
    d = fixed_inputs(repo)
    for k in (1, 2, 3):
        for name in (f"static32_{k}", f"noise32_{k}"):
            e = estimate_velocity(d[name].astype(float))
            assert not e.valid and e.quality < e.threshold - 0.05, name


def test_white_noise_is_not_valid_on_small_images():
    rng = np.random.default_rng(64)
    assert not any(estimate_velocity(rng.standard_normal((64, 64))).valid for _ in range(10))


BAD_OPTIONS = [
    ("L", 0), ("L", 2.5), ("L", "4"), ("L", True), ("L", np.inf), ("L", np.nan), ("L", [4, 4]),
    ("L", None), ("L", 4 + 1j),
    ("pad_factor", 1.5), ("pad_factor", "2"), ("pad_factor", 0), ("pad_factor", False),
    ("min_quality", np.nan), ("min_quality", None), ("min_quality", "0.8"), ("min_quality", 1.5),
    ("min_quality", -0.1), ("min_quality", np.inf),
    ("exposures", 1), ("exposures", 8.5), ("exposures", "8"), ("exposures", True),
    ("duty", 1), ("duty", -0.1), ("duty", np.nan), ("duty", "0.2"), ("duty", None),
    ("roi", 8), ("roi", 100.5), ("roi", "128"), ("roi", np.inf),
    ("clip_outliers", 2), ("clip_outliers", "yes"), ("clip_outliers", None),
    ("diagnostics", "no"), ("diagnostics", 2),
]


@pytest.mark.parametrize("name, value", BAD_OPTIONS)
def test_invalid_options_are_rejected(name, value):
    # same rules as tests/matlab/test_options.m; the message names the option
    bp = np.random.default_rng(0).standard_normal((32, 32))
    with pytest.raises(ValueError, match=name):
        estimate_velocity(bp, **{name: value})


@pytest.mark.parametrize("args, name", [(("2",), "pad_factor"), ((2, "64"), "roi"),
                                        ((2, None, 2), "clip_outliers")])
def test_log_spectrum_checks_its_options(args, name):
    with pytest.raises(ValueError, match=name):
        log_spectrum(np.zeros((32, 32)), *args)


def test_options_and_errors():
    bp = simulate_pulsed_image(64, (2, 1), 4, rng=3)
    with pytest.raises(ValueError):
        estimate_velocity(np.zeros((8, 8)))
    with pytest.raises(ValueError, match=r"\(32, 32, 3\).*pulsed_from_stack"):
        estimate_velocity(np.zeros((32, 32, 3)))
    nan_img = bp.copy()
    nan_img[3, 4] = np.nan
    with pytest.raises(ValueError):
        estimate_velocity(nan_img)
    e1 = estimate_velocity(bp, L=np.int32(4), pad_factor=np.uint8(2), min_quality=np.float32(0.5),
                           exposures=np.int64(9), duty=np.float32(0), roi=np.uint16(64),
                           clip_outliers=np.bool_(False))
    e2 = estimate_velocity(bp, L=4.0, min_quality=0.5)
    assert e1 == e2 and e1.exposures == 9 and e1.L == 4
    e = estimate_velocity(bp, min_quality=0.99)
    assert e.valid == (e.quality >= max(0.99, e.threshold) and e.in_range and not e.isotropic)
    assert e.diag is None
    e = estimate_velocity(bp, diagnostics=True)
    assert e.diag["spectrum"].shape == (128, 128)
    blank = estimate_velocity(np.ones((64, 64)))
    assert not blank.valid and blank.quality == 0 and blank.geba == 1
    u8 = np.round(255 * (bp - bp.min()) / np.ptp(bp)).astype(np.uint8)
    assert abs(angle_error(estimate_velocity(u8).angle_deg, estimate_velocity(bp).angle_deg)) < 1


def test_helpers():
    assert np.allclose(hann_window(5), [0, 0.5, 1, 0.5, 0], atol=1e-15)
    assert np.allclose(hann_window(64), np.hanning(64))
    img = np.arange(24).reshape(6, 4).T          # 4 x 6
    assert np.array_equal(crop_square(img), img[:, 1:5])
    assert np.array_equal(crop_square(img.T), img.T[1:5, :])
    M, j, _ = mean_projection(3 * np.ones((64, 64)), [0, 33, 90], 31)
    assert np.allclose(M, 3)
    assert np.array_equal(j, np.arange(-31, 32))
    y = -((np.arange(9) - 3.3) ** 2)
    i = int(np.argmax(y))
    assert abs(i + parabolic_offset(y, i, periodic=False) - 3.3) < 1e-12
    assert parabolic_offset(y, 0, periodic=False) == 0
