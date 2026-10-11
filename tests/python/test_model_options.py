"""Model options (exposures, duty), ROI, outlier clipping, isotropy check and frame stacks.

The same checks are in tests/matlab/test_model_options.m.
"""

import numpy as np
import pytest
from scipy.io import loadmat
from scipy.signal import fftconvolve

from fftvel import crop_square, estimate_velocity, pulsed_from_stack, simulate_pulsed_image

from conftest import angle_error


@pytest.fixture(scope="module")
def texture(repo):
    import matplotlib.image as mpimg

    return mpimg.imread(repo / "data" / "surface_texture.png").astype(float)


def velocity(speed, angle):
    a = np.deg2rad(angle)
    return speed * np.cos(a), -speed * np.sin(a)


def assert_correct(e, speed, angle, tol=0.01):
    assert abs(angle_error(e.angle_deg, angle)) < 0.5, e
    assert abs(e.speed / speed - 1) < tol, e
    assert e.valid, e


@pytest.mark.parametrize("K", [8, 10])
@pytest.mark.parametrize("speed, angle", [(0.6, 25), (5.0, 110)])
def test_even_number_of_exposures(K, speed, angle):
    bp = simulate_pulsed_image(256, velocity(speed, angle), 4, exposures=K, rng=K)
    e = estimate_velocity(bp, exposures=K)
    assert e.exposures == K and e.L == 4
    assert_correct(e, speed, angle)
    # exposures wins over L for the model
    assert estimate_velocity(bp, L=1, exposures=K).speed == e.speed


def test_simulated_exposures_are_the_mean_of_the_frames():
    bp, fr = simulate_pulsed_image(64, (2.5, -1.0), 4, exposures=6, duty=0.3, rng=2, return_frames=True)
    assert fr.shape == (64, 64, 6)
    assert np.allclose(bp, fr.mean(axis=2), atol=1e-12)
    bp9 = simulate_pulsed_image(64, (2.5, -1.0), 4, rng=2)
    assert np.array_equal(bp9, simulate_pulsed_image(64, (2.5, -1.0), 4, exposures=9, rng=2))


@pytest.mark.parametrize("duty", [0.2, 0.5])
@pytest.mark.parametrize("speed, angle", [(34.6, 10.1), (12.0, 120.0)])
def test_duty_on_the_real_texture(texture, duty, speed, angle):
    bp = simulate_pulsed_image(256, velocity(speed, angle), 4, texture=texture, duty=duty)
    assert_correct(estimate_velocity(bp, duty=duty), speed, angle)


def test_duty_zero_changes_nothing():
    bp = simulate_pulsed_image(128, (3, -2), 4, rng=1)
    assert estimate_velocity(bp, duty=0) == estimate_velocity(bp)


def test_roi():
    bp = simulate_pulsed_image(300, (3.0, -1.0), 4, rng=4)[:, :280]
    e = estimate_velocity(bp, roi=128)
    assert e.n == 128 and e.npad == 256
    assert_correct(e, np.hypot(3, 1), np.rad2deg(np.arctan2(1, 3)))
    assert e == estimate_velocity(crop_square(bp, 128))
    assert estimate_velocity(bp, roi=1000) == estimate_velocity(bp)       # larger: largest square


def hot_pixel_image(seed):
    """8-bit texture (mean 100, std 15 DN) moving by (1.5, -0.5) with two static hot pixels."""
    bp = simulate_pulsed_image(96, (1.5, -0.5), 4, rng=seed)
    x = np.clip(np.round(100 + 15 * bp / bp.std()), 0, 255)
    clean = x.copy()
    x[30, 30] = 255
    x[50, 60] = 255
    return clean, x


@pytest.mark.parametrize("seed", [0, 1])
def test_clip_outliers(seed):
    clean, x = hot_pixel_image(seed)
    off = estimate_velocity(x)
    assert off.valid and abs(off.speed / np.hypot(20, 30) - 1) < 0.01     # reported as motion
    assert_correct(estimate_velocity(x, clip_outliers=True), np.hypot(1.5, 0.5),
                   np.rad2deg(np.arctan2(0.5, 1.5)))
    assert estimate_velocity(clean, clip_outliers=True) == estimate_velocity(clean)


def test_defocused_static_texture_is_not_valid(texture):
    r = 6
    y, x = np.mgrid[-r:r + 1, -r:r + 1]
    disk = (x**2 + y**2 <= 2.5**2).astype(float)
    img = fftconvolve(texture[300:300 + 256 + 2 * r, 300:300 + 256 + 2 * r], disk / disk.sum(), mode="valid")
    e = estimate_velocity(img)
    assert e.quality >= e.threshold and e.in_range                    # passed before the check
    assert e.isotropic and not e.valid


def test_repeating_texture_is_not_valid(repo):
    img = loadmat(repo / "tests" / "parity" / "fixed_inputs.mat")["tile16"].astype(float)
    e = estimate_velocity(img)
    assert e.quality >= e.threshold and e.in_range and abs(e.speed - 16) < 0.1
    assert e.isotropic and not e.valid


def test_motion_is_not_isotropic():
    e = estimate_velocity(simulate_pulsed_image(128, (3, -2), 4, rng=1))
    assert not e.isotropic and e.quality - e.quality_perp > 0.15


@pytest.mark.parametrize("pad", [6, 8])
@pytest.mark.parametrize("n", [32, 64])
def test_padding_does_not_lower_the_threshold(repo, n, pad):
    d = loadmat(repo / "tests" / "parity" / "fixed_inputs.mat")
    rng = np.random.default_rng(100 * n + pad)
    images = [d[f"noise{n}_{k}"].astype(float) for k in (1, 2, 3)]
    images += [rng.standard_normal((n, n)) for _ in range(5)]
    e2 = estimate_velocity(images[0], pad_factor=2)
    for img in images:
        e = estimate_velocity(img, pad_factor=pad)
        assert not e.valid
        assert e.threshold > e2.threshold - 0.02              # about the same independent samples as pad 2


def test_flat_model_scores_zero(repo):
    # with pad_factor > 2 the speed grid reaches speeds at which exposures one step apart do not
    # overlap within the window; the model is flat there and must not score rounding errors
    d = loadmat(repo / "tests" / "parity" / "fixed_inputs.mat")
    e = estimate_velocity(d["noise32_1"].astype(float), pad_factor=8, diagnostics=True)
    a = np.deg2rad(e.angle_deg)
    flat = e.diag["speed_grid"] * max(abs(np.cos(a)), abs(np.sin(a))) >= e.n - 1
    assert flat.any() and np.all(e.diag["speed_score"][flat] == 0)


def test_pulsed_from_stack():
    _, fr = simulate_pulsed_image(64, (3.0, -2.0), 5, rng=3, return_frames=True)     # 11 frames
    u8 = np.round(255 * (fr - fr.min()) / np.ptp(fr)).astype(np.uint8)
    bp, K = pulsed_from_stack(u8)
    assert K == 11 and bp.dtype == float and bp.shape == (64, 64)
    assert np.array_equal(bp, u8.astype(float).mean(axis=2))
    camera = u8[:, :, None, :10]                       # rows x cols x 1 x frames, as from getdata
    bp, K = pulsed_from_stack(camera, L=4)
    assert K == 9 and np.array_equal(bp, u8[:, :, :9].astype(float).mean(axis=2))
    bp, K = pulsed_from_stack(camera, L=2, exposures=8)            # exposures wins
    assert K == 8 and np.array_equal(bp, u8[:, :, :8].astype(float).mean(axis=2))
    first = np.moveaxis(u8, 2, 0)                                  # frames first
    assert np.array_equal(pulsed_from_stack(first, axis=0)[0], pulsed_from_stack(u8)[0])
    with pytest.raises(ValueError, match="even"):
        pulsed_from_stack(camera)
    with pytest.raises(ValueError, match="11 frames"):
        pulsed_from_stack(camera, L=5)
    with pytest.raises(ValueError):
        pulsed_from_stack(u8[:, :, 0])                             # one 2-D frame
    with pytest.raises(ValueError):
        pulsed_from_stack(u8[:, :, :1])                            # a single frame
    with pytest.raises(ValueError, match="L"):
        pulsed_from_stack(u8, L="2")
    with pytest.raises(ValueError, match="exposures"):
        pulsed_from_stack(u8, exposures=1)


def test_simulate_with_rgb_texture(texture):
    rgb = np.stack([texture, 0.5 * texture, texture**2], axis=2)
    grey = 0.2989 * rgb[:, :, 0] + 0.5870 * rgb[:, :, 1] + 0.1140 * rgb[:, :, 2]
    a = simulate_pulsed_image(64, (2, 1), 2, texture=rgb)
    assert np.array_equal(a, simulate_pulsed_image(64, (2, 1), 2, texture=grey))
    rgba = np.concatenate([rgb, np.ones_like(texture)[:, :, None]], axis=2)
    assert np.array_equal(a, simulate_pulsed_image(64, (2, 1), 2, texture=rgba))
    with pytest.raises(ValueError, match="2-D"):
        simulate_pulsed_image(64, (2, 1), 2, texture=rgb[:, :, :2])


@pytest.mark.parametrize("kwargs", [{"L": "2"}, {"L": 0}, {"exposures": 1}, {"exposures": 6.5},
                                    {"duty": 1}, {"duty": "0.2"}])
def test_simulate_checks_its_options(kwargs):
    # same checks as tests/matlab/test_model_options.m
    args = {"L": 2, **kwargs}
    with pytest.raises(ValueError, match=next(iter(kwargs))):
        simulate_pulsed_image(64, (2, 1), args.pop("L"), rng=0, **args)
