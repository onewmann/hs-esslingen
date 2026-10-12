"""Model options (exposures, duty), ROI, outlier clipping, isotropy check, blur check, saturation
and frame stacks.

The same checks are in tests/matlab/test_model_options.m.
"""

from dataclasses import replace

import numpy as np
import pytest
from scipy.io import loadmat
from scipy.signal import fftconvolve

from fftvel import angle_difference, estimate_velocity, pulsed_from_stack, simulate_pulsed_image
from fftvel.core import (BLUR_STEP, CLIP_RADIUS, box_template, crop_square, dilate, dirichlet_template,
                         log_spectrum, saturated_share)


@pytest.fixture(scope="module")
def texture(repo):
    import matplotlib.image as mpimg

    return mpimg.imread(repo / "data" / "surface_texture.png").astype(float)


def velocity(speed, angle):
    a = np.deg2rad(angle)
    return speed * np.cos(a), -speed * np.sin(a)


def assert_correct(e, speed, angle, tol=0.01):
    assert abs(angle_difference(e.angle_deg, angle)) < 0.5, e
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


def test_clip_outliers_takes_the_rim_along():
    # the clipped pixels and all pixels within CLIP_RADIUS of them are set to the median
    _, x = hot_pixel_image(0)
    y = x.copy()
    rows, cols = np.mgrid[0:96, 0:96]
    for r, c in ((30, 30), (50, 60)):
        y[(rows - r) ** 2 + (cols - c) ** 2 <= CLIP_RADIUS**2] = np.median(x)
    assert np.array_equal(log_spectrum(x, clip_outliers=True)[0], log_spectrum(y)[0])


def test_box_template_is_the_dirichlet_template_of_close_exposures():
    j = np.arange(2, 230)
    for K in (2, 9, 40):
        a = box_template(j, [K * BLUR_STEP], 512, 256, 30.0)
        b = dirichlet_template(j, [BLUR_STEP], 4, 512, 256, 30.0, exposures=K)
        assert np.abs(a - b).max() < 1e-12, K


def blurred(n, length, angle, seed):
    """One continuous exposure: the texture moves by ``length`` px (10 exposures of duty 0.99, no gaps),
    with 10 % noise as a camera adds."""
    return simulate_pulsed_image(n, velocity(length / 10, angle), 4, exposures=10, duty=0.99, noise=0.1,
                                 rng=seed)


@pytest.mark.parametrize("length, angle", [(5.0, 30.0), (9.0, 100.0)])
def test_single_exposure_is_blur_like(length, angle):
    # a blur of d px looks like K flashes d/K apart; the blur check flags it
    bp = blurred(128, length, angle, 3)
    e = estimate_velocity(bp)
    assert e.valid and abs(e.speed * e.exposures / length - 1) < 0.05
    assert e.blur_like and abs(e.blur_length / length - 1) < 0.05
    r = estimate_velocity(bp, reject_blur=True)
    assert not r.valid and r == replace(e, valid=False)


def test_fast_flashes_are_not_blur_like():
    bp = simulate_pulsed_image(128, (7.0, 3.0), 4, rng=1)
    e = estimate_velocity(bp)
    assert not e.blur_like and e.quality - e.quality_blur > 0.3
    assert estimate_velocity(bp, reject_blur=True) == e


def test_dilate():
    m = np.zeros((7, 9), dtype=bool)
    m[3, 4] = True
    d = dilate(m, 2)
    assert d.sum() == 13 and d[3, 2] and d[1, 4] and d[2, 3] and not d[1, 3]     # a disc of radius 2
    m[3, 4], m[0, 0] = False, True
    assert dilate(m, 2).sum() == 6                        # a quarter disc at the corner, nothing wraps
    assert np.array_equal(dilate(m, 0), m)
    assert not dilate(np.zeros((5, 5), dtype=bool), 3).any()


def test_saturated_share():
    # same rules as tests/matlab/test_model_options.m
    x = np.full((20, 20), 7, dtype=np.uint8)
    x[0, :2] = 255
    assert saturated_share(x) == 2 / 400                  # at the uint8 maximum: any share counts
    x[0, :2] = 254
    assert saturated_share(x) == 2 / 400                  # below it: counts from 0.5 %
    x[0, 1] = 7
    assert saturated_share(x) == 0                        # 0.25 %
    assert saturated_share(x.astype(float)) == 0
    x[0, 1] = 254
    assert saturated_share(x.astype(np.float32)) == 2 / 400
    big = np.full((100, 100), 50, dtype=np.uint16)
    big[0, 0] = 65535
    assert saturated_share(big) == 1e-4
    assert saturated_share(big.astype(float)) == 0        # one pixel at a float maximum is no evidence
    assert saturated_share(np.full((8, 8), 255, dtype=np.uint8)) == 1
    assert saturated_share(np.full((8, 8), 255.0)) == 0   # constant: nothing to compare with
    assert saturated_share(np.full((8, 8), 7, dtype=np.int16)) == 0
    i16 = np.zeros((10, 10), dtype=np.int16)
    i16[0, 0] = 32767
    assert saturated_share(i16) == 0.01
    b = np.zeros((10, 10), dtype=bool)
    b[0, 0] = True
    assert saturated_share(b) == 0.01


def test_saturated_counts_the_analysed_square_before_clipping():
    bp = simulate_pulsed_image(160, (3.0, -1.0), 4, rng=5)
    x = np.clip(np.round(100 + 15 * bp / bp.std()), 0, 255).astype(np.uint8)
    x[:16, :16] = 255                                     # outside the central 128 px
    assert estimate_velocity(x, roi=128).saturated == 0
    assert estimate_velocity(x).saturated == 256 / 160**2
    x[60:70, 60:70] = 255                                 # inside
    e = estimate_velocity(x, roi=128, clip_outliers=True)
    assert e.saturated == 100 / 128**2 and estimate_velocity(x, roi=128).saturated == e.saturated
    assert estimate_velocity(x.astype(float), roi=128).saturated == 100 / 128**2   # float: >= 0.5 %


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
