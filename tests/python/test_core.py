import numpy as np
import pytest

from fftvel import (crop_square, estimate_velocity, hann_window, mean_projection,
                    simulate_pulsed_image)
from fftvel.core import parabolic_offset

from conftest import angle_error


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


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_static_scene_is_not_valid(seed):
    e = estimate_velocity(simulate_pulsed_image(256, (0, 0), 4, rng=seed))
    assert not e.valid


def test_options_and_errors():
    bp = simulate_pulsed_image(64, (2, 1), 4, rng=3)
    for kwargs in ({"L": 0}, {"L": 2.5}, {"pad_factor": 1.5}):
        with pytest.raises(ValueError):
            estimate_velocity(bp, **kwargs)
    with pytest.raises(ValueError):
        estimate_velocity(np.zeros((8, 8)))
    with pytest.raises(ValueError):
        estimate_velocity(np.zeros((32, 32, 3)))
    e = estimate_velocity(bp, min_quality=0.99)
    assert e.valid == (e.quality >= 0.99)
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
