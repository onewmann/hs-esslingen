import numpy as np
import pytest

from fftvel import angle_difference, simulate_pulsed_image
from fftvel.stream import FrameWindow, estimate_frame, estimate_stream, estimate_window, frame_step


def _frames(times, v=(4.0, -1.0), L=4, seed=0):
    _, fr = simulate_pulsed_image(160, v, L, rng=seed, return_frames=True)
    for k, t in enumerate(times):
        yield fr[:, :, k], t


def test_regular_stream_gives_valid_estimate():
    times = np.arange(9) / 30.0
    out = list(estimate_stream(_frames(times), L=4))
    assert len(out) == 1
    r = out[0]
    assert r.regular and r.valid
    assert abs(r.speed / np.hypot(4, 1) - 1) < 0.01
    assert abs(angle_difference(r.angle_deg, np.rad2deg(np.arctan2(1, 4)))) < 0.5
    assert abs(r.speed_per_s / (30 * np.hypot(4, 1)) - 1) < 0.01


def test_uneven_frame_steps_are_flagged():
    times = np.array([0, 1, 2, 3, 5, 6, 7, 8, 9]) / 30.0     # one frame interval doubled
    r = next(estimate_stream(_frames(times), L=4))
    assert not r.regular and not r.valid


def test_window_crop_and_colour():
    w = FrameWindow(L=1, roi=50)
    rgb = np.zeros((80, 120, 3))
    rgb[..., 1] = 1.0
    for t in (0.0, 0.1, 0.2):
        w.push(rgb, t)
    bp, step, regular = w.pulsed_image()
    assert bp.shape == (50, 50)
    assert np.allclose(bp, 0.5870)
    assert abs(step - 0.1) < 1e-12 and regular


def test_single_channel_frames_and_bursts():
    times = np.arange(9) / 25.0
    burst = [(f[:, :, None], t) for f, t in _frames(times)]      # (H, W, 1) frames
    r = estimate_window(burst, L=4)
    assert r.regular and r.valid
    assert abs(r.speed / np.hypot(4, 1) - 1) < 0.01


def test_alternating_intervals_are_flagged():
    # intervals 1, 1.3, 1, 1.3, ...: the median of the 8 intervals lies in the middle, so a test of
    # the largest deviation from it would pass them; the spread (max - min) does not
    times = np.cumsum([0, 1, 1.3, 1, 1.3, 1, 1.3, 1, 1.3]) / 30.0
    r = next(estimate_stream(_frames(times), L=4))
    assert not r.regular and not r.valid
    assert frame_step([0, 1.0, 2.45])[1] is False            # L = 1: intervals 1 and 1.45
    assert frame_step(np.cumsum([0] + [1, 1.09] * 4))[1] is True
    assert frame_step(np.zeros(9)) == (0.0, False)
    assert frame_step(np.cumsum([0, 1, 1.3, 1, 1.3]), tol=0.35)[1] is True


def test_frame_step_with_nominal_period():
    rng = np.random.default_rng(0)
    host = 0.002 * np.arange(9) + rng.uniform(-0.0004, 0.0004, 9)     # 500 fps, host jitter 0.4 ms
    assert frame_step(host)[1] is False
    step, regular = frame_step(host, period=0.002)
    assert regular and step == 0.002
    for drop in (np.r_[0:5, 6:10], np.r_[0:9, 10:18]):               # one dropped frame, K = 9 and 17
        step, regular = frame_step(0.002 * drop, period=0.002)
        assert not regular
    with pytest.raises(ValueError, match="period"):
        frame_step(host, period=0)
    with pytest.raises(ValueError, match="tol"):
        frame_step(host, tol="0.1")


def test_flash_frame():
    # one camera frame with 9 strobe flashes, 3.61 px apart, flashes 1 ms apart
    bp = simulate_pulsed_image(200, (3.0, -2.0), 4, rng=2)
    r = estimate_frame(np.repeat(bp[:, :, None], 3, axis=2), t=0.5, roi=160, flash_period=1e-3)
    assert r.valid and r.regular and r.t == 0.5 and r.step == 1e-3
    assert abs(r.speed / np.hypot(3, 2) - 1) < 0.01
    assert abs(r.speed_per_s / (1e3 * np.hypot(3, 2)) - 1) < 0.01
    assert r.estimate.n == 160 and r.estimate.exposures == 9
    r = estimate_frame(bp)
    assert r.valid and np.isnan(r.speed_per_s) and np.isnan(r.step)
    bp8 = simulate_pulsed_image(160, (1.0, 2.0), 4, rng=3, exposures=8, duty=0.3)
    r = estimate_frame(bp8, exposures=8, duty=0.3)
    assert r.valid and abs(r.speed / np.hypot(1, 2) - 1) < 0.01
    with pytest.raises(ValueError, match="flash_period"):
        estimate_frame(bp, flash_period=-1)


def test_saturation_warning():
    bp = simulate_pulsed_image(128, (3.0, -2.0), 4, rng=2)
    frame = np.round(100 + 15 * bp / bp.std()).astype(np.uint8)
    r = estimate_frame(frame)
    assert r.warning == "" and r.estimate.saturated == 0
    frame[:20, :20] = 255                                   # 2.4 % of the pixels
    with pytest.warns(UserWarning, match="saturated"):
        r = estimate_frame(frame)
    assert r.warning == "2.4 % of the analysed pixels are saturated"
    assert r.estimate.saturated == pytest.approx(400 / 128**2)
    burst = [(frame, k / 50.0) for k in range(9)]           # average mode: the mean keeps the 255s
    with pytest.warns(UserWarning, match="saturated"):
        r = estimate_window(burst, L=4)
    assert r.warning.startswith("2.4 %")


def test_window_options():
    _, fr = simulate_pulsed_image(160, (4.0, -1.0), 4, rng=1, return_frames=True, exposures=8, duty=0.2)
    burst = [(fr[:, :, k], k / 50.0) for k in range(8)]
    r = estimate_window(burst, L=4, exposures=8, duty=0.2, clip_outliers=True)
    assert r.regular and r.valid and r.estimate.exposures == 8
    assert abs(r.speed / np.hypot(4, 1) - 1) < 0.01 and abs(r.step - 0.02) < 1e-12
    out = list(estimate_stream(iter(burst * 2), L=4, exposures=8))
    assert len(out) == 2                                         # windows of K = 8 frames
    for bad in (dict(L="4"), dict(exposures=1), dict(roi=8), dict(tol=-0.1), dict(period="0.01")):
        with pytest.raises(ValueError):
            FrameWindow(**bad)
