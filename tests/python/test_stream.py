import numpy as np

from fftvel import simulate_pulsed_image
from fftvel.stream import FrameWindow, estimate_stream, estimate_window

from conftest import angle_error


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
    assert abs(angle_error(r.angle_deg, np.rad2deg(np.arctan2(1, 4)))) < 0.5
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
