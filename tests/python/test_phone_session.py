"""data/measurements/2026-10-16_phone_120fps: the video is read upright with its timestamps, and meta.csv matches it."""

import csv
import importlib.util

import numpy as np
import pytest

from conftest import REPO

SESSION = REPO / "data" / "measurements" / "2026-10-16_phone_120fps"


@pytest.fixture
def ef():
    spec = importlib.util.spec_from_file_location("extract_frames", SESSION / "extract_frames.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_video_frames_and_timestamps(ef):
    pytest.importorskip("cv2")
    frames, t = ef.read_video(SESSION / "desk_120fps.mp4", max_frames=12)
    assert frames.shape == (12, 512, 512) and frames.dtype == np.uint8
    assert np.all(np.diff(t) > 0)
    assert 0.177 < t[1] < 0.179                          # the gap after the first frame
    assert np.allclose(np.diff(t[1:]), 1 / 120.2, rtol=0.01)
    # the square below the lamp's spot: mean about 89 and no saturation (a wrongly turned frame gives
    # a different square with a mean near 118 that holds part of the saturated spot)
    assert 87 < frames.mean() < 92 and frames.max() < 200 and frames.std() > 5


def test_meta_csv_lists_the_windows_of_the_video(ef, tmp_path):
    assert ef.write_meta(tmp_path / "meta.csv", 447) == 89
    assert (tmp_path / "meta.csv").read_text(encoding="utf-8") == (SESSION / "meta.csv").read_text(encoding="utf-8")
    with open(SESSION / "results.csv", newline="", encoding="utf-8") as f:
        res = list(csv.DictReader(f))
    assert len(res) == 89 and {r["mode"] for r in res} == {"average"}
