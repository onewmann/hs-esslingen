"""The phone sessions in data/measurements: the videos are read upright with their timestamps, and meta.csv and the
results match them."""

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


GRID = REPO / "data" / "measurements" / "2026-10-17_phone_grid"


@pytest.fixture
def grid():
    spec = importlib.util.spec_from_file_location("grid_extract", GRID / "extract_frames.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_grid_period_of_a_drawn_grid(grid):
    pytest.importorskip("scipy")
    img = np.full((512, 512), 120.0)
    for x in np.arange(10.0, 512, 124.5):              # thin dark lines 124.5 px apart, both ways
        img[:, int(round(x))] = 60
        img[int(round(x)), :] = 60
    assert grid.grid_period(img) == pytest.approx(124.5, abs=1.0)


def test_grid_scale_skips_misread_windows(grid, monkeypatch):
    periods = iter([125, 125, 126, 160, 126, 127])     # one misreading among a slow drift
    monkeypatch.setattr(grid, "grid_period", lambda img: next(periods))
    gsd, _ = grid.window_gsd(np.zeros((6, 4, 4)), window=1)
    assert gsd[3] is None and gsd[0] == pytest.approx(5 / 125) and gsd[5] == pytest.approx(5 / 127)


def test_grid_line_angle_of_a_turned_grid(grid):
    pytest.importorskip("scipy")
    rows, cols = np.mgrid[0:1024, 0:1024].astype(float)
    x, y = cols, -rows                                  # y up, as angle_deg in results.csv
    for deg in (-3.0, 0.0, 4.0):
        a = np.radians(deg)
        u = (-x * np.sin(a) + y * np.cos(a)) % 130      # across the lines that run at deg
        v = (x * np.cos(a) + y * np.sin(a)) % 130       # across the lines at right angles to them
        img = np.where((np.minimum(u, 130 - u) < 1) | (np.minimum(v, 130 - v) < 1), 60.0, 120.0)
        assert grid.line_angle(img) == pytest.approx(deg, abs=0.1)


def _read(name):
    with open(GRID / name, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def test_grid_meta_and_results(grid):
    meta = _read("meta.csv")
    assert len(meta) == 198 and [r["static"] for r in meta].count("1") == 89
    moving = [r for r in meta if r["static"] == "0"]
    assert [i for i, r in enumerate(moving) if r["gsd_mm_per_px"] == ""] == [107, 108]   # the paper's edge
    assert all(0.036 < float(r["gsd_mm_per_px"]) < 0.041 for r in moving[:107])
    lines = _read("grid_lines.csv")
    assert [r["file"] for r in lines] == [r["file"] + f"#{r['first_frame']}-{int(r['first_frame']) + 4}"
                                          for r in meta]
    res = _read("results.csv")
    assert len(res) == 198
    even = [r for r in res if r["file"].startswith("frames_moving") and r["uneven_steps"] == "0"]
    valid = [r for r in even if r["valid"] == "1"]
    off = [abs((float(r["angle_deg"]) - float(r["ref_angle_deg"]) + 90) % 180 - 90) for r in valid]
    assert len(even) == 91 and len(valid) == 83 and sum(d > 2 for d in off) == 10
    assert [r["file"] for r in res if r["static"] == "1" and r["valid"] == "1"] == [
        f"frames_still.npz#{s}-{s + 4}" for s in (306, 311, 316)]


def test_grid_directions_against_the_lines(capsys):
    spec = importlib.util.spec_from_file_location("grid_directions", GRID / "directions.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.main(GRID) == 0
    out = capsys.readouterr().out
    assert "grid lines at -0.1 to 4.1 deg" in out and "closer to the lines than the step: 10" in out
    assert "in 62 of 83" in out


def test_grid_lines_use_the_frames_of_the_window(grid, monkeypatch, tmp_path):
    pytest.importorskip("cv2")
    field = np.arange(12.0)[:, None, None] * np.ones((12, 4, 4))     # frame i (from 0) holds the value i
    monkeypatch.setattr(grid, "_reader", lambda: lambda *args: (field, np.arange(12) / 120))
    monkeypatch.setattr(grid, "grid_period", lambda img: 130.0)
    monkeypatch.setattr(grid, "line_angle", lambda img: float(img[0, 0]))
    assert grid.main(["--out", str(tmp_path)]) == 0
    with open(tmp_path / "grid_lines.csv", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    # window s (from 1) holds the frames s to s + 4, the values s - 1 to s + 3, mean s + 1
    assert [float(r["line_angle_deg"]) for r in rows][:2] == [2.0, 7.0]
