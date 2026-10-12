"""experiments/real_measurements.py: meta.csv checks, stacks, the summary and a run on data/real/example."""

import csv
import math
import shutil
import subprocess
import sys

import numpy as np
import pytest

from fftvel import simulate_pulsed_image


@pytest.fixture
def rm(repo, monkeypatch):
    """The module real_measurements, imported from experiments/ and removed again afterwards."""
    monkeypatch.setattr(sys, "path", [str(repo / "experiments"), *sys.path])
    for name in ("real_measurements", "common"):
        monkeypatch.delitem(sys.modules, name, raising=False)
    import real_measurements
    yield real_measurements
    for name in ("real_measurements", "common"):
        sys.modules.pop(name, None)


def write_meta(path, header, *rows, delimiter=",", encoding="utf-8"):
    with open(path, "w", newline="", encoding=encoding) as f:
        w = csv.writer(f, delimiter=delimiter, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    return path


def read_results(path):
    with open(path, newline="") as f:
        return {r["file"]: r for r in csv.DictReader(f)}


def test_example_folder(repo):
    folder = repo / "data" / "real" / "example"
    assert sum(p.stat().st_size for p in folder.iterdir()) < 1_000_000
    from PIL import Image
    with open(folder / "meta.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 3 and sum(r["static"] == "1" for r in rows) == 1
    for name in [r["file"] for r in rows] + ["ref_a.png", "ref_b.png"]:
        with Image.open(folder / name) as im:
            assert im.mode == "I;16" and im.size == (256, 256), name
    assert not (folder / "results.csv").exists()


def test_run_on_example(repo, tmp_path):
    example = shutil.copytree(repo / "data" / "real" / "example", tmp_path / "example")
    figure = tmp_path / "fig.png"
    r = subprocess.run([sys.executable, str(repo / "experiments" / "real_measurements.py"),
                        str(example / "meta.csv"), "--figure", str(figure), "--jobs", "2"],
                       capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr
    res = read_results(example / "results.csv")
    assert set(res) == {"flash_a.png", "flash_b.png", "static.png"}
    for name, angle in (("flash_a.png", 30), ("flash_b.png", 120)):
        row = res[name]
        assert row["valid"] == "1" and row["error"] == "" and row["warnings"] == "", row
        assert abs(float(row["speed_err"])) < 0.01
        assert abs(float(row["angle_err_deg"])) < 0.5
        assert float(row["speed_mm_s"]) == pytest.approx(
            float(row["speed_px"]) * float(row["gsd_mm_per_px"]) / float(row["flash_period_s"]))
        assert abs(float(row["angle_deg"]) - angle) < 0.5
    a = res["flash_a.png"]
    assert abs(float(a["ref_speed_diff"])) < 0.01 and abs(float(a["ref_angle_diff_deg"])) < 0.5
    assert float(a["ref_speed_mm_s"]) == pytest.approx(200, rel=0.01)
    assert res["static.png"]["valid"] == "0" and res["static.png"]["reason"]
    out = r.stdout
    assert "Moving images marked valid: 2 of 2" in out
    assert "Static images marked valid: 0 of 1 (should be 0; 95 % upper bound of the rate: 95.0 %)" in out
    assert "Valid but off by more than 5 % or 2 deg: 0 of 2" in out
    assert "Phase correlation of the reference frames (1 valid image)" in out
    assert figure.stat().st_size > 10_000


def test_upper_bound(rm):
    for n in (1, 10, 30, 100):
        assert rm.upper_bound(0, n) == pytest.approx(1 - 0.05 ** (1 / n))
    assert rm.upper_bound(3, 3) == 1.0
    assert math.isnan(rm.upper_bound(0, 0))
    assert rm.upper_bound(0, 20) < rm.upper_bound(1, 20) < rm.upper_bound(2, 20)
    assert rm.upper_bound(1, 20) == pytest.approx(0.2161, abs=1e-4)     # Clopper-Pearson table value


def test_bad_rows_are_reported_and_skipped(rm, tmp_path, capsys):
    folder = tmp_path / "session 1"                       # a path with a space, as on Windows desktops
    folder.mkdir()
    np.save(folder / "a b.npy", simulate_pulsed_image(128, (3.0, 0.0), 4, rng=1))
    header = ["file", "exposures", "duty", "set_speed_mm_s", "static", "roi", "ref_frame_a", "note"]
    meta = write_meta(folder / "meta.csv", header,
                      ["a b.npy", "1", "0", "", "0", "", "", ""],
                      ["a b.npy", "9", "1.2", "", "0", "", "", ""],
                      ["missing.png", "9", "", "", "", "", "", ""],
                      ["a b.npy", "9", "", "0", "0", "", "", ""],
                      ["a b.npy", "9", "", "", "2", "8", "", ""],
                      ["a b.npy", "", "", "", "", "", "a b.npy", ""],
                      ["a b.npy", "9", "0", "2", "", "0", "", "", "", ""],   # decimal comma: 0,2 in two cells
                      ["a b.npy", "9", "", "", "", "", "", "\u0394t 1 ms \u2192 3 px"])
    assert rm.main([str(meta)]) == 1                      # some rows not analysed
    out = capsys.readouterr().out
    assert "7 of 8 rows of meta.csv skipped because of problems; the other rows are analysed:" in out
    for text in ("line 2: exposures '1' must be an integer >= 2", "line 3: duty '1.2' must be a number in [0, 1)",
                 "line 4: file not found: missing.png", "line 5: set_speed_mm_s is 0 but static is 0",
                 "line 6: static '2' must be 0 or 1", "line 6: roi '8' must be 0",
                 "line 7: exposures is empty", "line 7: give both ref_frame_a and ref_frame_b",
                 "line 8: 10 cells but the header has 8 columns; a decimal comma"):
        assert text in out, text
    assert "a b.npy: " in out and "Not analysed (7;" in out
    with open(folder / "results.csv", newline="", encoding="utf-8") as f:
        header_out = next(csv.reader(f))
        f.seek(0)
        rows = list(csv.DictReader(f))
    assert tuple(header_out) == rm.FIELDS and header_out[:3] == ["file", "line", "valid"]
    assert [r["line"] for r in rows] == [str(k) for k in range(2, 10)]
    assert all(r["error"].startswith("row skipped: ") for r in rows[:7])
    assert rows[7]["valid"] == "1" and rows[7]["error"] == "" and rows[7]["note"] == "\u0394t 1 ms \u2192 3 px"


def test_unusable_meta_or_output_stops_before_any_estimate(rm, tmp_path, capsys):
    np.save(tmp_path / "a.npy", np.zeros((32, 32)))
    meta = write_meta(tmp_path / "meta.csv", ["name", "K"], ["a.npy", "9"])
    assert rm.main([str(meta)]) == 2
    assert "meta.csv has no column file, exposures" in capsys.readouterr().err
    meta = write_meta(tmp_path / "meta.csv", ["file", "exposures"], ["a.npy", "9"])
    assert rm.main([str(meta), "--out", str(tmp_path / "a.npy" / "results.csv")]) == 2
    assert "cannot write" in capsys.readouterr().err
    assert not (tmp_path / "results.csv").exists()


def test_meta_formats(rm, tmp_path):
    (tmp_path / "a.png").write_bytes(b"")
    header = ["file", "exposures", "gsd_mm_per_px", "set_angle_deg", "note", "operator"]
    meta = write_meta(tmp_path / "meta.csv", header,
                      ["a.png", "9", "0,05", "-30", "Flash \u00fcber Rolle", "me"],
                      ["#a.png", "9", "", "", "excluded: axis still accelerating", ""],
                      ["", "", "", "", "", ""],
                      delimiter=";", encoding="cp1252")
    rows, commented, problems, unknown = rm.read_meta(meta)
    assert problems == [] and unknown == ["operator"]
    assert commented == ["line 3: #a.png"]
    (row,) = rows
    assert row["gsd_mm_per_px"] == 0.05 and row["set_angle_deg"] == -30 and row["note"] == "Flash \u00fcber Rolle"
    assert row["roi"] == rm.DEFAULT_ROI and row["static"] == 0 and row["duty"] == 0.0
    assert row["flash_period_s"] is None and row["ref"] is None
    (tmp_path / "sub").mkdir()
    meta = write_meta(tmp_path / "comma.csv", ["file", "exposures"], ["sub\\..\\a.png", "8"])
    rows, _, problems, _ = rm.read_meta(meta, default_roi=0)
    assert problems == [] and rows[0]["roi"] == 0 and rows[0]["path"].is_file()


def test_stack_gives_one_image_per_frame(rm, tmp_path):
    frames = [simulate_pulsed_image(128, (3.0 * k, -1.0), 2, rng=k) for k in (1, 2, 3)]
    np.save(tmp_path / "stack.npy", np.stack(frames, axis=2).astype(np.float32))
    meta = write_meta(tmp_path / "meta.csv", ["file", "exposures", "flash_period_s", "gsd_mm_per_px"],
                      ["stack.npy", "5", "0.001", "0.1"], ["stack.npy#2", "5", "0.001", "0.1"],
                      ["stack.npy#4", "5", "", ""])
    rows, _, problems, _ = rm.read_meta(meta)
    assert problems == []
    tasks, failed = rm.build_tasks(rows)
    assert [t["name"] for t in tasks] == ["stack.npy#1", "stack.npy#2", "stack.npy#3", "stack.npy#2"]
    assert all(t["image"].shape == (128, 128) for t in tasks)
    assert len(failed) == 1 and "has 3 frames, not frame 4" in failed[0]["error"]
    assert tasks[1]["warnings"] == [] and tasks[3]["warnings"] == [
        "the same image as in line 2, so it counts twice in the summary"]
    assert tasks[0]["image"].base is None                 # a copy, not a view of the whole frame
    res = rm.analyse(tasks[3])
    assert res["valid"] == 1 and res["speed_px"] == pytest.approx(np.hypot(6, 1), rel=0.01)
    assert res["speed_mm_s"] == pytest.approx(res["speed_px"] * 0.1 / 0.001)


def test_lowest_speed_is_the_estimator_limit(rm):
    from fftvel import estimate_velocity
    from fftvel.core import GUARD_FACTOR, SPEED_STEP
    kext = int(np.floor(np.log(GUARD_FACTOR) / SPEED_STEP))     # the speed grid starts at vmin / GUARD_FACTOR
    for n, K in ((64, 9), (96, 4)):
        e = estimate_velocity(simulate_pulsed_image(n, (2.0, 0.0), 4, rng=3), exposures=K, diagnostics=True)
        assert rm.lowest_speed(K, n) == pytest.approx(e.diag["speed_grid"][kext], rel=1e-12)
    assert rm.lowest_speed(9, 512) == pytest.approx(3.35 / 9, rel=0.01)


def test_warnings_and_summary(rm, tmp_path):
    img = simulate_pulsed_image(128, (2.0, 0.0), 4, rng=5)
    img = np.clip(np.round(128 + 30 * img), 0, 254).astype(np.uint8)
    img[:10, :10] = 255                                   # two saturated patches, 1.2 % of the pixels
    img[-10:, -10:] = 255
    np.save(tmp_path / "a.npy", img)
    meta = write_meta(tmp_path / "meta.csv",
                      ["file", "exposures", "flash_period_s", "gsd_mm_per_px", "set_speed_mm_s",
                       "set_angle_deg", "exposure_s"],
                      ["a.npy", "9", "0.001", "0.05", "100", "0", "0.005"],
                      ["a.npy", "9", "0.001", "0.05", "1", "0", ""])
    rows, _, problems, _ = rm.read_meta(meta)
    assert problems == []
    tasks, failed = rm.build_tasks(rows)
    results = [rm.analyse(t) for t in tasks]
    w0, w1 = results[0]["warnings"], results[1]["warnings"]
    assert "shorter than (K - 1 + duty) * T" in w0 and "1.2 % of the analysed pixels are saturated" in w0
    assert "at the maximum 255" in w0 and results[0]["saturated"] == pytest.approx(200 / 128**2)
    assert "the set speed is 0.02 px per flash, outside" in w1
    assert results[0]["expected_px"] == pytest.approx(2.0)
    lines = []
    s = rm.summary(results, log=lines.append)
    assert s["moving"] == s["moving_valid"] == s["checked"] == 2 and s["static"] == 0
    assert s["off"] == 1 and s["speed_n"] == 2 and s["speed_max"] > 1000      # the second row: set 1 mm/s
    assert s["off_bound"] == pytest.approx(rm.upper_bound(1, 2))
    assert any(ln.startswith("Static images: none.") for ln in lines)
    assert any(ln.startswith("Valid but off by more than 5 % or 2 deg: 1 of 2") for ln in lines)
