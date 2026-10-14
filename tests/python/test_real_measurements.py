"""experiments/real_measurements.py: meta.csv checks, stacks, the summary, average mode, the files it reads and runs
on the synthetic sessions in data/measurements."""

import csv
import math
import shutil
import subprocess
import sys

import numpy as np
import pytest

from fftvel import make_texture, simulate_pulsed_image


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


def test_synthetic_flash_folder(repo):
    folder = repo / "data" / "measurements" / "synthetic_flash"
    assert sum(p.stat().st_size for p in folder.iterdir()) < 1_000_000
    from PIL import Image
    with open(folder / "meta.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 3 and sum(r["static"] == "1" for r in rows) == 1
    for name in [r["file"] for r in rows] + ["ref_a.png", "ref_b.png"]:
        with open(folder / name, "rb") as f:
            head = f.read(26)
        assert head[12:16] == b"IHDR" and head[24:26] == bytes([16, 0]), name     # 16-bit grey PNG
        with Image.open(folder / name) as im:
            # Pillow >= 10.3 opens 16-bit grey as I;16, older versions as I (32-bit)
            assert im.mode in ("I;16", "I") and im.size == (256, 256), name
            assert np.asarray(im).max() > 255, name
    assert not (folder / "results.csv").exists()


def test_run_on_synthetic_flash(repo, tmp_path):
    session = shutil.copytree(repo / "data" / "measurements" / "synthetic_flash", tmp_path / "synthetic_flash")
    figure = tmp_path / "fig.png"
    r = subprocess.run([sys.executable, str(repo / "experiments" / "real_measurements.py"),
                        str(session / "meta.csv"), "--figure", str(figure), "--jobs", "2"],
                       capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr
    res = read_results(session / "results.csv")
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


def test_files_are_read_with_fftvel_io(rm, repo, tmp_path):
    # the evaluation reads files with the public fftvel.io, not with private names of the command line, so
    # its messages name no command line options
    source = (repo / "experiments" / "real_measurements.py").read_text(encoding="utf-8")
    assert "fftvel.__main__" not in source and "SystemExit" not in source
    bp = simulate_pulsed_image(96, (3.0, -1.0), 4, rng=2)
    np.savez(tmp_path / "image.npz", image=bp)                                  # one array, not named frames
    np.savez(tmp_path / "two.npz", a=bp, b=bp)                                  # two images: which one?
    np.save(tmp_path / "cube.npy", np.repeat(bp[:48, :48, None], 40, axis=2))  # 48 x 48 x 40: frames axis?
    meta = write_meta(tmp_path / "meta.csv", ["file", "exposures"], ["image.npz", "9"], ["two.npz", "9"],
                      ["cube.npy", "9"])
    rows, _, problems, _ = rm.read_meta(meta)
    assert problems == []
    tasks, failed = rm.build_tasks(rows)
    assert [t["name"] for t in tasks] == ["image.npz"] and rm.analyse(tasks[0])["valid"] == 1
    errors = [r["error"] for r in failed]
    assert "two.npz holds several numeric arrays" in errors[0] and "--var" not in errors[0]
    assert "give the name" not in errors[0] and "save only the image" in errors[0] and "np.savez" in errors[0]
    assert "cannot tell which axis of the 48 x 48 x 40 stack" in errors[1] and "--frames-axis" not in errors[1]


def test_a_file_that_cannot_be_read_fails_only_its_own_row(rm, tmp_path, monkeypatch):
    bp = simulate_pulsed_image(96, (3.0, -1.0), 4, rng=2)
    np.save(tmp_path / "good.npy", bp)
    (tmp_path / "empty.npy").write_bytes(b"")                                  # 0 bytes: EOFError in NumPy
    (tmp_path / "odd.npy").write_bytes(b"x")
    meta = write_meta(tmp_path / "meta.csv", ["file", "exposures"], ["empty.npy", "9"], ["odd.npy", "9"],
                      ["good.npy", "9"])
    rows, _, problems, _ = rm.read_meta(meta)
    real = rm.fftvel.io.load

    def load(path, **kw):                       # a loader error of a type that fftvel.io does not promise
        if path.endswith("odd.npy"):
            raise RuntimeError("damaged")
        return real(path, **kw)

    monkeypatch.setattr(rm.fftvel.io, "load", load)
    rm._read.cache_clear()
    tasks, failed = rm.build_tasks(rows)
    assert [t["name"] for t in tasks] == ["good.npy"]
    assert "cannot read" in failed[0]["error"] and "empty.npy" in failed[0]["error"]
    assert failed[1]["error"].endswith("odd.npy (RuntimeError: damaged)")


def test_colour_arrays_get_a_warning(rm, tmp_path):
    from scipy.io import savemat
    bp = simulate_pulsed_image(96, (3.0, -1.0), 4, rng=2)
    rgb = np.round(np.stack([bp, bp, bp], axis=2) * 20 + 128).astype(np.uint8)
    np.save(tmp_path / "colour.npy", rgb)
    savemat(tmp_path / "getdata.mat", {"recording": rgb[:, :, None, :]})      # rows x cols x 1 x frames
    from PIL import Image
    Image.fromarray(rgb).save(tmp_path / "colour.png")                        # the loader makes it grey
    meta = write_meta(tmp_path / "meta.csv", ["file", "exposures"], ["colour.npy#2", "9"], ["getdata.mat", "9"],
                      ["colour.png", "9"])
    rows, _, _, _ = rm.read_meta(meta)
    tasks, failed = rm.build_tasks(rows)
    assert failed == [] and [t["name"] for t in tasks] == ["colour.npy#2", "getdata.mat#1", "getdata.mat#2",
                                                           "getdata.mat#3", "colour.png"]
    assert tasks[0]["warnings"] == ["the file holds 96 x 96 x 3: if it is a colour image, convert it to grey; "
                                    "otherwise its 3 frames count as separate images"]
    assert all(t["warnings"] == [] for t in tasks[1:])


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
    colour = ("the file holds 128 x 128 x 3: if it is a colour image, convert it to grey; otherwise its 3 frames "
              "count as separate images")
    assert tasks[1]["warnings"] == [colour] and tasks[3]["warnings"] == [
        colour, "the same image as in line 2, so it counts twice in the summary"]
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


# Average mode

# values of the run on data/measurements/synthetic_flash before average mode was added (flash mode must
# give the same)
OLD_FLASH_ROWS = {
    "flash_a.png": dict(valid=1, speed_px=4.002197201325098, angle_deg=30.017069240380067,
                        quality=0.9455707004007519, speed_err=0.0005493003312744271,
                        ref_speed_diff=0.0004216654838808953),
    "flash_b.png": dict(valid=1, speed_px=5.596782005853841, angle_deg=120.10808733669346,
                        quality=0.9099130304993869, speed_err=-0.0005746418118142138),
}
AVERAGE_COLUMNS = ("frames", "first_frame", "frame_period_s", "measured_period_s", "period_source", "step_dx_px",
                   "step_dy_px", "step_speed_px", "step_spread", "step_angle_spread_deg", "uneven_steps")


def shifted_frames(positions, n=96, seed=0):
    """Frames of a periodic random texture at the given (x, y) positions in px (exact Fourier shifts)."""
    m = n + 32
    tex = np.fft.fft2(make_texture(m, rng=seed))
    f = np.fft.fftfreq(m)
    fx, fy = np.meshgrid(f, f)
    frames = [np.real(np.fft.ifft2(tex * np.exp(-2j * np.pi * (fx * x + fy * y))))[16:16 + n, 16:16 + n]
              for x, y in positions]
    return np.stack(frames, axis=2)


def counts(frames, rng=0, noise=1.0):
    """Texture frames as 8-bit counts with noise."""
    g = np.random.default_rng(rng)
    return np.clip(np.round(100 + 25 * frames + noise * g.standard_normal(frames.shape)), 0, 255).astype(np.uint8)


def test_synthetic_average_folder(repo):
    folder = repo / "data" / "measurements" / "synthetic_average"
    assert sum(p.stat().st_size for p in folder.iterdir()) < 1_000_000
    with open(folder / "meta.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    assert [r["mode"] for r in rows] == ["average"] * 3 and sum(r["static"] == "1" for r in rows) == 1
    for r in rows:
        with np.load(folder / r["file"]) as z:
            assert z["frames"].shape == (192, 192, 9) and z["frames"].dtype == np.uint8, r["file"]
            assert z["t"].shape == (9,) and float(z["frame_rate"]) == 200 and float(z["exposure_s"]) == 0.0005
            assert str(z["clock"]) == "camera"
    assert not (folder / "results.csv").exists()


def test_run_on_synthetic_average(repo, tmp_path):
    session = shutil.copytree(repo / "data" / "measurements" / "synthetic_average", tmp_path / "synthetic_average")
    r = subprocess.run([sys.executable, str(repo / "experiments" / "real_measurements.py"),
                        str(session / "meta.csv"), "--jobs", "2"], capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr
    out = r.stdout
    res = read_results(session / "results.csv")
    assert set(res) == {"moving_regular.npz", "moving_dropped.npz", "static.npz"}
    reg, drop, static = res["moving_regular.npz"], res["moving_dropped.npz"], res["static.npz"]
    # regular: valid, close to the true 6 px per frame at 30 deg, and close to the mean step
    assert reg["valid"] == "1" and reg["error"] == "" and reg["warnings"] == "", reg
    assert abs(float(reg["speed_err"])) < 0.005 and abs(float(reg["angle_err_deg"])) < 0.3
    assert reg["mode"] == "average" and reg["exposures"] == reg["frames"] == "9" and reg["first_frame"] == "1"
    assert float(reg["duty"]) == pytest.approx(0.1)                      # exposure_s / frame_period_s
    assert float(reg["step_dx_px"]) == pytest.approx(6 * np.cos(np.pi / 6), abs=0.02)
    assert float(reg["step_dy_px"]) == pytest.approx(-3, abs=0.02)
    assert float(reg["step_spread"]) < 0.02 and float(reg["step_angle_spread_deg"]) < 0.5
    assert reg["uneven_steps"] == "0" and reg["period_source"] == "frame_period_s"
    assert float(reg["measured_period_s"]) == pytest.approx(0.005)
    assert abs(float(reg["ref_speed_diff"])) < 0.005 and abs(float(reg["ref_angle_diff_deg"])) < 0.3
    assert float(reg["ref_speed_mm_s"]) == pytest.approx(60, rel=0.005)
    # one frame missing: flagged by the timestamps and by the steps; the period comes from the timestamps
    assert drop["uneven_steps"] == "1" and float(drop["step_spread"]) == pytest.approx(1, abs=0.02)
    assert "uneven frame intervals in the timestamps: 5 to 10 ms" in drop["warnings"]
    assert "uneven steps between consecutive frames (phase correlation)" in drop["warnings"]
    assert drop["period_source"] == "timestamps" and float(drop["measured_period_s"]) == pytest.approx(0.005)
    assert float(drop["speed_mm_s"]) == pytest.approx(float(drop["speed_px"]) * 0.05 / 0.005)
    # stopped surface: not valid; steps too short to check, but the timestamps are regular
    assert static["valid"] == "0" and static["reason"] and static["uneven_steps"] == "0"
    assert float(static["step_speed_px"]) < 0.05 and static["step_spread"] == "" and static["ref_speed_diff"] == ""
    assert "Average mode (averaged frames): 3 windows, 1 with uneven steps (column uneven_steps)" in out
    assert "  moving_dropped.npz  uneven steps" in out
    assert "Frame period from the timestamps (frame_period_s empty) for 2 windows" in out
    assert "Static images marked valid: 0 of 1" in out
    # the step comparison leaves out the window with uneven steps, whose mean step is 12.5 % too long
    assert ("Mean step between consecutive frames, phase correlation (1 valid window with even steps): angle "
            "difference median") in out
    assert abs(float(drop["ref_speed_diff"]) + 1 / 9) < 0.01
    assert "moving_dropped.npz: 60 mm/s (6.000 px/frame, frame period 5 ms from the timestamps)" in out


def test_old_synthetic_flash_rows_unchanged(rm, repo, tmp_path):
    folder = repo / "data" / "measurements" / "synthetic_flash"
    rows, _, problems, _ = rm.read_meta(folder / "meta.csv")
    assert problems == []
    results = [rm.analyse(t) for t in rm.build_tasks(rows)[0]]
    by_file = {r["file"]: r for r in results}
    for name, expected in OLD_FLASH_ROWS.items():
        for key, value in expected.items():
            assert by_file[name][key] == pytest.approx(value, rel=1e-7), (name, key)
    assert by_file["static.png"]["valid"] == 0 and by_file["static.png"]["reason"].startswith("quality 0.29")
    for r in results:
        assert r["mode"] == "flash" and all(r[k] is None for k in AVERAGE_COLUMNS), r["file"]
    # the same meta.csv with the new columns, mode flash and the average-mode columns empty: identical rows
    with open(folder / "meta.csv", newline="") as f:
        old = list(csv.reader(f))
    meta = write_meta(tmp_path / "meta.csv", old[0] + ["mode", "frames", "first_frame", "frame_period_s"],
                      *[r + ["flash", "", "", ""] for r in old[1:]])
    for name in ("flash_a.png", "flash_b.png", "static.png", "ref_a.png", "ref_b.png"):
        shutil.copy(folder / name, tmp_path / name)
    rows2, _, problems, _ = rm.read_meta(meta)
    assert problems == []
    results2 = [rm.analyse(t) for t in rm.build_tasks(rows2)[0]]
    assert results2 == results


def test_average_mat_from_getdata_and_host_times(rm, tmp_path):
    from scipy.io import savemat
    v = (2.5, -1.0)
    _, fr = simulate_pulsed_image(96, v, 1, exposures=11, rng=7, return_frames=True)
    frames = counts(fr, rng=8)[:, :, None, :]                       # getdata: H x W x 1 x F, uint8
    jitter = np.array([0, 1.5, -1.2, 0.8, -1.5, 1.0, 0.3, -0.9, 1.4, -0.4, 0.6]) * 1e-3
    t = (0.01 * np.arange(11) + jitter)[:, None]                    # host times, up to 1.5 ms off
    for name, clock in (("host.mat", "host"), ("camera.mat", "camera")):
        savemat(tmp_path / name, dict(frames=frames, t=t, exposure_s=1e-4, frame_rate=100.0, clock=clock))
    savemat(tmp_path / "manual.mat", dict(data=frames, time=t))     # [data, time] = getdata(...), saved by hand
    meta = write_meta(tmp_path / "meta.csv",
                      ["file", "mode", "frames", "first_frame", "frame_period_s", "gsd_mm_per_px", "set_speed_mm_s",
                       "set_angle_deg"],
                      ["host.mat", "average", "9", "2", "", "0.1", f"{np.hypot(*v) * 10:g}",
                       f"{np.degrees(np.arctan2(1, 2.5)):g}"],
                      ["camera.mat", "average", "9", "2", "", "", "", ""],
                      ["host.mat", "average", "", "", "", "", "", ""],               # 11 frames, all averaged
                      ["manual.mat", "average", "9", "2", "0.01", "0.1", "", ""],
                      ["manual.mat", "average", "9", "", "", "0.1", "", ""])
    rows, _, problems, _ = rm.read_meta(meta)
    assert problems == []
    tasks, failed = rm.build_tasks(rows)
    assert failed == [] and [t["name"] for t in tasks] == ["host.mat#2-10", "camera.mat#2-10", "host.mat",
                                                           "manual.mat#2-10", "manual.mat#1-9"]
    assert tasks[0]["frames"].shape == (96, 96, 9) and tasks[0]["frames"].dtype == np.uint8
    assert np.array_equal(tasks[0]["frames"], frames[:, :, 0, 1:10]) and tasks[2]["K"] == 11
    assert np.allclose(tasks[0]["times"], t[1:10, 0]) and tasks[0]["clock"] == "host"
    assert np.allclose(tasks[3]["times"], t[1:10, 0]) and tasks[3]["clock"] == "host"    # a .mat without clock
    host, camera, whole, manual, manual_free = (rm.analyse(task) for task in tasks)
    # host times: accepted against the frame rate in the file, whose period gives the speed, as in live_basler
    assert host["valid"] == 1 and host["uneven_steps"] == 0 and host["period_source"] == "frame_rate"
    assert host["speed_px"] == pytest.approx(np.hypot(*v), rel=0.01) and abs(host["angle_err_deg"]) < 0.5
    assert host["speed_mm_s"] == pytest.approx(host["speed_px"] * 0.1 / 0.01)
    assert host["measured_period_s"] == pytest.approx((t[9, 0] - t[1, 0]) / 8)       # the mean interval
    assert abs(host["speed_err"]) < 0.01 and host["warnings"] == ""
    assert host["exposure_s"] == pytest.approx(1e-4) and host["duty"] == pytest.approx(0.01)   # from the file
    # the same times in a .mat without clock count as host times; frame_period_s gives the nominal period
    assert manual["uneven_steps"] == 0 and manual["period_source"] == "frame_period_s"
    assert manual["speed_mm_s"] == pytest.approx(host["speed_mm_s"])
    assert "frames 2 to 9 are also averaged in line 5" in manual_free["warnings"]
    # without a nominal period the spread of the host times (3 ms of 10 ms) is flagged, and the warning says why
    assert manual_free["uneven_steps"] == 1 and manual_free["period_source"] == "timestamps"
    assert "uneven frame intervals in the host times" in manual_free["warnings"]
    assert "give the camera's frame period in frame_period_s" in manual_free["warnings"]
    assert manual_free["measured_period_s"] == pytest.approx((t[8, 0] - t[0, 0]) / 8)
    assert host["step_dx_px"] == pytest.approx(v[0], abs=0.03) and host["step_dy_px"] == pytest.approx(v[1], abs=0.03)
    # the same times taken as exact camera times: their spread (3 ms of 10 ms) is flagged, the steps are even
    assert camera["uneven_steps"] == 1 and camera["step_spread"] < 0.05
    assert "uneven frame intervals" in camera["warnings"] and "uneven steps between" not in camera["warnings"]
    assert whole["exposures"] == 11 and whole["valid"] == 1


def test_average_host_times_with_a_dropped_frame(rm, tmp_path):
    from scipy.io import savemat
    positions = [k for k in range(15) if k != 6]                    # 14 frames, the 7th lost
    frames = counts(shifted_frames([(2.0 * k, 1.0 * k) for k in positions]))[:, :, None, :]
    jitter = 1e-3 * np.array([0, 0.6, -0.5, 0.4, -0.7, 0.2, 0.5, -0.3, 0.6, -0.4, 0.1, 0.7, -0.6, 0.3])
    savemat(tmp_path / "rec.mat", dict(frames=frames, t=0.01 * np.array(positions, float) + jitter))   # no clock
    meta = write_meta(tmp_path / "meta.csv", ["file", "mode", "frames", "first_frame", "frame_period_s"],
                      ["rec.mat", "average", "9", "1", "0.01"],       # 1 to 9: the mean interval 11 % long
                      ["rec.mat", "average", "13", "2", "0.01"],      # 2 to 14: 8 % long, one interval 20 ms
                      ["rec.mat", "average", "5", "10", "0.01"])      # 10 to 14: no gap (frames of line 3)
    rows, _, problems, _ = rm.read_meta(meta)
    assert problems == []
    short, long_, after = (rm.analyse(t) for t in rm.build_tasks(rows)[0])
    assert short["uneven_steps"] == 1 and long_["uneven_steps"] == 1 and after["uneven_steps"] == 0
    assert "uneven frame intervals in the host times" in short["warnings"]
    assert "differs from the frame period 10 ms by more than 10 %: a frame was dropped" in short["warnings"]
    assert long_["warnings"].startswith("frames 2 to 9 are also averaged in line 2")
    assert "a frame interval in the host times is off by half the frame period (10 ms) or more" in long_["warnings"]
    assert after["period_source"] == "frame_period_s" and after["valid"] == 1
    assert after["warnings"] == "frames 10 to 14 are also averaged in line 3, so the two windows are not independent"


def test_average_mat73_timestamps(rm, tmp_path):
    h5py = pytest.importorskip("h5py")
    _, fr = simulate_pulsed_image(96, (2.0, -1.5), 4, rng=2, return_frames=True)
    frames = counts(fr)[:, :, None, :]
    with h5py.File(tmp_path / "rec.mat", "w") as f:            # MATLAB v7.3 stores arrays transposed
        for name, value, cls in (("frames", frames, "uint8"), ("t", 0.01 * np.arange(9)[:, None], "double"),
                                 ("clock", np.array([[ord(c)] for c in "host"], dtype=np.uint16), "char"),
                                 ("frame_rate", np.array([[100.0]]), "double")):
            f.create_dataset(name, data=value.T).attrs["MATLAB_class"] = np.bytes_(cls)
    assert rm._extras(str(tmp_path / "rec.mat")) == dict(t=pytest.approx(0.01 * np.arange(9)), clock="host",
                                                         frame_rate=100.0)
    meta = write_meta(tmp_path / "meta.csv", ["file", "mode"], ["rec.mat", "average"])
    rows, _, problems, _ = rm.read_meta(meta)
    (r,) = [rm.analyse(t) for t in rm.build_tasks(rows)[0]]
    assert r["valid"] == 1 and r["period_source"] == "frame_rate" and r["measured_period_s"] == pytest.approx(0.01)
    assert r["speed_px"] == pytest.approx(2.5, rel=0.01) and r["warnings"] == ""


def test_average_stack_options_and_warnings(rm, tmp_path):
    # 10 frames, frames first, 3 px per frame along x; a .npy holds no timestamps
    stack = counts(shifted_frames([(3.0 * k, 0.0) for k in range(10)]))
    np.save(tmp_path / "even.npy", np.moveaxis(stack, 2, 0))
    np.save(tmp_path / "one.npy", stack[:, :, 0])
    # equal step lengths, but the direction turns by 10 deg after four steps
    turn = [(3.0 * k, 0.0) for k in range(5)] + [(12 + 3 * np.cos(np.radians(10)) * k,
                                                  -3 * np.sin(np.radians(10)) * k) for k in range(1, 5)]
    np.save(tmp_path / "turn.npy", counts(shifted_frames(turn, seed=1)))
    header = ["file", "mode", "frames", "first_frame", "frame_period_s", "exposure_s", "duty", "gsd_mm_per_px"]
    meta = write_meta(tmp_path / "meta.csv", header,
                      ["even.npy", "average", "", "", "", "", "", ""],               # even number of frames
                      ["even.npy", "average", "9", "3", "", "", "", ""],             # frames 3 to 11 of 10
                      ["one.npy", "average", "", "", "", "", "", ""],                # not a stack
                      ["even.npy", "average", "9", "2", "0.01", "0.001", "", "0.05"],
                      ["even.npy", "average", "9", "", "0.01", "0.005", "", ""],     # duty 0.5
                      ["even.npy", "average", "9", "", "", "0.003", "", ""],         # no period: duty unknown
                      ["even.npy", "average", "7", "", "0.01", "0.003", "0.02", ""],  # duty given
                      ["turn.npy", "average", "", "", "", "", "", ""])
    rows, _, problems, _ = rm.read_meta(meta)
    assert problems == []
    tasks, failed = rm.build_tasks(rows)
    errors = [r["error"] for r in failed]
    assert "the stack has 10 frames, an even number; give frames" in errors[0]
    assert "needs frames 3 to 11, the stack has 10" in errors[1]
    assert "average mode needs a stack of frames, the file holds an array of 96 x 96" in errors[2]
    res = [rm.analyse(t) for t in tasks]
    a, half, unknown, given, turned = res
    assert a["file"] == "even.npy#2-10" and a["valid"] == 1 and a["period_source"] == "frame_period_s"
    assert a["duty"] == pytest.approx(0.1) and a["speed_px"] == pytest.approx(3, rel=0.01)
    assert a["speed_mm_s"] == pytest.approx(a["speed_px"] * 0.05 / 0.01)
    assert a["ref_speed_mm_s"] == pytest.approx(15, rel=0.01)
    assert a["measured_period_s"] is None and a["uneven_steps"] == 0 and a["warnings"] == ""
    assert half["duty"] == pytest.approx(0.5) and "Between 0.4 and 0.6" in half["warnings"]
    assert unknown["duty"] == 0 and unknown["period_source"] is None and unknown["speed_mm_s"] is None
    assert given["duty"] == pytest.approx(0.02) and given["exposures"] == 7
    assert given["warnings"] == "frames 2 to 7 are also averaged in line 5, so the two windows are not independent"
    assert unknown["warnings"].startswith("the same frames as in line 6, so they count twice in the summary")
    assert turned["uneven_steps"] == 1 and turned["step_spread"] < 0.05
    assert 9.5 < turned["step_angle_spread_deg"] < 10.5 and "directions spread over 10." in turned["warnings"]
    lines = []
    s = rm.summary(res, log=lines.append)
    assert s["average"] == 5 and s["uneven"] == 1
    lines = []
    rm.summary([res[0], dict(res[2], duty=0.3)], log=lines.append)
    assert "Speed error |est/set - 1|: no valid moving image with a set speed, gsd and frame period" in lines


def test_average_duty_warning_from_exposure(rm, tmp_path):
    np.save(tmp_path / "s.npy", counts(shifted_frames([(2.0 * k, 1.0 * k) for k in range(9)])))
    meta = write_meta(tmp_path / "meta.csv", ["file", "mode", "frame_period_s", "exposure_s"],
                      ["s.npy", "average", "0.01", "0.003"])
    rows, _, problems, _ = rm.read_meta(meta)
    assert problems == []
    (r,) = [rm.analyse(t) for t in rm.build_tasks(rows)[0]]
    assert r["duty"] == pytest.approx(0.3) and "0.30 of the frame period, more than a quarter" in r["warnings"]


def test_average_row_problems(rm, tmp_path):
    np.save(tmp_path / "s.npy", np.zeros((3, 32, 32), dtype=np.uint8))
    header = ["file", "mode", "exposures", "frames", "flash_period_s", "frame_period_s", "exposure_s",
              "ref_frame_a", "ref_frame_b"]
    meta = write_meta(tmp_path / "meta.csv", header,
                      ["s.npy", "avg", "9", "", "", "", "", "", ""],
                      ["s.npy", "", "9", "3", "", "0.01", "", "", ""],
                      ["s.npy#2", "average", "", "", "", "", "", "", ""],
                      ["s.npy", "average", "", "", "", "", "", "s.npy", "s.npy"],
                      ["s.npy", "average", "9", "7", "", "", "", "", ""],
                      ["s.npy", "average", "", "", "0.02", "0.01", "", "", ""],
                      ["s.npy", "average", "", "", "", "0.01", "0.01", "", ""],
                      ["s.npy", "average", "", "1", "", "", "", "", ""],
                      ["s.npy", "AVERAGE", "", "", "", "", "", "", ""])
    rows, _, problems, _ = rm.read_meta(meta)
    for text in ("line 2: mode 'avg' must be flash or average",
                 "line 3: frames, frame_period_s are for mode average only",
                 "line 4: s.npy#2: file#k selects one frame; in average mode give first_frame and frames",
                 "line 5: ref_frame_a and ref_frame_b are for mode flash",
                 "line 6: exposures 9 and frames 7 differ",
                 "line 7: frame_period_s 0.01 and flash_period_s 0.02 differ",
                 "line 8: exposure_s 0.01 is not shorter than the frame period 0.01 s",
                 "line 9: frames '1' must be an integer >= 2"):
        assert any(p.startswith(text) for p in problems), text
    assert len(problems) == 8 and rows[-1]["mode"] == "average" and rows[-1]["problems"] == []
    # with a mode column, exposures may be left out of the header
    meta = write_meta(tmp_path / "meta2.csv", ["file", "mode"], ["s.npy", "average"])
    rows, _, problems, _ = rm.read_meta(meta)
    assert problems == [] and rows[0]["exposures"] is None
    meta = write_meta(tmp_path / "meta3.csv", ["file", "mode"], ["s.npy", "flash"])
    assert rm.read_meta(meta)[2] == ["line 2: exposures is empty; give the number of flashes in the image"]
