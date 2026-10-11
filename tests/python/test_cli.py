"""Command line ``python -m fftvel estimate``: options, file formats and error messages."""

import builtins
import json

import numpy as np
import pytest
from scipy.io import savemat

import fftvel.__main__ as cli
from fftvel import crop_square, estimate_velocity, simulate_pulsed_image
from fftvel.__main__ import main

FIELDS = ("angle_deg", "speed", "quality", "valid", "exposures", "n")


@pytest.fixture(scope="module")
def shot():
    """128 px pulsed image with v = (3, -2), L = 4, and its 9 frames."""
    return simulate_pulsed_image(128, (3.0, -2.0), 4, rng=1, return_frames=True)


def run(argv, capsys):
    """Run the CLI; return the JSON result (all of stdout) and the stderr text."""
    assert main([str(x) for x in argv]) == 0
    cap = capsys.readouterr()
    return json.loads(cap.out), cap.err


def fails(argv, capsys):
    """Run the CLI, expect it to stop; return its message (SystemExit text or argparse stderr)."""
    with pytest.raises(SystemExit) as exc:
        main([str(x) for x in argv])
    code = exc.value.code
    return code if isinstance(code, str) else capsys.readouterr().err


def same(r, est):
    """CLI result equals the library result (the CLI rounds to 4 digits)."""
    for k in FIELDS:
        v = getattr(est, k)
        assert r[k] == (round(v, 4) if isinstance(v, float) else v), k


def hide(monkeypatch, *names):
    """Make ``import name`` fail inside fftvel.__main__ only, as if the package were not installed."""
    real = builtins.__import__

    def fake(name, globals=None, locals=None, fromlist=(), level=0):
        if globals is not None and globals.get("__name__") == "fftvel.__main__" \
                and name.split(".")[0] in names:
            raise ImportError(f"No module named '{name}'")
        return real(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake)


def to_uint16(x):
    return np.round((x - x.min()) / np.ptp(x) * 65535).astype(np.uint16)


# --- flash count and options -----------------------------------------------------------------

def test_single_image_needs_the_flash_count(shot, tmp_path, capsys):
    bp, _ = shot
    np.save(tmp_path / "bp.npy", bp)
    msg = fails(["estimate", tmp_path / "bp.npy"], capsys)
    assert "--L or --exposures" in msg and "K_true/K" in msg
    r, err = run(["estimate", tmp_path / "bp.npy", "--L", "4"], capsys)
    same(r, estimate_velocity(bp, L=4))
    assert r["L"] == 4 and "9 exposures (L = 4)" in err
    r, _ = run(["estimate", tmp_path / "bp.npy", "--exposures", "9"], capsys)
    same(r, estimate_velocity(bp, L=4))
    assert r["L"] == 4
    assert "only one" in fails(["estimate", tmp_path / "bp.npy", "--L", "4", "--exposures", "8"], capsys)


def test_even_flash_count(tmp_path, capsys):
    bp = simulate_pulsed_image(128, (1.5, 1.0), 4, exposures=8, rng=2)
    np.save(tmp_path / "bp8.npy", bp)
    r, err = run(["estimate", tmp_path / "bp8.npy", "--exposures", "8"], capsys)
    same(r, estimate_velocity(bp, exposures=8))
    assert r["exposures"] == 8 and r["L"] is None and r["valid"]
    assert "model: 8 exposures," in err


def test_options_reach_the_estimator(shot, tmp_path, monkeypatch, capsys):
    bp, _ = shot
    np.save(tmp_path / "bp.npy", bp)
    calls = []

    def spy(*args, **kw):
        calls.append(kw)
        return estimate_velocity(*args, **kw)

    monkeypatch.setattr(cli, "estimate_velocity", spy)
    r, err = run(["estimate", tmp_path / "bp.npy", "--L", "4", "--duty", "0.3", "--roi", "64",
                  "--clip-outliers"], capsys)
    assert calls[-1] == {"L": 4, "diagnostics": False, "exposures": 9, "duty": 0.3, "roi": 64,
                         "clip_outliers": True}
    assert r["n"] == 64 and "outliers clipped" in err and "duty 0.3" in err
    run(["estimate", tmp_path / "bp.npy", "--L", "4"], capsys)
    assert calls[-1]["roi"] == 512 and calls[-1]["duty"] == 0 and calls[-1]["clip_outliers"] is False
    run(["estimate", tmp_path / "bp.npy", "--L", "4", "--roi", "0"], capsys)
    assert calls[-1]["roi"] is None


@pytest.mark.parametrize("bad", [["--L", "0"], ["--exposures", "1"], ["--duty", "1"], ["--duty", "-0.1"],
                                 ["--roi", "8"], ["--roi", "-1"], ["--frames-axis", "first"]])
def test_invalid_options(shot, tmp_path, capsys, bad):
    np.save(tmp_path / "bp.npy", shot[0])
    flash = [] if bad[0] == "--L" else ["--L", "4"]
    with pytest.raises(SystemExit) as exc:
        main(["estimate", str(tmp_path / "bp.npy")] + flash + bad)
    assert exc.value.code == 2
    assert bad[0] in capsys.readouterr().err


def test_roi_default_and_crop(tmp_path, capsys):
    bp = simulate_pulsed_image(200, (3.0, -2.0), 4, rng=4)[20:180, :]       # 160 x 200
    np.save(tmp_path / "wide.npy", bp)
    r, err = run(["estimate", tmp_path / "wide.npy", "--L", "4", "--roi", "64"], capsys)
    same(r, estimate_velocity(crop_square(bp, 64), L=4))
    assert "central 64 x 64 px of the 160 x 200 image (--roi 64)" in err
    r, err = run(["estimate", tmp_path / "wide.npy", "--L", "4"], capsys)
    same(r, estimate_velocity(bp, L=4))
    assert r["n"] == 160 and "(--roi 512)" in err


# --- frame stacks ---------------------------------------------------------------------------

def test_stack_frames_first_or_last(shot, tmp_path, capsys):
    bp, fr = shot
    expected = estimate_velocity(fr.mean(axis=2), L=4)
    np.save(tmp_path / "hwf.npy", fr)
    np.save(tmp_path / "fhw.npy", np.moveaxis(fr, 2, 0))          # np.stack(list of frames)
    for name, axis in (("hwf.npy", "last"), ("fhw.npy", "first")):
        r, err = run(["estimate", tmp_path / name, "--stack"], capsys)
        same(r, expected)
        assert f"mean of frames 1-9 of 9, frames on the {axis} axis" in err
        r, _ = run(["estimate", tmp_path / name, "--stack", "--frames-axis", axis], capsys)
        same(r, expected)
    np.save(tmp_path / "middle.npy", np.moveaxis(fr, 2, 1))
    assert "--frames-axis" in fails(["estimate", tmp_path / "middle.npy", "--stack"], capsys)


def test_stack_with_unclear_axis_needs_the_flag(tmp_path, capsys):
    _, fr = simulate_pulsed_image(64, (2.0, 1.0), 20, rng=5, return_frames=True)     # 64 x 64 x 41
    np.save(tmp_path / "deep.npy", fr)
    msg = fails(["estimate", tmp_path / "deep.npy", "--stack"], capsys)
    assert "64 x 64 x 41" in msg and "--frames-axis" in msg
    r, err = run(["estimate", tmp_path / "deep.npy", "--stack", "--frames-axis", "last"], capsys)
    assert r["exposures"] == 41 and r["L"] == 20


def test_stack_errors_name_the_flags(shot, tmp_path, capsys):
    _, fr = shot
    np.save(tmp_path / "even.npy", fr[:, :, :8])
    msg = fails(["estimate", tmp_path / "even.npy", "--stack"], capsys)
    assert "8 frames, an even number" in msg and "--L" in msg and "--exposures" in msg
    msg = fails(["estimate", tmp_path / "even.npy", "--stack", "--L", "5"], capsys)
    assert "11 frames" in msg and "128 x 128 x 8" in msg
    assert "--stack" in fails(["estimate", tmp_path / "even.npy", "--L", "4"], capsys)
    r, err = run(["estimate", tmp_path / "even.npy", "--stack", "--exposures", "8"], capsys)
    same(r, estimate_velocity(fr[:, :, :8].mean(axis=2), L=3, exposures=8))
    assert r["L"] is None and "mean of frames 1-8 of 8" in err


def test_getdata_stack_converts_only_the_frames_used(shot, tmp_path, capsys):
    _, fr = shot
    u8 = np.round(np.clip(fr, -3, 3) * 40 + 128).astype(np.uint8)
    u8 = np.concatenate([u8, u8[:, :, :1]], axis=2)[:, :, None, :]   # rows x cols x 1 x 10, as getdata
    savemat(tmp_path / "rec.mat", {"recording": u8})
    np.save(tmp_path / "rec.npy", u8)
    for name in ("rec.mat", "rec.npy"):
        a = cli._load(str(tmp_path / name))
        assert a.dtype == np.uint8 and a.shape == (128, 128, 1, 10)     # no float copy of all frames
        bp, K = cli._pulsed_image(a, True, L=4)
        assert K == 9 and np.array_equal(bp, u8[:, :, 0, :9].astype(float).mean(axis=2))
    assert isinstance(cli._load(str(tmp_path / "rec.npy")), np.memmap)  # read from disk when used
    r, err = run(["estimate", tmp_path / "rec.mat", "--stack", "--L", "4"], capsys)
    same(r, estimate_velocity(u8[:, :, 0, :9].astype(float).mean(axis=2), L=4))
    assert "variable 'recording', 128 x 128 x 1 x 10 uint8" in err
    assert "mean of frames 1-9 of 10" in err


# --- .mat files -----------------------------------------------------------------------------

def test_mat_variable_choice(shot, tmp_path, capsys):
    bp, fr = shot
    savemat(tmp_path / "one.mat", {"L": 4.0, "bp": bp})          # one array of at least 2-D
    r, err = run(["estimate", tmp_path / "one.mat", "--L", "4"], capsys)
    same(r, estimate_velocity(bp, L=4))
    assert "variable 'bp'" in err
    savemat(tmp_path / "ts.mat", {"t": np.arange(9.0), "frames": fr})
    r, err = run(["estimate", tmp_path / "ts.mat", "--stack"], capsys)
    assert "variable 'frames'" in err and r["exposures"] == 9
    savemat(tmp_path / "two.mat", {"bp": bp, "frames": fr})       # the 2-D one, or the 3-D one for --stack
    r, err = run(["estimate", tmp_path / "two.mat", "--L", "4"], capsys)
    same(r, estimate_velocity(bp, L=4))
    assert "variable 'bp'" in err
    r, err = run(["estimate", tmp_path / "two.mat", "--stack"], capsys)
    same(r, estimate_velocity(fr.mean(axis=2), L=4))
    assert "variable 'frames'" in err
    savemat(tmp_path / "ws.mat", {"bg": fr[:, :, 0], "bp": bp, "frames": fr})
    msg = fails(["estimate", tmp_path / "ws.mat", "--L", "4"], capsys)
    assert "--var" in msg and all(f"  {v}" in msg for v in ("bg ", "bp ", "frames "))
    assert "128 x 128 x 9" in msg
    _, err = run(["estimate", tmp_path / "ws.mat", "--stack"], capsys)
    assert "variable 'frames'" in err
    r, err = run(["estimate", tmp_path / "ws.mat", "--var", "bp", "--L", "4"], capsys)
    same(r, estimate_velocity(bp, L=4))
    msg = fails(["estimate", tmp_path / "ws.mat", "--var", "frame", "--stack"], capsys)
    assert "no variable 'frame'" in msg and "frames" in msg


def test_mat_struct_cell_and_bad_files(shot, tmp_path, capsys):
    bp, fr = shot
    savemat(tmp_path / "s.mat", {"s": {"frames": fr, "t": 1.0}})
    msg = fails(["estimate", tmp_path / "s.mat", "--stack"], capsys)
    assert "no numeric array" in msg and "s (struct)" in msg and "s.frames" in msg
    cell = np.empty((1, 3), dtype=object)
    for k in range(3):
        cell[0, k] = fr[:, :, k]
    savemat(tmp_path / "c.mat", {"c": cell, "bp": bp})
    msg = fails(["estimate", tmp_path / "c.mat", "--var", "c", "--stack"], capsys)
    assert "class cell" in msg and "cat(3, c{:})" in msg and "  bp" in msg
    r, _ = run(["estimate", tmp_path / "c.mat", "--L", "4"], capsys)        # bp is the only numeric array
    same(r, estimate_velocity(bp, L=4))
    (tmp_path / "oct.mat").write_text("# Created by Octave 8.4.0\n# name: frames\n# type: matrix\n")
    assert "save('-v7', file, 'frames')" in fails(["estimate", tmp_path / "oct.mat", "--stack"], capsys)
    (tmp_path / "junk.mat").write_bytes(bytes(range(256)) * 4)
    assert "cannot read" in fails(["estimate", tmp_path / "junk.mat", "--stack"], capsys)
    savemat(tmp_path / "cut.mat", {"bp": bp})
    data = (tmp_path / "cut.mat").read_bytes()
    (tmp_path / "cut.mat").write_bytes(data[:len(data) // 2])   # header readable, data cut off
    assert "cannot read variable 'bp'" in fails(["estimate", tmp_path / "cut.mat", "--L", "4"], capsys)
    assert "no such file" in fails(["estimate", tmp_path / "none.mat", "--L", "4"], capsys)
    np.save(tmp_path / "bp.npy", bp)
    assert ".mat files only" in fails(["estimate", tmp_path / "bp.npy", "--var", "bp", "--L", "4"], capsys)


def test_mat_needs_scipy(shot, tmp_path, monkeypatch, capsys):
    savemat(tmp_path / "one.mat", {"bp": shot[0]})
    hide(monkeypatch, "scipy")
    assert 'python[data]' in fails(["estimate", tmp_path / "one.mat", "--L", "4"], capsys)


def write_v73(path, variables):
    """MATLAB v7.3 file as MATLAB writes it: HDF5 behind a 512-byte header, arrays transposed."""
    h5py = pytest.importorskip("h5py")
    with h5py.File(path, "w", userblock_size=512) as f:
        for name, value in variables.items():
            if isinstance(value, dict):
                g = f.create_group(name)
                g.attrs["MATLAB_class"] = np.bytes_("struct")
                for k, v in value.items():
                    g.create_dataset(k, data=np.asarray(v).T)
            else:
                d = f.create_dataset(name, data=np.ascontiguousarray(np.asarray(value).T))
                d.attrs["MATLAB_class"] = np.bytes_({"float64": "double", "uint8": "uint8"}[value.dtype.name])
        f.create_group("#refs#")
    head = b"MATLAB 7.3 MAT-file, Platform: GLNXA64, Created on: Sun Oct 11 12:00:00 2026 HDF5 schema 1.00 ."
    with open(path, "r+b") as f:
        f.write((head.ljust(116, b" ") + b"\x00" * 8 + b"\x00\x02IM").ljust(512, b"\x00"))


def test_mat_v73(shot, tmp_path, capsys):
    _, fr = shot
    fr = fr[:, :120, :]                                   # 128 x 120 x 9: checks the axis order
    write_v73(tmp_path / "v73.mat", {"frames": fr, "info": {"t": np.arange(9.0)}})
    a = cli._load(str(tmp_path / "v73.mat"))
    assert np.array_equal(a, fr)
    r, err = run(["estimate", tmp_path / "v73.mat", "--stack"], capsys)
    same(r, estimate_velocity(fr.mean(axis=2), L=4))
    assert "variable 'frames' (MATLAB v7.3), 128 x 120 x 9" in err
    msg = fails(["estimate", tmp_path / "v73.mat", "--var", "info", "--stack"], capsys)
    assert "class struct" in msg and "frames  128 x 120 x 9 double" in msg
    h5py = pytest.importorskip("h5py")
    with h5py.File(tmp_path / "v73.mat", "r+") as f:     # a sparse matrix is a group of class double
        g = f.create_group("sp")
        g.attrs["MATLAB_class"] = np.bytes_("double")
        g.attrs["MATLAB_sparse"] = np.uint64(9)
    assert "class sparse" in fails(["estimate", tmp_path / "v73.mat", "--var", "sp", "--stack"], capsys)
    _, err = run(["estimate", tmp_path / "v73.mat", "--stack"], capsys)
    assert "variable 'frames'" in err


def test_mat_v73_without_h5py(tmp_path, monkeypatch, capsys):
    head = b"MATLAB 7.3 MAT-file, Platform: PCWIN64, Created on: Sun Oct 11 12:00:00 2026 HDF5 schema 1.00 ."
    (tmp_path / "v73.mat").write_bytes((head.ljust(116, b" ") + b"\x00" * 8 + b"\x00\x02IM").ljust(512, b"\x00")
                                       + b"\x89HDF\r\n\x1a\n" + bytes(64))
    hide(monkeypatch, "h5py")
    msg = fails(["estimate", tmp_path / "v73.mat", "--stack"], capsys)
    assert "v7.3" in msg and "h5py" in msg and 'python[data]' in msg and "-v7" in msg


# --- image files ----------------------------------------------------------------------------

def test_png_and_tiff_keep_the_bit_depth(shot, tmp_path, capsys):
    Image = pytest.importorskip("PIL.Image")
    bp, _ = shot
    u16 = to_uint16(bp)
    i32 = np.round((bp - bp.min()) / np.ptp(bp) * 4000).astype(np.int32)      # 12-bit values
    f32 = ((bp - bp.min()) / np.ptp(bp)).astype(np.float32)                   # values in [0, 1]
    for name, img in (("bp16.png", u16), ("bp16.tif", u16), ("bpI.tif", i32), ("bpF.tif", f32)):
        Image.fromarray(img).save(tmp_path / name)
        a = cli._load(str(tmp_path / name))
        assert np.array_equal(a, img), name              # not scaled or clipped to 8 bit
        r, _ = run(["estimate", tmp_path / name, "--L", "4"], capsys)
        same(r, estimate_velocity(img, L=4))
        assert r["valid"] and abs(r["speed"] / np.hypot(3, 2) - 1) < 0.05, name


def test_multipage_tiff_is_a_stack(shot, tmp_path, capsys):
    Image = pytest.importorskip("PIL.Image")
    _, fr = shot
    pages = to_uint16(fr)
    ims = [Image.fromarray(pages[:, :, k]) for k in range(9)]
    ims[0].save(tmp_path / "pages.tif", save_all=True, append_images=ims[1:])
    a = cli._load(str(tmp_path / "pages.tif"))
    assert np.array_equal(a, pages)
    r, err = run(["estimate", tmp_path / "pages.tif", "--stack"], capsys)
    same(r, estimate_velocity(pages.astype(float).mean(axis=2), L=4))
    assert "9 pages stacked as 128 x 128 x 9 uint16" in err
    ims[0].save(tmp_path / "thumb.tif", save_all=True, append_images=[ims[1].resize((32, 32))])
    assert np.array_equal(cli._load(str(tmp_path / "thumb.tif")), pages[:, :, 0])
    assert "skipped 1 pages of another size" in capsys.readouterr().err


def test_rgb_image_becomes_grey(shot, tmp_path, capsys):
    Image = pytest.importorskip("PIL.Image")
    bp, _ = shot
    g = np.round((bp - bp.min()) / np.ptp(bp) * 255).astype(np.uint8)
    rgb = np.stack([g, np.roll(g, 1, axis=0), g // 2], axis=2)
    for name, img in (("rgb.png", rgb), ("rgba.png", np.dstack([rgb, np.full_like(g, 255)]))):
        Image.fromarray(img).save(tmp_path / name)
        a = cli._load(str(tmp_path / name))
        assert np.allclose(a, rgb @ np.array([0.2989, 0.5870, 0.1140]), rtol=0, atol=1e-12)


def test_float64_tiff_through_tifffile(shot, tmp_path, capsys):
    tifffile = pytest.importorskip("tifffile")
    pytest.importorskip("PIL")
    bp, _ = shot
    tifffile.imwrite(tmp_path / "bp64.tif", bp)                # Pillow cannot read float64
    assert np.array_equal(cli._load(str(tmp_path / "bp64.tif")), bp)
    r, err = run(["estimate", tmp_path / "bp64.tif", "--L", "4"], capsys)
    same(r, estimate_velocity(bp, L=4))
    assert "tifffile" in err


def test_jpeg_warning(shot, tmp_path, capsys):
    Image = pytest.importorskip("PIL.Image")
    bp, _ = shot
    g = np.round((bp - bp.min()) / np.ptp(bp) * 255).astype(np.uint8)
    Image.fromarray(g).save(tmp_path / "bp.jpg", quality=95)
    _, err = run(["estimate", tmp_path / "bp.jpg", "--L", "4"], capsys)
    assert "warning: JPEG" in err and "lossless" in err
    Image.fromarray(g).save(tmp_path / "bp.png")
    _, err = run(["estimate", tmp_path / "bp.png", "--L", "4"], capsys)
    assert "JPEG" not in err


def test_image_readers_missing(shot, tmp_path, monkeypatch, capsys):
    Image = pytest.importorskip("PIL.Image")
    pytest.importorskip("matplotlib.image")
    u16 = to_uint16(shot[0])
    Image.fromarray(u16).save(tmp_path / "bp16.png")
    hide(monkeypatch, "PIL")                                   # falls back to matplotlib
    a = cli._load(str(tmp_path / "bp16.png"))
    assert np.allclose(a, u16 / 65535, rtol=1e-6, atol=0)
    hide(monkeypatch, "PIL", "matplotlib")
    msg = fails(["estimate", tmp_path / "bp16.png", "--L", "4"], capsys)
    assert "Pillow" in msg and "pip install pillow" in msg


def test_unreadable_files(tmp_path, capsys):
    pytest.importorskip("PIL")
    (tmp_path / "bad.png").write_bytes(b"not an image")
    assert "cannot read" in fails(["estimate", tmp_path / "bad.png", "--L", "4"], capsys)
    np.save(tmp_path / "obj.npy", np.array([{"a": 1}], dtype=object), allow_pickle=True)
    assert "cannot read" in fails(["estimate", tmp_path / "obj.npy", "--L", "4"], capsys)
    np.save(tmp_path / "cplx.npy", np.ones((32, 32), complex))
    assert "complex" in fails(["estimate", tmp_path / "cplx.npy", "--L", "4"], capsys)
    np.save(tmp_path / "tiny.npy", np.ones((8, 8)))           # errors of the estimator end in one line
    assert "16x16" in fails(["estimate", tmp_path / "tiny.npy", "--L", "4"], capsys)
    nan = np.ones((32, 32))
    nan[3, 3] = np.nan
    np.save(tmp_path / "nan.npy", nan)
    assert "NaN" in fails(["estimate", tmp_path / "nan.npy", "--L", "4"], capsys)


# --- figure ---------------------------------------------------------------------------------

def test_plot_without_display(shot, tmp_path, monkeypatch, capsys):
    matplotlib = pytest.importorskip("matplotlib")
    monkeypatch.delenv("DISPLAY", raising=False)
    np.save(tmp_path / "bp.npy", shot[0])
    png = tmp_path / "fig.png"
    r, _ = run(["estimate", tmp_path / "bp.npy", "--L", "4", "--roi", "96", "--plot", png], capsys)
    assert png.stat().st_size > 10_000 and r["n"] == 96
    assert matplotlib.get_backend().lower() == "agg"


def test_plot_without_matplotlib(shot, tmp_path, monkeypatch, capsys):
    np.save(tmp_path / "bp.npy", shot[0])
    hide(monkeypatch, "matplotlib")
    msg = fails(["estimate", tmp_path / "bp.npy", "--L", "4", "--plot", tmp_path / "fig.png"], capsys)
    assert "--plot needs matplotlib" in msg and 'python[plot]' in msg
    assert capsys.readouterr().out == ""                       # stopped before the estimate
    assert not (tmp_path / "fig.png").exists()
