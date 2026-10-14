"""Command line ``fftvel`` (``python -m fftvel``): help, simulate, estimate, file formats and messages.

The files are read by fftvel.io (more checks of it in test_io.py); the
command line adds the notes on stderr and the options in its messages.
"""

import builtins
import json
import re
import shutil
import subprocess
import sys

import numpy as np
import pytest
from scipy.io import savemat

import fftvel
import fftvel.__main__ as cli
from fftvel import estimate_velocity, simulate_pulsed_image
from fftvel.core import crop_square
from fftvel.io import load
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


def hide(monkeypatch, *names, module="fftvel.io"):
    """Make ``import name`` fail inside ``module`` only, as if the package were not installed.

    The readers (SciPy, h5py, Pillow, tifffile, matplotlib.image) are imported
    in fftvel.io, matplotlib for --plot in fftvel.__main__.
    """
    real = builtins.__import__

    def fake(name, globals=None, locals=None, fromlist=(), level=0):
        if globals is not None and globals.get("__name__") == module and name.split(".")[0] in names:
            raise ImportError(f"No module named '{name}'")
        return real(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake)


def to_uint16(x):
    return np.round((x - x.min()) / np.ptp(x) * 65535).astype(np.uint16)


# --- help, version, simulate and the console script ------------------------------------------

def test_help_is_plain_text_and_leads_with_the_flash_count(capsys):
    for argv in ([], ["simulate"], ["estimate"]):
        with pytest.raises(SystemExit) as exc:
            main(argv + ["--help"])
        assert exc.value.code == 0
        out = capsys.readouterr().out
        assert "``" not in out and "(default: None)" not in " ".join(out.split()), argv
        examples = out[out.index("examples:"):].splitlines()
        assert "--exposures 9" in examples[1], argv                 # the first example gives the flash count


def test_every_option_has_help_and_shows_its_default(capsys):
    _, s, e = cli._parser()
    for sub in (s, e):
        for action in sub._actions:
            assert action.help, (sub.prog, action.dest)
    with pytest.raises(SystemExit):
        main(["simulate", "--help"])
    out = " ".join(capsys.readouterr().out.split())             # the help is wrapped to the terminal width
    for text in ("(default: 3.0)", "(default: -2.0)", "(default: 256)", "(default: 0.0)", "(default: 0)",
                 "(default 4, so K = 9)"):
        assert text in out, text
    with pytest.raises(SystemExit):
        main(["estimate", "--help"])
    out = " ".join(capsys.readouterr().out.split())
    assert "(default: 512)" in out and "(default: auto)" in out


def test_version(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0 and capsys.readouterr().out.split()[-1] == fftvel.__version__


def test_program_name():
    # the console script names itself in the help, python -m fftvel otherwise
    assert cli._prog("/usr/local/bin/fftvel") == "fftvel"
    assert cli._prog("fftvel.exe") == "fftvel"
    assert cli._prog("/repo/python/src/fftvel/__main__.py") == "python -m fftvel"
    assert cli._prog("/usr/bin/pytest") == "python -m fftvel"
    assert cli._prog("") == cli._prog(None) == "python -m fftvel"


def test_simulate_writes_only_json_to_stdout(capsys, tmp_path):
    png = tmp_path / "fig.png"
    r, err = run(["simulate", "--vx", "4", "--vy", "1", "--plot", png], capsys)    # stdout is all JSON
    assert abs(r["speed"] - 4.1231) < 0.05 and r["exposures"] == 9 and r["L"] == 4
    assert "true: angle 165.96 deg, speed 4.123 px per step" in err
    assert png.stat().st_size > 10_000


def test_simulate_flash_count(capsys):
    r, _ = run(["simulate", "--exposures", "9", "--vx", "5", "--vy", "-3"], capsys)
    assert r == run(["simulate", "--L", "4", "--vx", "5", "--vy", "-3"], capsys)[0]
    assert r == run(["simulate", "--vx", "5", "--vy", "-3"], capsys)[0]          # default 9 flashes
    same(r, estimate_velocity(simulate_pulsed_image(256, (5, -3), 4, rng=0), exposures=9))
    r, _ = run(["simulate", "--exposures", "8", "--vx", "2", "--vy", "1", "--n", "128", "--seed", "3"], capsys)
    assert r["exposures"] == 8 and r["L"] is None
    same(r, estimate_velocity(simulate_pulsed_image(128, (2, 1), 3, exposures=8, rng=3), L=3, exposures=8))


@pytest.mark.parametrize("bad, text", [(["--L", "4", "--exposures", "8"], "only one"), (["--L", "0"], "--L"),
                                       (["--exposures", "1"], "--exposures"), (["--n", "8"], "--n"),
                                       (["--noise", "-1"], "--noise"), (["--vx", "nan"], "--vx"),
                                       (["--vy", "-inf"], "--vy"), (["--seed", "-1"], "--seed"),
                                       (["--vx", "1e6"], "more than 4 times --n (1024 px)"),
                                       (["--L", "20", "--vx", "26"], "move by 1040 px")])
def test_simulate_invalid_options(capsys, bad, text):
    with pytest.raises(SystemExit) as exc:
        main(["simulate"] + bad)
    cap = capsys.readouterr()
    assert exc.value.code == 2 and text in cap.err and cap.out == ""       # stopped before any work


def test_simulate_travel_limit(capsys):
    r, _ = run(["simulate", "--vx", "128", "--vy", "0", "--exposures", "9"], capsys)   # 8 * 128 = 4 * 256 px
    assert r["exposures"] == 9


def test_console_script(shot, tmp_path):
    from importlib.metadata import entry_points
    (ep,) = entry_points(group="console_scripts", name="fftvel")    # after pip install -e ./python
    assert ep.value == "fftvel.__main__:main" and ep.load() is main
    exe = shutil.which("fftvel")
    if exe is None:
        pytest.skip("the console script fftvel is not on PATH")
    r = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=120)
    assert r.returncode == 0 and r.stdout.strip() == f"fftvel {fftvel.__version__}"
    np.save(tmp_path / "bp.npy", shot[0])
    r = subprocess.run([exe, "estimate", str(tmp_path / "bp.npy"), "--exposures", "9"], capture_output=True,
                       text=True, timeout=120)
    assert r.returncode == 0, r.stderr
    same(json.loads(r.stdout), estimate_velocity(shot[0], exposures=9))
    assert "model: 9 exposures (L = 4)" in r.stderr
    r = subprocess.run([exe, "estimate", str(tmp_path / "bp.npy")], capture_output=True, text=True, timeout=120)
    assert r.returncode == 1 and "--exposures or --L" in r.stderr
    r = subprocess.run([sys.executable, "-m", "fftvel", "estimate", "--help"], capture_output=True, text=True,
                       timeout=120)
    assert r.returncode == 0 and re.search(r"^usage: python -m fftvel estimate", r.stdout, re.MULTILINE)


# --- flash count and options -----------------------------------------------------------------

def test_single_image_needs_the_flash_count(shot, tmp_path, capsys):
    bp, _ = shot
    np.save(tmp_path / "bp.npy", bp)
    msg = fails(["estimate", tmp_path / "bp.npy"], capsys)
    assert "--exposures or --L" in msg and "K_true/K" in msg
    r, err = run(["estimate", tmp_path / "bp.npy", "--L", "4"], capsys)
    same(r, estimate_velocity(bp, L=4))
    assert r["L"] == 4 and "9 exposures (L = 4)" in err
    r, _ = run(["estimate", tmp_path / "bp.npy", "--exposures", "9"], capsys)
    same(r, estimate_velocity(bp, L=4))
    assert r["L"] == 4
    with pytest.raises(SystemExit) as exc:      # checked with the other options, before the file is read
        main(["estimate", str(tmp_path / "none.npy"), "--L", "4", "--exposures", "8"])
    assert exc.value.code == 2 and "--L 4 means 9 exposures but --exposures is 8" in capsys.readouterr().err


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
        a = load(tmp_path / name)
        assert a.dtype == np.uint8 and a.shape == (128, 128, 1, 10)     # no float copy of all frames
        bp, K = cli._pulsed_image(a, True, L=4)
        assert K == 9 and np.array_equal(bp, u8[:, :, 0, :9].astype(float).mean(axis=2))
    assert isinstance(load(tmp_path / "rec.npy"), np.memmap)            # read from disk when used
    r, err = run(["estimate", tmp_path / "rec.mat", "--stack", "--L", "4"], capsys)
    same(r, estimate_velocity(u8[:, :, 0, :9].astype(float).mean(axis=2), L=4))
    assert "variable 'recording', 128 x 128 x 1 x 10 uint8" in err
    assert "mean of frames 1-9 of 10" in err


def test_stack_files_with_any_case_extension(tmp_path, capsys):
    _, fr = simulate_pulsed_image(256, (3.0, -2.0), 2, rng=3, return_frames=True)   # 5 frames
    savemat(tmp_path / "stack.mat", {"frames": fr})
    with open(tmp_path / "STACK.NPY", "wb") as f:     # np.save would append .npy
        np.save(f, fr)
    for name in ("stack.mat", "STACK.NPY"):
        r, _ = run(["estimate", tmp_path / name, "--stack"], capsys)
        assert r["L"] == 2 and abs(r["speed"] / np.hypot(3, 2) - 1) < 0.01
    assert "--stack" in fails(["estimate", tmp_path / "stack.mat"], capsys)          # 3-D without --stack
    np.save(tmp_path / "even.npy", fr[:, :, :4])
    assert "even number" in fails(["estimate", tmp_path / "even.npy", "--stack"], capsys)


# --- .npz files -----------------------------------------------------------------------------

def test_npz_as_the_live_script_writes_it(shot, tmp_path, capsys):
    bp, fr = shot
    u8 = np.round(np.clip(fr, -3, 3) * 40 + 128).astype(np.uint8)
    np.savez(tmp_path / "burst_0001.npz", frames=u8, t=np.arange(9) / 200, clock="camera", frame_rate=200.0,
             exposure_s=5e-4)
    r, err = run(["estimate", tmp_path / "burst_0001.npz", "--stack"], capsys)
    same(r, estimate_velocity(u8.astype(float).mean(axis=2), L=4))
    assert "array 'frames', 128 x 128 x 9 uint8" in err and "mean of frames 1-9 of 9" in err
    g = np.round((bp - bp.min()) / np.ptp(bp) * 255).astype(np.uint8)
    np.savez(tmp_path / "frame_0001.npz", frames=g[:, :, None], t=np.zeros(1), clock="host")    # flash mode
    r, err = run(["estimate", tmp_path / "frame_0001.npz", "--exposures", "9"], capsys)
    same(r, estimate_velocity(g, exposures=9))


def test_npz_array_choice(shot, tmp_path, capsys):
    bp, fr = shot
    np.savez(tmp_path / "two.npz", image=bp, stack=fr)
    r, err = run(["estimate", tmp_path / "two.npz", "--exposures", "9"], capsys)
    same(r, estimate_velocity(bp, exposures=9))
    assert "array 'image'" in err
    r, err = run(["estimate", tmp_path / "two.npz", "--stack"], capsys)
    assert "array 'stack'" in err and r["exposures"] == 9
    np.savez(tmp_path / "three.npz", bg=fr[:, :, 0], image=bp, stack=fr)
    msg = fails(["estimate", tmp_path / "three.npz", "--exposures", "9"], capsys)
    assert "choose one with --var" in msg and all(f"  {v}" in msg for v in ("bg ", "image ", "stack "))
    r, _ = run(["estimate", tmp_path / "three.npz", "--var", "image", "--exposures", "9"], capsys)
    same(r, estimate_velocity(bp, exposures=9))
    assert "no array 'img'" in fails(["estimate", tmp_path / "three.npz", "--var", "img", "--L", "4"], capsys)
    (tmp_path / "bad.npz").write_bytes(b"not a zip file")
    assert "cannot read" in fails(["estimate", tmp_path / "bad.npz", "--L", "4"], capsys)


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


def test_parity_case_from_mat(repo, capsys):
    # bp01 of the parity cases is a pulsed image with v = (3, -2) and 9 flashes
    r, err = run(["estimate", repo / "tests" / "parity" / "cases.mat", "--var", "bp01", "--exposures", "9"],
                 capsys)
    assert r["valid"] and "variable 'bp01'" in err


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
    assert ".mat and .npz files only" in fails(["estimate", tmp_path / "bp.npy", "--var", "bp", "--L", "4"],
                                               capsys)


def test_mat_needs_scipy(shot, tmp_path, monkeypatch, capsys):
    savemat(tmp_path / "one.mat", {"bp": shot[0]})
    hide(monkeypatch, "scipy")
    msg = fails(["estimate", tmp_path / "one.mat", "--L", "4"], capsys)
    assert "needs SciPy" in msg and "pip install scipy" in msg and '"data"' in msg


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
    a = load(tmp_path / "v73.mat")
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
    assert "v7.3" in msg and "pip install h5py" in msg and '"data"' in msg and "-v7" in msg


# --- image files ----------------------------------------------------------------------------

def test_png_and_tiff_keep_the_bit_depth(shot, tmp_path, capsys):
    Image = pytest.importorskip("PIL.Image")
    bp, _ = shot
    u16 = to_uint16(bp)
    i32 = np.round((bp - bp.min()) / np.ptp(bp) * 4000).astype(np.int32)      # 12-bit values
    f32 = ((bp - bp.min()) / np.ptp(bp)).astype(np.float32)                   # values in [0, 1]
    for name, img in (("bp16.png", u16), ("bp16.tif", u16), ("bpI.tif", i32), ("bpF.tif", f32)):
        Image.fromarray(img).save(tmp_path / name)
        a = load(tmp_path / name)
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
    a = load(tmp_path / "pages.tif")
    assert np.array_equal(a, pages)
    r, err = run(["estimate", tmp_path / "pages.tif", "--stack"], capsys)
    same(r, estimate_velocity(pages.astype(float).mean(axis=2), L=4))
    assert "9 pages stacked as 128 x 128 x 9 uint16" in err
    ims[0].save(tmp_path / "thumb.tif", save_all=True, append_images=[ims[1].resize((32, 32))])
    r, err = run(["estimate", tmp_path / "thumb.tif", "--L", "4"], capsys)
    same(r, estimate_velocity(pages[:, :, 0], L=4))
    assert "skipped 1 pages of another size" in err


def test_rgb_image_becomes_grey(shot, tmp_path, capsys):
    Image = pytest.importorskip("PIL.Image")
    bp, _ = shot
    g = np.round((bp - bp.min()) / np.ptp(bp) * 255).astype(np.uint8)
    rgb = np.stack([g, np.roll(g, 1, axis=0), g // 2], axis=2)
    for name, img in (("rgb.png", rgb), ("rgba.png", np.dstack([rgb, np.full_like(g, 255)]))):
        Image.fromarray(img).save(tmp_path / name)
        a = load(tmp_path / name)
        assert np.allclose(a, rgb @ np.array([0.2989, 0.5870, 0.1140]), rtol=0, atol=1e-12)


def test_float64_tiff_through_tifffile(shot, tmp_path, capsys):
    tifffile = pytest.importorskip("tifffile")
    pytest.importorskip("PIL")
    bp, _ = shot
    tifffile.imwrite(tmp_path / "bp64.tif", bp)                # Pillow cannot read float64
    assert np.array_equal(load(tmp_path / "bp64.tif"), bp)
    r, err = run(["estimate", tmp_path / "bp64.tif", "--L", "4"], capsys)
    same(r, estimate_velocity(bp, L=4))
    assert "tifffile" in err


def test_saturation_warning(shot, tmp_path, capsys):
    bp, _ = shot
    g = np.round(100 + 15 * bp / bp.std()).astype(np.uint8)
    np.save(tmp_path / "ok.npy", g)
    r, err = run(["estimate", tmp_path / "ok.npy", "--L", "4"], capsys)
    assert r["saturated"] == 0 and "saturated" not in err
    g[:12, :12] = 255                                       # 144 of 16384 px, 0.9 %: no warning yet
    np.save(tmp_path / "a.npy", g)
    r, err = run(["estimate", tmp_path / "a.npy", "--L", "4"], capsys)
    assert r["saturated"] == round(144 / 128**2, 4) and "saturated" not in err
    g[-12:, -12:] = 255                                     # 1.8 %
    np.save(tmp_path / "b.npy", g)
    r, err = run(["estimate", tmp_path / "b.npy", "--L", "4"], capsys)
    assert r["saturated"] == round(288 / 128**2, 4)
    assert "warning: 1.8 % of the analysed pixels are saturated" in err


def test_blur_warning(shot, tmp_path, capsys):
    # one exposure blurred over 6 px (10 exposures of duty 0.99), read as 9 flashes 0.67 px apart
    blur = simulate_pulsed_image(128, (6 / 10, 0.0), 4, exposures=10, duty=0.99, noise=0.1, rng=4)
    np.save(tmp_path / "blur.npy", blur)
    r, err = run(["estimate", tmp_path / "blur.npy", "--L", "4"], capsys)
    assert r["valid"] and r["blur_like"] and "fits better than the flashes" in err
    np.save(tmp_path / "blur_stack.npy", np.repeat(blur[:, :, None], 9, axis=2))   # mean = the blurred image
    r, err = run(["estimate", tmp_path / "blur_stack.npy", "--stack"], capsys)
    assert r["valid"] and r["blur_like"] and "fits better than 9 separate frames" in err and "flash" not in err
    np.save(tmp_path / "flashes.npy", shot[0])
    r, err = run(["estimate", tmp_path / "flashes.npy", "--L", "4"], capsys)
    assert r["valid"] and not r["blur_like"] and "fits better" not in err


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
    a = load(tmp_path / "bp16.png")
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


def test_plot_path_is_checked_before_any_work(shot, tmp_path, capsys):
    pytest.importorskip("matplotlib")
    np.save(tmp_path / "bp.npy", shot[0])
    for target, text in ((tmp_path / "none" / "fig.png", "no such folder"), (tmp_path, "not the folder"),
                         (tmp_path / "fig.xyz", "cannot write .xyz files")):
        for argv in (["simulate", "--n", "64"], ["estimate", tmp_path / "bp.npy", "--L", "4"]):
            with pytest.raises(SystemExit) as exc:
                main([str(x) for x in argv] + ["--plot", str(target)])
            cap = capsys.readouterr()
            assert exc.value.code == 2 and text in cap.err and cap.out == "", (target, argv)
    _, err = run(["simulate", "--n", "64", "--plot", tmp_path / "fig"], capsys)    # matplotlib adds .png
    assert (tmp_path / "fig.png").is_file() and f"figure written to {tmp_path / 'fig.png'}" in err


def test_plot_that_cannot_be_written_gives_no_result(shot, tmp_path, monkeypatch, capsys):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    def refuse(self, *args, **kw):
        raise PermissionError(13, "Permission denied")

    monkeypatch.setattr(Figure, "savefig", refuse)
    msg = fails(["simulate", "--n", "64", "--plot", tmp_path / "fig.png"], capsys)
    assert msg == f"cannot write the figure {tmp_path / 'fig.png'}: Permission denied"
    assert capsys.readouterr().out == ""                       # no JSON before the error


def test_plot_without_matplotlib(shot, tmp_path, monkeypatch, capsys):
    np.save(tmp_path / "bp.npy", shot[0])
    hide(monkeypatch, "matplotlib", module="fftvel.__main__")
    msg = fails(["estimate", tmp_path / "bp.npy", "--L", "4", "--plot", tmp_path / "fig.png"], capsys)
    assert "--plot needs matplotlib" in msg and 'python[plot]' in msg
    assert capsys.readouterr().out == ""                       # stopped before the estimate
    assert not (tmp_path / "fig.png").exists()
