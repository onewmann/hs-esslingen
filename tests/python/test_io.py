"""fftvel.io: reading .npz files, the frames axis, grey conversion, notes, warnings and errors.

The other formats are checked through the command line in test_cli.py.
"""

import logging

import numpy as np
import pytest
from scipy.io import savemat

import fftvel
from fftvel import io
from fftvel.core import GREY_WEIGHTS


@pytest.fixture(scope="module")
def frames():
    """64 px stack of 9 frames, uint8, with v = (2, -1)."""
    _, fr = fftvel.simulate_pulsed_image(64, (2.0, -1.0), 4, rng=3, return_frames=True)
    return np.round(np.clip(fr, -3, 3) * 40 + 128).astype(np.uint8)


def test_public_names():
    assert sorted(io.__all__) == ["HDF5_SIGNATURE", "frames_axis", "load", "to_grey"]
    assert io.HDF5_SIGNATURE == b"\x89HDF\r\n\x1a\n"
    assert not hasattr(fftvel, "load") and "io" not in fftvel.__all__      # not re-exported


def test_npz_frames_first(frames, tmp_path):
    # the format of python/examples/live_camera.py --save: frames, t, clock, frame_rate, exposure_s
    np.savez(tmp_path / "burst.npz", frames=frames, t=np.arange(9) / 200, clock="camera", frame_rate=200.0,
             exposure_s=5e-4)
    for stack in (False, True):
        a = io.load(tmp_path / "burst.npz", stack=stack)
        assert a.dtype == np.uint8 and np.array_equal(a, frames)
    np.savez(tmp_path / "bg.npz", frames=frames, background=frames[:, :, 0])    # frames wins without stack too
    assert np.array_equal(io.load(tmp_path / "bg.npz"), frames)
    assert np.array_equal(io.load(tmp_path / "bg.npz", var="background"), frames[:, :, 0])


def test_npz_without_frames(frames, tmp_path):
    img = frames[:, :, 0]
    np.savez(tmp_path / "one.npz", image=img, t=np.zeros(1))                # the only array with 2 axes > 1
    assert np.array_equal(io.load(tmp_path / "one.npz"), img)
    np.savez(tmp_path / "two.npz", image=img, stack=frames)                 # 2-D, or 3-D with stack
    assert np.array_equal(io.load(tmp_path / "two.npz"), img)
    assert np.array_equal(io.load(tmp_path / "two.npz", stack=True), frames)
    np.savez(tmp_path / "pos.npz", frames)                                  # np.savez(file, frames): arr_0
    assert np.array_equal(io.load(tmp_path / "pos.npz"), frames)


def test_npz_errors_are_neutral(frames, tmp_path):
    np.savez(tmp_path / "three.npz", a=frames[:, :, 0], b=frames[:, :, 1], c=frames)
    with pytest.raises(ValueError) as exc:
        io.load(tmp_path / "three.npz")
    msg = str(exc.value)
    assert "several numeric arrays" in msg and "under the name frames" in msg and "--" not in msg
    assert "  a  64 x 64 uint8" in msg and "  c  64 x 64 x 9 uint8" in msg
    with pytest.raises(ValueError, match="no array 'd'"):
        io.load(tmp_path / "three.npz", var="d")
    np.savez(tmp_path / "text.npz", frames=np.array(["a", "b"]))
    with pytest.raises(ValueError, match="no numeric array"):
        io.load(tmp_path / "text.npz")
    np.savez(tmp_path / "cplx.npz", frames=np.ones((16, 16), complex))
    with pytest.raises(ValueError, match="not a numeric array|no numeric array"):
        io.load(tmp_path / "cplx.npz")
    np.savez(tmp_path / "obj.npz", frames=np.array([{"a": 1}], dtype=object), image=frames[:, :, 0])
    assert np.array_equal(io.load(tmp_path / "obj.npz"), frames[:, :, 0])  # entries that need pickle are skipped
    (tmp_path / "bad.npz").write_bytes(b"not a zip file")
    with pytest.raises(ValueError, match="cannot read .* as a .npz file"):
        io.load(tmp_path / "bad.npz")


def test_errors(frames, tmp_path):
    with pytest.raises(FileNotFoundError, match="no such file"):
        io.load(tmp_path / "none.npy")
    np.save(tmp_path / "a.npy", frames)
    with pytest.raises(ValueError, match=".mat and .npz files only"):
        io.load(tmp_path / "a.npy", var="frames")
    savemat(tmp_path / "ab.mat", {"a": frames[:, :, 0], "b": frames[:, :, 1]})
    with pytest.raises(ValueError) as exc:
        io.load(tmp_path / "ab.mat")
    assert "several numeric arrays" in str(exc.value) and "--" not in str(exc.value)
    assert np.array_equal(io.load(tmp_path / "ab.mat", var="b"), frames[:, :, 1])


def test_notes_go_to_the_logger(frames, tmp_path, caplog):
    np.savez(tmp_path / "burst.npz", frames=frames)
    savemat(tmp_path / "rec.mat", {"recording": frames})
    with caplog.at_level(logging.INFO, logger="fftvel.io"):
        io.load(tmp_path / "burst.npz")
        io.load(tmp_path / "rec.mat")
    assert [r.name for r in caplog.records] == ["fftvel.io", "fftvel.io"]
    assert caplog.records[0].getMessage().endswith("burst.npz: array 'frames', 64 x 64 x 9 uint8")
    assert caplog.records[1].getMessage().endswith("rec.mat: variable 'recording', 64 x 64 x 9 uint8")


def test_jpeg_gives_a_warning(frames, tmp_path):
    Image = pytest.importorskip("PIL.Image")
    Image.fromarray(frames[:, :, 0]).save(tmp_path / "f.jpg", quality=95)
    with pytest.warns(UserWarning, match="JPEG .* lossless"):
        a = io.load(tmp_path / "f.jpg")
    assert a.shape == (64, 64) and a.dtype == np.uint8


def test_frames_axis():
    assert io.frames_axis((64, 64, 9)) == -1 and io.frames_axis((9, 64, 64)) == 0
    assert io.frames_axis((64, 64, 32)) == -1 and io.frames_axis((32, 64, 64)) == 0
    assert io.frames_axis((64, 9, 64), "first") == 0 and io.frames_axis((64, 9, 64), "last") == -1
    with pytest.raises(ValueError) as exc:
        io.frames_axis((64, 64, 41))
    assert "64 x 64 x 41" in str(exc.value) and "--" not in str(exc.value)
    with pytest.raises(ValueError, match="choice"):
        io.frames_axis((64, 64, 9), "middle")
    with pytest.raises(ValueError, match="3-D"):
        io.frames_axis((64, 64))


def test_to_grey():
    rng = np.random.default_rng(0)
    rgb = rng.integers(0, 256, (8, 10, 3)).astype(np.uint8)
    cr, cg, cb = GREY_WEIGHTS
    grey = cr * rgb[..., 0] + cg * rgb[..., 1] + cb * rgb[..., 2]
    assert np.array_equal(io.to_grey(rgb), grey)
    assert np.array_equal(io.to_grey(np.dstack([rgb, np.zeros((8, 10), np.uint8)])), grey)    # alpha ignored
    g = io.to_grey(rgb[..., :1])
    assert g.shape == (8, 10) and g.dtype == float and np.array_equal(g, rgb[..., 0])
    assert io.to_grey(rgb[..., 0]).dtype == float
    with pytest.raises(ValueError, match="shape"):
        io.to_grey(rgb[..., :2])
    from fftvel.stream import to_grey                     # the same function, for the live script
    assert to_grey is io.to_grey
