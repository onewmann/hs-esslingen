"""fftvel.io: reading .npz files, damaged files, 16-bit colour, the frames axis, grey conversion, notes,
warnings and errors.

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


def test_damaged_files_give_value_error(frames, tmp_path):
    # the promise of load: ValueError for a file that cannot be read, whatever NumPy or h5py raise
    np.save(tmp_path / "a.npy", frames)
    data = (tmp_path / "a.npy").read_bytes()
    for name, content in (("empty.npy", b""), ("cut.npy", data[:60]), ("head.npy", data[:20] + b"\x00" * 200)):
        (tmp_path / name).write_bytes(content)
        with pytest.raises(ValueError, match=f"cannot read .*{name}"):
            io.load(tmp_path / name)
    h5py = pytest.importorskip("h5py")
    with h5py.File(tmp_path / "v73.mat", "w", userblock_size=512) as f:
        for name in ("a", "b"):
            f[name] = np.ascontiguousarray(frames.T)
            f[name].attrs["MATLAB_class"] = np.bytes_("uint8")
    assert np.array_equal(io.load(tmp_path / "v73.mat", var="b"), frames)
    with pytest.raises(io._SeveralArrays):                  # the messages of the loader itself stay as they are
        io.load(tmp_path / "v73.mat")
    with pytest.raises(ValueError, match="no variable 'c'"):
        io.load(tmp_path / "v73.mat", var="c")
    (tmp_path / "cut.mat").write_bytes((tmp_path / "v73.mat").read_bytes()[:3000])   # h5py: OSError
    with pytest.raises(ValueError, match="cannot read .*cut.mat as a MATLAB v7.3 file"):
        io.load(tmp_path / "cut.mat", var="a")


def write_png16(path, a):
    """16-bit PNG written by hand (Pillow cannot write 16-bit colour): grey, grey and alpha, RGB or RGBA."""
    import struct
    import zlib
    colour_type = {1: 0, 2: 4, 3: 2, 4: 6}[1 if a.ndim == 2 else a.shape[2]]
    raw = b"".join(b"\x00" + row.astype(">u2").tobytes() for row in a)

    def chunk(kind, data):
        return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data))

    head = struct.pack(">IIBBBBB", a.shape[1], a.shape[0], 16, colour_type, 0, 0, 0)
    path.write_bytes(b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", head) + chunk(b"IDAT", zlib.compress(raw))
                     + chunk(b"IEND", b""))


def test_16_bit_colour_is_not_cut_to_8_bits(tmp_path):
    pytest.importorskip("PIL")
    rgb = np.random.default_rng(0).integers(0, 65536, (32, 40, 3)).astype(np.uint16)
    grey = io.to_grey(rgb)
    write_png16(tmp_path / "grey16.png", rgb[..., 0])
    assert np.array_equal(io.load(tmp_path / "grey16.png"), rgb[..., 0])
    for name, a in (("rgb16.png", rgb), ("rgba16.png", np.dstack([rgb, rgb[..., :1]])),
                    ("la16.png", rgb[..., :2])):
        write_png16(tmp_path / name, a)                   # Pillow would keep only the high byte
        with pytest.raises(ValueError, match="only 8 bits .* 16-bit colour"):
            io.load(tmp_path / name)
    (tmp_path / "rgb16.ppm").write_bytes(b"P6 40 32 65535\n" + rgb.astype(">u2").tobytes())
    with pytest.raises(ValueError, match="only 8 bits .* 16-bit colour"):    # Pillow scales it to 8 bits
        io.load(tmp_path / "rgb16.ppm")
    (tmp_path / "rgb8.ppm").write_bytes(b"P6 40 32 255\n" + (rgb >> 8).astype(np.uint8).tobytes())
    assert np.allclose(io.load(tmp_path / "rgb8.ppm"), io.to_grey(rgb >> 8), rtol=0, atol=1e-9)
    tifffile = pytest.importorskip("tifffile")
    tifffile.imwrite(tmp_path / "rgb16.tif", rgb, photometric="rgb")
    tifffile.imwrite(tmp_path / "planar.tif", np.moveaxis(rgb, 2, 0), photometric="rgb", planarconfig="separate")
    tifffile.imwrite(tmp_path / "pages.tif", np.stack([rgb, rgb[::-1]]), photometric="rgb")
    for name in ("rgb16.tif", "planar.tif"):
        assert np.allclose(io.load(tmp_path / name), grey, rtol=0, atol=1e-9), name
    a = io.load(tmp_path / "pages.tif")
    assert a.shape == (32, 40, 2) and np.allclose(a[..., 1], grey[::-1], rtol=0, atol=1e-9)
    tifffile.imwrite(tmp_path / "rgb8.tif", (rgb >> 8).astype(np.uint8), photometric="rgb")
    assert np.allclose(io.load(tmp_path / "rgb8.tif"), io.to_grey(rgb >> 8), rtol=0, atol=1e-9)


def test_16_bit_bmp_is_not_taken_for_16_bit_colour(tmp_path):
    # 16 bits per pixel (5-6-5, raw mode BGR;16) are fewer than 8 bits per channel: read as usual
    import struct
    Image = pytest.importorskip("PIL.Image")
    w, h = 8, 4
    pixels = (np.arange(w * h, dtype=np.uint16) * 0x0841).astype("<u2").tobytes()
    head = struct.pack("<IiiHHIIiiII", 40, w, h, 1, 16, 3, len(pixels), 2835, 2835, 0, 0)
    masks = struct.pack("<III", 0xF800, 0x07E0, 0x001F)
    offset = 14 + len(head) + len(masks)
    (tmp_path / "565.bmp").write_bytes(b"BM" + struct.pack("<IHHI", offset + len(pixels), 0, 0, offset) + head
                                       + masks + pixels)
    with Image.open(tmp_path / "565.bmp") as im:
        assert im.mode == "RGB" and im.tile[0][3][0] == "BGR;16"
        expected = io.to_grey(np.asarray(im))
    assert np.array_equal(io.load(tmp_path / "565.bmp"), expected)


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
