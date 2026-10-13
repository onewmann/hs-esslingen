"""Reading images and frame stacks from files.

* :func:`load` reads one image or frame stack from an image file (PNG,
  TIFF, BMP, JPEG and others), a ``.npy``, ``.npz`` or ``.mat`` file, in
  its own data type.
* :func:`frames_axis` tells which axis of a 3-D stack holds the frames.
* :func:`to_grey` converts a colour image to grey with the weights of
  MATLAB ``rgb2gray``.

``.npy`` and ``.npz`` files need only NumPy. Image files need Pillow (or,
without it, matplotlib), float64 and 16-bit colour TIFF files tifffile,
MAT-files SciPy and MATLAB v7.3 files also h5py; the extra ``data`` of the
package brings SciPy, h5py and Pillow.

Errors are raised as FileNotFoundError (no such file), ImportError (a
reader that the file needs is not installed) and ValueError (the file
cannot be read, or holds no fitting numeric array). What was read, such
as the variable, the shape and the data type, is logged at level INFO on
the logger ``fftvel.io``; the command line prints these lines on stderr.
"""

from __future__ import annotations

import contextlib
import logging
import os
import warnings

import numpy as np

from .core import GREY_WEIGHTS

__all__ = ["HDF5_SIGNATURE", "frames_axis", "load", "to_grey"]

HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"             # at byte 0, 512, 1024 or 2048 of an HDF5 file
NUMERIC = ("double", "single", "int8", "uint8", "int16", "uint16", "int32", "uint32",
           "int64", "uint64", "logical")          # MATLAB classes read as numeric arrays
PIL_KEEP = ("1", "L", "I", "F", "RGB", "RGBA")     # Pillow modes read as they are (also I;16*)
MAT_HINT = "frames = s.frames (struct) or frames = cat(3, c{:}) (cell); save(file, 'frames', '-v7')"
NPZ_HINT = "np.savez(file, frames=frames, t=t)"
DEEP_COLOUR = "Pillow reads only 8 bits of each channel of a 16-bit colour (or grey and alpha) image"
JPEG_WARNING = ("JPEG compression adds an 8 x 8 px block grid, which can give wrong results that are "
                "marked valid; use a lossless format (PNG, TIFF, BMP) for measurements")

_log = logging.getLogger(__name__)


class _SeveralArrays(ValueError):
    """Several numeric arrays of a .mat or .npz file fit and none was named.

    ``path`` is the file, ``listing`` lists the numeric arrays, one per
    line. The command line catches it to name its option ``--var``.
    """

    def __init__(self, path: str, listing: str, npz: bool, stack: bool):
        self.path, self.listing = path, listing
        how = "give the name of the one to read, or save only that one"
        if npz:
            how += ", or save it under the name frames"
        super().__init__(f"{path} holds several numeric arrays that could be the "
                         f"{'stack' if stack else 'image'}; {how}:\n{listing}")


def _size(shape) -> str:
    return " x ".join(str(d) for d in shape) if len(shape) else "scalar"


def to_grey(image) -> np.ndarray:
    """Grey-scale float image.

    An RGB or RGBA image (rows x cols x 3 or 4, channels in RGB order)
    becomes ``0.2989 R + 0.5870 G + 0.1140 B``, the weights of MATLAB
    ``rgb2gray``; an alpha channel is ignored. A rows x cols x 1 array
    loses its last axis, and a 2-D array is only converted to float. Other
    shapes raise ValueError.
    """
    a = np.asarray(image, dtype=float)
    if a.ndim == 3 and a.shape[2] == 1:
        return a[..., 0]
    if a.ndim == 3 and a.shape[2] in (3, 4):
        cr, cg, cb = GREY_WEIGHTS
        return cr * a[..., 0] + cg * a[..., 1] + cb * a[..., 2]
    if a.ndim == 2:
        return a
    raise ValueError(f"expected a grey image (rows x cols) or a colour image (rows x cols x 3 or 4), "
                     f"got shape {a.shape}")


def _colour_to_grey(page: np.ndarray) -> np.ndarray:
    """RGB(A) pages to grey; other pages unchanged, in their own data type."""
    return to_grey(page) if page.ndim == 3 and page.shape[2] in (3, 4) else page


def frames_axis(shape, choice: str = "auto") -> int:
    """Axis of a 3-D stack that holds the frames: 0 (the first) or -1 (the last).

    ``shape`` is the shape of the stack. ``choice`` is ``'first'``,
    ``'last'`` or ``'auto'``, which takes the last axis if it is at most
    half as long as each of the other two, else the first axis if that one
    is, and raises ValueError if neither is.
    """
    if choice not in ("auto", "first", "last"):
        raise ValueError(f"choice must be 'auto', 'first' or 'last', got {choice!r}")
    shape = tuple(shape)
    if len(shape) != 3:
        raise ValueError(f"expected the shape of a 3-D stack (rows, cols, frames), got {shape}")
    if choice == "first":
        return 0
    if choice == "last":
        return -1
    if 2 * shape[2] <= min(shape[0], shape[1]):
        return -1
    if 2 * shape[0] <= min(shape[1], shape[2]):
        return 0
    raise ValueError(f"cannot tell which axis of the {_size(shape)} stack holds the frames: neither the "
                     f"first nor the last axis is at most half as long as the other two")


def load(path, *, var: str | None = None, stack: bool = False) -> np.ndarray:
    """Image or frame stack in a file, in its own data type.

    ``path`` is an image file, a ``.npy``, ``.npz`` or ``.mat`` file (the
    extension decides). Image files keep 16-bit, 32-bit integer and float32
    data; colour images become grey (see :func:`to_grey`), several pages of
    the size of the first one a stack (rows, cols, pages). 16-bit colour
    TIFF is read through tifffile; 16-bit colour in other formats raises
    ValueError, because Pillow reads only 8 bits of it. ``.npy`` files
    are memory-mapped, so only the frames used are read from disk. MAT-files
    are read through SciPy (v5 and v7) or h5py (v7.3, axes in MATLAB order).

    ``var`` names the array to read in a ``.mat`` or ``.npz`` file.
    Without it, a ``.npz`` entry named ``frames`` is read if there is one;
    otherwise the only numeric array, else the only one with 2 axes longer
    than 1 (with ``stack`` True: 3), else the only one with at least 2.
    ``stack`` changes nothing else; the frames axis is found with
    :func:`frames_axis`.

    Raises FileNotFoundError if there is no such file, ImportError if a
    reader that the file needs is not installed, and ValueError if the file
    cannot be read or holds no fitting real numeric array. A JPEG file gives
    a UserWarning: its 8 x 8 px blocks can give wrong results that are
    marked valid. What was read is logged at level INFO on the logger
    ``fftvel.io``.
    """
    path = os.fspath(path)
    low = path.lower()
    if not os.path.isfile(path):
        raise FileNotFoundError(f"no such file: {path}")
    if var is not None and not low.endswith((".mat", ".npz")):
        raise ValueError(f"a variable name applies to .mat and .npz files only, not to {path}")
    jpeg = False
    if low.endswith(".npy"):
        try:
            a = np.load(path, mmap_mode="r")
        except Exception as err:  # noqa: BLE001 - objects, an empty or damaged file (EOFError, TokenError, ...)
            raise ValueError(f"cannot read {path} as a .npy file ({err})") from None
        _log.info("%s: %s %s", path, _size(a.shape), a.dtype)
    elif low.endswith(".npz"):
        a = _load_npz(path, var, stack)
    elif low.endswith(".mat"):
        a = _load_mat(path, var, stack)
    else:
        a, jpeg = _load_image(path)
    if a.dtype.kind not in "biuf":
        raise ValueError(f"{path} holds {a.dtype} data, not a real numeric array")
    if jpeg:
        warnings.warn(JPEG_WARNING, UserWarning, stacklevel=2)
    return a


def _pick_variable(path: str, info, var: str | None, stack: bool = False, npz: bool = False) -> str:
    """Name of the variable (or .npz entry) to read.

    That is ``var``, else the only numeric array, else the only one with 2
    (or with ``stack`` 3) dimensions longer than 1, else the only one with
    at least 2. ``info`` lists ``(name, shape, class, numeric)`` for every
    variable in the file.
    """
    numeric = [(n, s, c) for n, s, c, ok in info if ok]
    width = max([len(n) for n, s, c in numeric] + [1])
    listing = "\n".join(f"  {n:<{width}}  {_size(s)} {c}" for n, s, c in numeric) or "  (none)"
    what = "array" if npz else "variable"
    hint = NPZ_HINT if npz else MAT_HINT
    if var is not None:
        found = [(c, ok) for n, s, c, ok in info if n == var]
        if not found:
            raise ValueError(f"{path} has no {what} '{var}'. Numeric arrays in the file:\n{listing}")
        if not found[0][1]:
            raise ValueError(f"{what} '{var}' in {path} is of class {found[0][0]}, not a numeric array. "
                             f"Save the array itself, e.g. {hint}. Numeric arrays in the file:\n{listing}")
        return var
    if len(numeric) == 1:
        return numeric[0][0]
    fits = [n for n, s, c in numeric if sum(d > 1 for d in s) == (3 if stack else 2)]
    if len(fits) == 1:
        return fits[0]
    large = [n for n, s, c in numeric if sum(d > 1 for d in s) >= 2]
    if len(large) == 1:
        return large[0]
    if not numeric:
        others = ", ".join(f"{n} ({c})" for n, s, c, ok in info) or "nothing"
        raise ValueError(f"{path} holds no numeric array, only: {others}. Save the array itself, e.g. {hint}")
    raise _SeveralArrays(path, listing, npz, stack)


def _load_npz(path: str, var: str | None, stack: bool = False) -> np.ndarray:
    """One numeric array of a .npz file: ``var``, else ``frames``, else as in :func:`_pick_variable`."""
    arrays, info = {}, []
    try:
        with np.load(path, allow_pickle=False) as z:
            for name in z.files:
                try:
                    a = z[name]
                except ValueError:              # an array of Python objects, which needs pickle
                    info.append((name, (), "object", False))
                    continue
                arrays[name] = a
                info.append((name, a.shape, a.dtype.name, a.dtype.kind in "biuf"))
    except Exception as err:  # noqa: BLE001 - not a zip file, a damaged one, or a damaged entry
        raise ValueError(f"cannot read {path} as a .npz file ({err})") from None
    if var is None and "frames" in arrays and arrays["frames"].dtype.kind in "biuf":
        name = "frames"
    else:
        name = _pick_variable(path, info, var, stack, npz=True)
    a = arrays[name]
    _log.info("%s: array '%s', %s %s", path, name, _size(a.shape), a.dtype)
    return a


def _load_mat(path: str, var: str | None, stack: bool = False) -> np.ndarray:
    """One numeric variable of a MAT-file (v5/v7 through SciPy, v7.3 through h5py)."""
    with open(path, "rb") as f:
        head = f.read(2056)
    if any(head[o:o + 8] == HDF5_SIGNATURE for o in (0, 512, 1024, 2048)):
        return _load_mat73(path, var, stack)
    if head.startswith((b"# Created by Octave", b"Octave-1-")):
        raise ValueError(f"{path} is in Octave's own format (its save writes text by default), not a "
                         f"MAT-file; save it with save('-v7', file, 'frames') instead")
    try:
        from scipy.io import loadmat, whosmat
    except ImportError:
        raise ImportError('reading .mat files needs SciPy: pip install scipy (the extra "data" of '
                          'fftvel brings it)') from None
    try:
        info = [(n, s, c, c in NUMERIC) for n, s, c in whosmat(path)]
    except NotImplementedError:                 # v7.3 that the signature check missed
        return _load_mat73(path, var, stack)
    except Exception as err:  # noqa: BLE001 - SciPy raises ValueError, TypeError or MatReadError
        raise ValueError(f"cannot read {path} as a MAT-file ({err}); save it with "
                         f"save(file, 'frames', '-v7') in MATLAB or save('-v7', file, 'frames') "
                         f"in Octave") from None
    name = _pick_variable(path, info, var, stack)
    try:
        a = loadmat(path, variable_names=[name])[name]
    except Exception as err:  # noqa: BLE001 - e.g. a damaged or truncated file
        raise ValueError(f"cannot read variable '{name}' of {path} ({err})") from None
    _log.info("%s: variable '%s', %s %s", path, name, _size(a.shape), a.dtype)
    return a


def _load_mat73(path: str, var: str | None, stack: bool = False) -> np.ndarray:
    """One numeric variable of a MATLAB v7.3 (HDF5) file.

    MATLAB stores arrays column-major, so h5py returns them with the axes in
    reverse order; they are transposed back here.
    """
    try:
        import h5py
    except ImportError:
        raise ImportError(f"{path} is a MATLAB v7.3 (HDF5) file; reading it needs h5py: pip install h5py "
                          f"(the extra \"data\" of fftvel brings it). Or save it in MATLAB with "
                          f"save(file, 'frames', '-v7')") from None
    with _reading_mat73(path):
        with h5py.File(path, "r") as f:
            info = _mat73_info(h5py, f)
    name = _pick_variable(path, info, var, stack)   # outside: its ValueError has a message of its own
    with _reading_mat73(path):
        with h5py.File(path, "r") as f:
            a = np.asarray(f[name][()])
            if "MATLAB_class" in f[name].attrs:
                a = a.T
    _log.info("%s: variable '%s' (MATLAB v7.3), %s %s", path, name, _size(a.shape), a.dtype)
    return a


@contextlib.contextmanager
def _reading_mat73(path: str):
    """Turn any error of h5py into ValueError: OSError (truncated file), RuntimeError (damaged file), ..."""
    try:
        yield
    except Exception as err:  # noqa: BLE001 - only h5py calls run inside
        raise ValueError(f"cannot read {path} as a MATLAB v7.3 file ({err})") from None


def _mat73_info(h5py, f) -> list:
    """``(name, shape, class, numeric)`` of every variable in an open MATLAB v7.3 file."""
    info = []
    for name, obj in f.items():
        if name.startswith("#"):                # #refs#, #subsystem#: data behind cells and objects
            continue
        cls = obj.attrs.get("MATLAB_class", b"")
        cls = cls.decode() if isinstance(cls, bytes) else str(cls)
        if isinstance(obj, h5py.Group):         # struct, object or sparse matrix
            info.append((name, (), "sparse" if "MATLAB_sparse" in obj.attrs else cls or "group", False))
            continue
        numeric = obj.dtype.kind in "biuf" and "MATLAB_empty" not in obj.attrs
        if cls and obj.dtype.kind == "V":
            cls = "complex " + cls              # stored as (real, imag) pairs
        if cls:
            info.append((name, obj.shape[::-1], cls, numeric and cls in NUMERIC))
        else:                                   # HDF5 data not written by MATLAB: keep the axes
            info.append((name, obj.shape, obj.dtype.name if numeric else "non-numeric", numeric))
    return info


def _load_image(path: str):
    """``(array, jpeg)``: an image file in its own bit depth (8/16 bit, int32, float32), through Pillow."""
    try:
        from PIL import Image, ImageSequence
    except ImportError:
        return _load_image_mpl(path)
    deep = False
    try:
        with Image.open(path) as im:
            fmt = im.format
            pages = []
            for page in ImageSequence.Iterator(im):
                if _deep_colour(page, fmt):
                    deep = True
                    break
                if page.mode not in PIL_KEEP and not page.mode.startswith("I;16"):
                    page = page.convert("L" if page.mode in ("LA", "La") else "RGB")
                pages.append(_colour_to_grey(np.asarray(page)))
    except (OSError, ValueError, SyntaxError) as err:     # UnidentifiedImageError is an OSError
        if path.lower().endswith((".tif", ".tiff")):
            return _load_tiff(path, err), False
        raise ValueError(f"cannot read {path} as an image ({err})") from None
    if deep and fmt == "TIFF":
        return _load_tiff(path, DEEP_COLOUR), False
    if deep:
        raise ValueError(f"cannot read {path}: {DEEP_COLOUR}; save it as 16-bit grey, as TIFF or as .npy")
    return _stack_pages(path, pages, fmt), fmt in ("JPEG", "MPO")


def _deep_colour(page, fmt) -> bool:
    """True for a colour page with more than 8 bits per sample, which Pillow would read as 8 bits.

    PNG shows it in the raw mode of its tiles (``RGB;16B``; but ``BGR;16``
    of a BMP is 5-6-5 bits in 16), PPM in the maximum value that its tiles
    carry, TIFF in its BitsPerSample tag (a planar TIFF has the raw modes
    ``R``, ``G``, ``B``).
    """
    if page.mode not in ("RGB", "RGBA", "RGBX", "LA", "La"):     # I;16, I and F keep their depth
        return False
    if fmt == "TIFF":
        bits = page.tag_v2.get(258, 8)          # BitsPerSample, one value per channel
        return max(bits if isinstance(bits, tuple) else (bits,)) > 8
    for tile in page.tile:
        args = tile[3] if len(tile) > 3 else None
        raw = args[0] if isinstance(args, tuple) and args else args
        if isinstance(raw, str) and raw.endswith((";16B", ";16L", ";16N")):
            return True
        if tile[0] == "ppm" and isinstance(args, tuple) and len(args) > 1 and args[1] > 255:
            return True                         # Pillow scales the values to 0..255
    return False


def _load_tiff(path: str, err) -> np.ndarray:
    """TIFF that Pillow cannot read (e.g. float64) or not in full (16-bit colour), through tifffile."""
    try:
        import tifffile
    except ImportError:
        raise ImportError(f"cannot read {path} with Pillow ({err}); pip install tifffile, or save it as a "
                          f"grey TIFF of 8 or 16 bit integers or float32, or as .npy") from None
    try:
        with tifffile.TiffFile(path) as tf:
            pages = []
            for p in tf.pages:
                a = p.asarray()
                if a.ndim == 3 and p.axes.startswith("S"):     # planar: samples first (SYX)
                    a = np.moveaxis(a, 0, -1)
                pages.append(_colour_to_grey(a))
    except (OSError, ValueError) as err2:
        raise ValueError(f"cannot read {path} as a TIFF ({err2})") from None
    return _stack_pages(path, pages, "TIFF (tifffile)")


def _load_image_mpl(path: str):
    """Fallback without Pillow: ``matplotlib.image.imread`` (PNG scaled to [0, 1])."""
    try:
        import matplotlib.image as mpimg
    except ImportError:
        raise ImportError(f"reading {path} needs Pillow: pip install pillow (the extras \"data\" and "
                          f"\"plot\" of fftvel bring it). Without it, save the image as .npy") from None
    jpeg = path.lower().endswith((".jpg", ".jpeg", ".jfif"))
    return _stack_pages(path, [_colour_to_grey(np.asarray(mpimg.imread(path)))], "image (matplotlib)"), jpeg


def _stack_pages(path: str, pages, fmt) -> np.ndarray:
    """One page as it is, several pages of the size of the first one as a stack (rows, cols, pages)."""
    same = [pg for pg in pages if pg.shape == pages[0].shape]
    if len(same) < len(pages):                  # e.g. a thumbnail stored as a second page
        _log.info("%s: skipped %d pages of another size than the first", path, len(pages) - len(same))
    if len(same) == 1:
        a = same[0]
        _log.info("%s: %s, %s %s", path, fmt, _size(a.shape), a.dtype)
        return a
    a = np.stack(same, axis=2)
    _log.info("%s: %s, %d pages stacked as %s %s", path, fmt, len(same), _size(a.shape), a.dtype)
    return a
