"""Command line: ``python -m fftvel simulate`` or ``python -m fftvel estimate FILE``."""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

from .core import (DEFAULT_L, GREY_WEIGHTS, MIN_EXPOSURES, MIN_SIZE, SAT_WARN, estimate_velocity,
                   pulsed_from_stack, to_physical)
from .simulate import simulate_pulsed_image

NUMERIC = ("double", "single", "int8", "uint8", "int16", "uint16", "int32", "uint32",
           "int64", "uint64", "logical")          # MATLAB classes read as numeric arrays
HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"             # at byte 0, 512, 1024 or 2048 of an HDF5 file
PIL_KEEP = ("1", "L", "I", "F", "RGB", "RGBA")     # Pillow modes read as they are (also I;16*)
GREY = np.array(GREY_WEIGHTS)                     # MATLAB rgb2gray weights


def _note(msg: str) -> None:
    """Report what was read and used on stderr, so stdout stays plain JSON."""
    print(msg, file=sys.stderr)


def _size(shape) -> str:
    return " x ".join(str(d) for d in shape) if len(shape) else "scalar"


def _grey(a: np.ndarray) -> np.ndarray:
    """RGB(A) to grey with the MATLAB rgb2gray weights; other arrays unchanged."""
    if a.ndim == 3 and a.shape[2] in (3, 4):
        return a[..., :3] @ GREY
    return a


def _load(path: str, var: str | None = None, stack: bool = False) -> np.ndarray:
    """Read an image or array file in its own data type (.npy, .mat or an image file).

    .npy files are memory-mapped, so only the frames used are read from disk.
    Colour images become grey, multi-page images a stack (rows, cols, pages).
    ``stack`` tells which variable of a .mat file fits (see ``_pick_variable``).
    """
    low = path.lower()
    if not os.path.isfile(path):
        raise SystemExit(f"no such file: {path}")
    if var is not None and not low.endswith(".mat"):
        raise SystemExit("--var applies to .mat files only")
    if low.endswith(".npy"):
        try:
            a = np.load(path, mmap_mode="r")
        except ValueError as err:               # e.g. an array of Python objects
            raise SystemExit(f"cannot read {path}: {err}") from None
        _note(f"{path}: {_size(a.shape)} {a.dtype}")
    elif low.endswith(".mat"):
        a = _load_mat(path, var, stack)
    else:
        a = _load_image(path)
    if a.dtype.kind not in "biuf":
        raise SystemExit(f"{path} holds {a.dtype} data, not a real numeric array")
    return a


def _pick_variable(path: str, info, var: str | None, stack: bool = False) -> str:
    """Name of the variable to read.

    That is ``var``, else the only numeric array, else the only one with 2
    (or with ``stack`` 3) dimensions longer than 1, else the only one with at
    least 2. ``info`` lists ``(name, shape, class, numeric)`` for every
    variable in the file.
    """
    numeric = [(n, s, c) for n, s, c, ok in info if ok]
    width = max([len(n) for n, s, c in numeric] + [1])
    listing = "\n".join(f"  {n:<{width}}  {_size(s)} {c}" for n, s, c in numeric) or "  (none)"
    hint = "frames = s.frames (struct) or frames = cat(3, c{:}) (cell); save(file, 'frames', '-v7')"
    if var is not None:
        found = [(c, ok) for n, s, c, ok in info if n == var]
        if not found:
            raise SystemExit(f"{path} has no variable '{var}'. Numeric arrays in the file:\n{listing}")
        if not found[0][1]:
            raise SystemExit(f"variable '{var}' in {path} is of class {found[0][0]}, not a numeric "
                             f"array. Save the array itself, e.g. {hint}. Numeric arrays in the file:\n"
                             f"{listing}")
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
        raise SystemExit(f"{path} holds no numeric array, only: {others}. Save the array itself, "
                         f"e.g. {hint}")
    raise SystemExit(f"{path} holds several numeric arrays; choose one with --var:\n{listing}")


def _load_mat(path: str, var: str | None, stack: bool = False) -> np.ndarray:
    """One numeric variable of a MAT-file (v5/v7 through scipy, v7.3 through h5py)."""
    with open(path, "rb") as f:
        head = f.read(2056)
    if any(head[o:o + 8] == HDF5_SIGNATURE for o in (0, 512, 1024, 2048)):
        return _load_mat73(path, var, stack)
    if head.startswith((b"# Created by Octave", b"Octave-1-")):
        raise SystemExit(f"{path} is in Octave's own format (its save writes text by default), "
                         f"not a MAT-file; save it with save('-v7', file, 'frames') instead")
    try:
        from scipy.io import loadmat, whosmat
    except ImportError:
        raise SystemExit('reading .mat files needs scipy: pip install -e "./python[data]"') from None
    try:
        info = [(n, s, c, c in NUMERIC) for n, s, c in whosmat(path)]
    except NotImplementedError:                 # v7.3 that the signature check missed
        return _load_mat73(path, var, stack)
    except Exception as err:                    # scipy raises ValueError, TypeError or MatReadError
        raise SystemExit(f"cannot read {path} as a MAT-file ({err}); save it with "
                         f"save(file, 'frames', '-v7') in MATLAB or save('-v7', file, 'frames') "
                         f"in Octave") from None
    name = _pick_variable(path, info, var, stack)
    try:
        a = loadmat(path, variable_names=[name])[name]
    except Exception as err:                    # e.g. a damaged or truncated file
        raise SystemExit(f"cannot read variable '{name}' of {path} ({err})") from None
    _note(f"{path}: variable '{name}', {_size(a.shape)} {a.dtype}")
    return a


def _load_mat73(path: str, var: str | None, stack: bool = False) -> np.ndarray:
    """One numeric variable of a MATLAB v7.3 (HDF5) file.

    MATLAB stores arrays column-major, so h5py returns them with the axes in
    reverse order; they are transposed back here.
    """
    try:
        import h5py
    except ImportError:
        raise SystemExit(f"{path} is a MATLAB v7.3 (HDF5) file; reading it needs h5py: "
                         f'pip install -e "./python[data]" (or pip install h5py). Or save it in '
                         f"MATLAB with save(file, 'frames', '-v7')") from None
    with h5py.File(path, "r") as f:
        info = []
        for name, obj in f.items():
            if name.startswith("#"):            # #refs#, #subsystem#: data behind cells and objects
                continue
            cls = obj.attrs.get("MATLAB_class", b"")
            cls = cls.decode() if isinstance(cls, bytes) else str(cls)
            if isinstance(obj, h5py.Group):     # struct, object or sparse matrix
                info.append((name, (), "sparse" if "MATLAB_sparse" in obj.attrs else cls or "group", False))
                continue
            numeric = obj.dtype.kind in "biuf" and "MATLAB_empty" not in obj.attrs
            if cls and obj.dtype.kind == "V":
                cls = "complex " + cls          # stored as (real, imag) pairs
            if cls:
                info.append((name, obj.shape[::-1], cls, numeric and cls in NUMERIC))
            else:                               # HDF5 data not written by MATLAB: keep the axes
                info.append((name, obj.shape, obj.dtype.name if numeric else "non-numeric", numeric))
        name = _pick_variable(path, info, var, stack)
        a = np.asarray(f[name][()])
        if "MATLAB_class" in f[name].attrs:
            a = a.T
    _note(f"{path}: variable '{name}' (MATLAB v7.3), {_size(a.shape)} {a.dtype}")
    return a


def _load_image(path: str) -> np.ndarray:
    """Image file in its own bit depth (8/16 bit, int32, float32), through Pillow."""
    try:
        from PIL import Image, ImageSequence
    except ImportError:
        return _load_image_mpl(path)
    try:
        with Image.open(path) as im:
            fmt = im.format
            pages = []
            for page in ImageSequence.Iterator(im):
                if page.mode not in PIL_KEEP and not page.mode.startswith("I;16"):
                    page = page.convert("L" if page.mode in ("LA", "La") else "RGB")
                pages.append(_grey(np.asarray(page)))
    except (OSError, ValueError, SyntaxError) as err:     # UnidentifiedImageError is an OSError
        if path.lower().endswith((".tif", ".tiff")):
            return _load_tiff(path, err)
        raise SystemExit(f"cannot read {path} as an image ({err})") from None
    if fmt in ("JPEG", "MPO"):
        _warn_jpeg()
    return _stack_pages(path, pages, fmt)


def _load_tiff(path: str, err: Exception) -> np.ndarray:
    """TIFF that Pillow cannot read (e.g. float64), through tifffile if it is installed."""
    try:
        import tifffile
    except ImportError:
        raise SystemExit(f"cannot read {path} with Pillow ({err}). Save it as 8 or 16 bit integer or "
                         f"float32 TIFF, or as .npy, or pip install tifffile") from None
    try:
        with tifffile.TiffFile(path) as tf:
            pages = [_grey(p.asarray()) for p in tf.pages]
    except (OSError, ValueError) as err2:
        raise SystemExit(f"cannot read {path} as a TIFF ({err2})") from None
    return _stack_pages(path, pages, "TIFF (tifffile)")


def _load_image_mpl(path: str) -> np.ndarray:
    """Fallback without Pillow: matplotlib.image.imread (PNG scaled to [0, 1])."""
    try:
        import matplotlib.image as mpimg
    except ImportError:
        raise SystemExit(f"reading {path} needs Pillow: pip install pillow (it also comes with "
                         f'pip install -e "./python[plot]"). Without it, save the image as .npy') from None
    if path.lower().endswith((".jpg", ".jpeg", ".jfif")):
        _warn_jpeg()
    return _stack_pages(path, [_grey(np.asarray(mpimg.imread(path)))], "image (matplotlib)")


def _stack_pages(path: str, pages, fmt) -> np.ndarray:
    """One page as it is, several pages of the size of the first one as a stack (rows, cols, pages)."""
    same = [pg for pg in pages if pg.shape == pages[0].shape]
    if len(same) < len(pages):                  # e.g. a thumbnail stored as a second page
        _note(f"{path}: skipped {len(pages) - len(same)} pages of another size than the first")
    if len(same) == 1:
        a = same[0]
        _note(f"{path}: {fmt}, {_size(a.shape)} {a.dtype}")
        return a
    a = np.stack(same, axis=2)
    _note(f"{path}: {fmt}, {len(same)} pages stacked as {_size(a.shape)} {a.dtype}")
    return a


def _warn_jpeg() -> None:
    _note("warning: JPEG compression adds an 8 x 8 px block grid, which can give wrong results "
          "that are marked valid; use a lossless format (PNG, TIFF, BMP) for measurements")


def _frames_axis(shape, choice: str) -> int:
    """0 or -1: the axis of a 3-D stack that holds the frames.

    ``auto`` takes the first or last axis if it is at most half as long as
    each of the other two, and gives up otherwise.
    """
    if choice == "first":
        return 0
    if choice == "last":
        return -1
    if 2 * shape[2] <= min(shape[0], shape[1]):
        return -1
    if 2 * shape[0] <= min(shape[1], shape[2]):
        return 0
    raise SystemExit(f"cannot tell which axis of the {_size(shape)} stack holds the frames; "
                     f"pass --frames-axis first or --frames-axis last")


def _pulsed_image(a: np.ndarray, stack: bool, L: int | None = None, exposures: int | None = None,
                  frames_axis: str = "auto"):
    """Return ``(bp, K)``: the pulsed image and its number of exposures.

    A single image needs ``L`` or ``exposures``. A stack is averaged with
    :func:`fftvel.pulsed_from_stack` over its first ``K`` frames (``exposures``,
    else ``2L+1``, else all frames, which must then be odd). Axes of length 1
    are dropped first, so a MATLAB ``getdata`` stack (rows x cols x 1 x frames)
    works as it is. Only the frames used are converted to float.
    """
    if L is not None and exposures is not None and exposures != 2 * L + 1:
        raise SystemExit(f"--L {L} means {2 * L + 1} exposures but --exposures is {exposures}; "
                         f"give only one of them")
    shape = a.shape
    a = np.squeeze(a)
    if not stack:
        if a.ndim != 2:
            raise SystemExit(f"expected a 2-D image, got shape {shape}; use --stack for a frame stack "
                             f"(a colour image in a .npy or .mat file must be converted to grey first)")
        if L is None and exposures is None:
            raise SystemExit("pass --L or --exposures: the number of exposures (flashes) in the image, "
                             "2L+1 or K. A wrong count scales the speed by K_true/K and the result "
                             "can still be marked valid")
        return a, (2 * L + 1 if exposures is None else exposures)
    if a.ndim != 3:
        raise SystemExit(f"--stack expects a 3-D array (rows, cols, frames), got shape {shape}")
    axis = _frames_axis(a.shape, frames_axis)
    where = "first" if axis == 0 else "last"
    if L is None and exposures is None and a.shape[axis] % 2 == 0:
        raise SystemExit(f"the stack has {a.shape[axis]} frames, an even number; pass --L to average "
                         f"2L+1 of them or --exposures K for K (stack of {_size(shape)}, frames on "
                         f"the {where} axis)")
    try:
        bp, K = pulsed_from_stack(a, L=L, exposures=exposures, axis=axis)
    except ValueError as err:
        raise SystemExit(f"{err} (stack of {_size(shape)}, frames on the {where} axis)") from None
    _note(f"stack: mean of frames 1-{K} of {a.shape[axis]}, frames on the {where} axis of {_size(a.shape)}")
    return bp, K


def _warnings(est) -> None:
    """Warn on stderr about saturated pixels (more than SAT_WARN) and valid but blur-like results."""
    if est.valid and est.blur_like:
        _note(f"warning: one continuous exposure blurred over {est.blur_length:.3g} px (no flashes) fits better "
              f"than the flashes. Below about 2.2 px per flash this is common for multi-flash images with noise "
              f"too, so it only means that the image cannot tell the two apart; if light other than the flashes "
              f"reaches the sensor, check with an image taken without flashes, which must come out not valid")
    if est.saturated > SAT_WARN:
        _note(f"warning: {100 * est.saturated:.1f} % of the analysed pixels are saturated (at the maximum "
              f"value); saturated areas and glints can give wrong results that are marked valid. Lower the "
              f"exposure or the flash energy, or crop them out with --roi")


def _summary(est) -> dict:
    d = {k: (round(v, 4) if isinstance(v, float) else v)
         for k, v in vars(est).items() if k != "diag"}
    if est.exposures != 2 * est.L + 1:
        d["L"] = None              # an even number of exposures has no L
    return d


def _physical(est, gsd: float, period: float) -> dict:
    """``gsd``, ``period`` and the velocity per second (6 significant digits) for the JSON output."""
    p = to_physical(est, gsd, period)
    return {"gsd": gsd, "period": period, "speed_per_s": float(f"{p.speed:.6g}"),
            "vx_per_s": float(f"{p.vx:.6g}"), "vy_per_s": float(f"{p.vy:.6g}")}


def _use_agg() -> None:
    """Check for matplotlib before any work and select its file-only backend (no display needed)."""
    try:
        import matplotlib
    except ImportError:
        raise SystemExit('--plot needs matplotlib: pip install -e "./python[plot]"') from None
    matplotlib.use("Agg")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="python -m fftvel", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("simulate", help="simulate a pulsed image with known v and estimate it")
    s.add_argument("--vx", type=float, default=3.0)
    s.add_argument("--vy", type=float, default=-2.0)
    s.add_argument("--L", type=int, default=DEFAULT_L)
    s.add_argument("--n", type=int, default=256)
    s.add_argument("--noise", type=float, default=0.0)
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--plot", metavar="PNG", help="save the diagnostic figure")

    e = sub.add_parser("estimate", help="estimate v from an image or array file",
                       description="Estimate v from a pulsed image or a frame stack. Files: .png, .tif, "
                                   ".bmp, .jpg and other images (Pillow; colour becomes grey, several pages "
                                   "a stack), .npy, .mat (v7.3 needs h5py). Use a lossless format for "
                                   "measurements: JPEG blocks can give wrong valid results. What was read "
                                   "and used is reported on stderr, the result as JSON on stdout.")
    e.add_argument("file", help="image, .npy or .mat file")
    e.add_argument("--L", type=int, help="half the number of exposures, 2L+1 in total. A single image "
                                         "needs --L or --exposures; for --stack the default is all frames")
    e.add_argument("--exposures", type=int, metavar="K",
                   help="exact number of exposures (flashes), K >= 2, also even; instead of --L")
    e.add_argument("--duty", type=float, default=0.0,
                   help="exposure (or flash) duration divided by the time step, in [0, 1); default 0")
    e.add_argument("--roi", type=int, default=512,
                   help="side of the analysed central square in px (default 512; 0: the largest "
                        "square). The run time grows with its square")
    e.add_argument("--clip-outliers", action="store_true",
                   help="set pixels more than 6 robust sigma from the median, and all pixels within "
                        "3 px of them, to the median (static glints, hot pixels)")
    e.add_argument("--stack", action="store_true",
                   help="the file holds frames; average the first K (or 2L+1) of them")
    e.add_argument("--frames-axis", choices=("auto", "first", "last"), default="auto",
                   help="axis of the stack that holds the frames, after dropping axes of length 1 "
                        "(default auto: the first or last axis if it is clearly the shortest)")
    e.add_argument("--gsd", type=float, metavar="LENGTH_PER_PX",
                   help="ground sample distance: the length on the surface that one pixel covers, in any "
                        "unit (for example mm per px). With --period the JSON also holds speed_per_s, "
                        "vx_per_s and vy_per_s in that unit per second (mm/s for mm per px)")
    e.add_argument("--period", type=float, metavar="S",
                   help="time step in seconds: between two flashes, or between two frames for --stack. "
                        "Goes with --gsd")
    e.add_argument("--var", help="variable name inside a .mat file (needed if it holds several arrays "
                                 "and none is the only image, or with --stack the only stack)")
    e.add_argument("--plot", metavar="PNG", help="save the diagnostic figure")

    a = p.parse_args(argv)
    if a.cmd == "estimate":
        if a.L is not None and a.L < 1:
            e.error("--L must be a positive integer")
        if a.exposures is not None and a.exposures < MIN_EXPOSURES:
            e.error(f"--exposures must be at least {MIN_EXPOSURES}")
        if not 0 <= a.duty < 1:
            e.error("--duty must be in [0, 1)")
        if a.roi < 0 or 0 < a.roi < MIN_SIZE:
            e.error(f"--roi must be 0 (the largest square) or at least {MIN_SIZE}")
        if a.frames_axis != "auto" and not a.stack:
            e.error("--frames-axis needs --stack")
        if (a.gsd is None) != (a.period is None):
            e.error("--gsd and --period go together; speed_per_s needs both")
        for name in ("gsd", "period"):
            value = getattr(a, name)
            if value is not None and not (np.isfinite(value) and value > 0):
                e.error(f"--{name} must be a finite number > 0")
    if a.plot:
        _use_agg()
    if a.cmd == "simulate":
        L = a.L
        bp = simulate_pulsed_image(a.n, (a.vx, a.vy), L, noise=a.noise, rng=a.seed)
        true_angle = float(np.rad2deg(np.arctan2(-a.vy, a.vx)) % 180)
        print(f"true: angle {true_angle:.2f} deg, speed {np.hypot(a.vx, a.vy):.3f} px/step")
        est = estimate_velocity(bp, L=L, diagnostics=bool(a.plot))
    else:
        bp, K = _pulsed_image(_load(a.file, a.var, a.stack), a.stack, a.L, a.exposures, a.frames_axis)
        L = a.L if a.L is not None else max(1, (K - 1) // 2)
        roi = a.roi or None
        h, w = bp.shape
        n = min(h, w) if roi is None else min(roi, h, w)
        _note(f"analysed: central {n} x {n} px of the {h} x {w} image (--roi {a.roi})")
        _note(f"model: {K} exposures" + (f" (L = {L})" if K == 2 * L + 1 else "")
              + f", duty {a.duty:g}" + (", outliers clipped" if a.clip_outliers else ""))
        try:
            est = estimate_velocity(bp, L=L, diagnostics=bool(a.plot), exposures=K, duty=a.duty,
                                    roi=roi, clip_outliers=a.clip_outliers)
        except ValueError as err:
            raise SystemExit(f"{a.file}: {err}") from None
    _warnings(est)
    out = _summary(est)
    if a.cmd == "estimate" and a.gsd is not None:
        _note(f"units: speed_per_s = speed * gsd / period, gsd {a.gsd:g} (length unit) per px, "
              f"period {a.period:g} s")
        out.update(_physical(est, a.gsd, a.period))
    print(json.dumps(out, indent=2))
    if a.plot:
        import matplotlib.pyplot as plt

        from .plotting import plot_estimate

        fig = plot_estimate(bp, est)
        fig.savefig(a.plot, dpi=120)
        plt.close(fig)
        _note(f"figure written to {a.plot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
