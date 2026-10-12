"""Command line: ``fftvel simulate`` and ``fftvel estimate FILE`` (also ``python -m fftvel``).

Argument parsing and output only: the files are read by :mod:`fftvel.io`,
whose errors become one-line messages that name the options. The result is
written to stdout as JSON, notes on what was read and used to stderr.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import warnings

import numpy as np

from . import __version__
from .core import (DEFAULT_L, MIN_EXPOSURES, MIN_SIZE, SAT_WARN, estimate_velocity, pulsed_from_stack,
                   to_physical)
from .io import _SeveralArrays, frames_axis, load
from .simulate import simulate_pulsed_image

DESCRIPTION = """\
Surface velocity from a single multi-flash image (2-D FFT and Radon
transform). "simulate" estimates a simulated image whose motion is known,
"estimate" an image file or a frame stack. The result is written to stdout
as JSON, notes on what was read and used to stderr.
"""
EXAMPLES = """\
examples:
  %(prog)s estimate image.png --exposures 9
  %(prog)s estimate image.png --exposures 9 --gsd 0.05 --period 0.001
  %(prog)s estimate frames.npy --stack --exposures 9
  %(prog)s simulate --exposures 9 --vx 5 --vy -3 --plot estimate.png
"""
SIMULATE = """\
Simulate a multi-flash image of a random texture that moves by (--vx, --vy)
px per flash, estimate it and write the result as JSON to stdout. The true
direction and speed are written to stderr.
"""
SIMULATE_EXAMPLES = """\
examples:
  %(prog)s --exposures 9 --vx 5 --vy -3
  %(prog)s --exposures 9 --vx 0.5 --vy 0 --noise 0.2 --plot estimate.png
"""
ESTIMATE = """\
Estimate direction and speed from a multi-flash image, or from a stack of
frames averaged into one image (--stack). A single image needs the number
of flashes, --exposures K (or --L L for K = 2L+1): a wrong count scales the
speed by K_true/K, and the result can still be marked valid.

Files: PNG, TIFF, BMP, JPEG and other images (Pillow; colour becomes grey,
several pages a stack), .npy, .npz and .mat (v7.3 needs h5py). Use a
lossless format for measurements: the 8 x 8 px blocks of JPEG can give
wrong results that are marked valid. What was read and used is reported on
stderr, the result as JSON on stdout, with the speed in px per flash (px
per frame for --stack).
"""
ESTIMATE_EXAMPLES = """\
examples:
  %(prog)s image.png --exposures 9
  %(prog)s data/measurements/synthetic_flash/flash_b.png --exposures 8 --duty 0.2 --gsd 0.05 --period 0.0005
  %(prog)s data/recordings/vid2.mat --stack --exposures 9
  %(prog)s burst_0001.npz --stack
"""


class _Formatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    """Descriptions as written, and '(default: X)' for options whose default is not None or False."""

    def _get_help_string(self, action):
        if action.default is None or action.default is False:
            return action.help
        return super()._get_help_string(action)


class _Notes(logging.Handler):
    """Writes the records of the logger fftvel.io (what was read) to stderr as notes."""

    def emit(self, record):
        _note(record.getMessage())


def _note(msg: str) -> None:
    """Report what was read and used on stderr, so stdout stays plain JSON."""
    print(msg, file=sys.stderr)


def _size(shape) -> str:
    return " x ".join(str(d) for d in shape) if len(shape) else "scalar"


def _prog(argv0: str | None) -> str:
    """Name of the program in the help: the console script 'fftvel', else 'python -m fftvel'."""
    name = os.path.basename(argv0 or "")
    stem = name[:-4] if name.lower().endswith(".exe") else name
    return stem if stem.startswith("fftvel") else "python -m fftvel"


def _read(path: str, var: str | None, stack: bool) -> np.ndarray:
    """:func:`fftvel.io.load` with its notes and warnings on stderr and its errors as messages."""
    if var is not None and not path.lower().endswith((".mat", ".npz")):
        raise SystemExit("--var applies to .mat and .npz files only")
    log = logging.getLogger("fftvel.io")
    handler, level, propagate = _Notes(), log.level, log.propagate
    log.addHandler(handler)
    log.setLevel(logging.INFO)
    log.propagate = False
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", UserWarning)
            a = load(path, var=var, stack=stack)
    except _SeveralArrays as err:
        raise SystemExit(f"{err.path} holds several numeric arrays; choose one with --var:\n{err.listing}") from None
    except (OSError, ValueError, ImportError) as err:
        raise SystemExit(str(err)) from None
    finally:
        log.removeHandler(handler)
        log.setLevel(level)
        log.propagate = propagate
    for w in caught:
        _note(f"warning: {w.message}")
    return a


def _pulsed_image(a: np.ndarray, stack: bool, L: int | None = None, exposures: int | None = None,
                  axis_choice: str = "auto"):
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
                             f"(a colour image in a .npy, .npz or .mat file must be converted to grey first)")
        if L is None and exposures is None:
            raise SystemExit("pass --exposures or --L: the number of exposures (flashes) in the image, "
                             "K or 2L+1. A wrong count scales the speed by K_true/K and the result "
                             "can still be marked valid")
        return a, (2 * L + 1 if exposures is None else exposures)
    if a.ndim != 3:
        raise SystemExit(f"--stack expects a 3-D array (rows, cols, frames), got shape {shape}")
    try:
        axis = frames_axis(a.shape, axis_choice)
    except ValueError as err:
        raise SystemExit(f"{err}; pass --frames-axis first or --frames-axis last") from None
    where = "first" if axis == 0 else "last"
    if L is None and exposures is None and a.shape[axis] % 2 == 0:
        raise SystemExit(f"the stack has {a.shape[axis]} frames, an even number; pass --exposures K to "
                         f"average K of them, or --L for 2L+1 (stack of {_size(shape)}, frames on the "
                         f"{where} axis)")
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


def _parser():
    """The parser and its subparsers ``simulate`` and ``estimate``."""
    p = argparse.ArgumentParser(prog=_prog(sys.argv[0] if sys.argv else None), description=DESCRIPTION,
                                epilog=EXAMPLES, formatter_class=_Formatter)
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True, metavar="{simulate,estimate}")

    s = sub.add_parser("simulate", help="simulate a multi-flash image with known motion and estimate it",
                       description=SIMULATE, epilog=SIMULATE_EXAMPLES, formatter_class=_Formatter)
    s.add_argument("--exposures", type=int, metavar="K",
                   help="number K >= 2 of flashes in the image (even or odd); instead of --L")
    s.add_argument("--L", type=int,
                   help=f"half-count of the thesis, K = 2L+1 (default {DEFAULT_L}, so K = {2 * DEFAULT_L + 1})")
    s.add_argument("--vx", type=float, default=3.0,
                   help="velocity along x (the columns, to the right) in px per flash")
    s.add_argument("--vy", type=float, default=-2.0,
                   help="velocity along y (the rows, pointing down) in px per flash")
    s.add_argument("--n", type=int, default=256, help="side of the square image in px")
    s.add_argument("--noise", type=float, default=0.0,
                   help="standard deviation of white noise added to the image, relative to the standard "
                        "deviation of the image without noise")
    s.add_argument("--seed", type=int, default=0, help="seed of the random texture and of the noise")
    s.add_argument("--plot", metavar="PNG", help="save the diagnostic figure as PNG")

    e = sub.add_parser("estimate", help="estimate the motion in an image or frame stack file",
                       description=ESTIMATE, epilog=ESTIMATE_EXAMPLES, formatter_class=_Formatter)
    e.add_argument("file", help="image, .npy, .npz or .mat file")
    e.add_argument("--exposures", type=int, metavar="K",
                   help="number K >= 2 of flashes in the image (even or odd), or of frames to average "
                        "with --stack. A single image needs --exposures or --L; for --stack the default is "
                        "all frames")
    e.add_argument("--L", type=int, help="half-count of the thesis, K = 2L+1; instead of --exposures")
    e.add_argument("--duty", type=float, default=0.0,
                   help="flash (or exposure) duration divided by the time step, in [0, 1)")
    e.add_argument("--roi", type=int, default=512,
                   help="side of the analysed central square in px, 0 for the largest square; the run "
                        "time grows with its square")
    e.add_argument("--clip-outliers", action="store_true",
                   help="set pixels more than 6 robust sigma from the median, and all pixels within "
                        "3 px of them, to the median (static glints, hot pixels)")
    e.add_argument("--stack", action="store_true",
                   help="the file holds frames; average the first K (or 2L+1) of them")
    e.add_argument("--frames-axis", choices=("auto", "first", "last"), default="auto",
                   help="axis of the stack that holds the frames, after dropping axes of length 1; auto "
                        "takes the first or last axis if it is clearly the shortest")
    e.add_argument("--gsd", type=float, metavar="LENGTH_PER_PX",
                   help="ground sample distance: the length on the surface that one pixel covers, in any "
                        "unit (for example mm per px). With --period the JSON also holds speed_per_s, "
                        "vx_per_s and vy_per_s in that unit per second (mm/s for mm per px)")
    e.add_argument("--period", type=float, metavar="S",
                   help="time step in seconds: between two flashes, or between two frames for --stack. "
                        "Goes with --gsd")
    e.add_argument("--var", help="array to read in a .mat or .npz file (needed if it holds several arrays "
                                 "and none is the only image, or with --stack the only stack; in a .npz "
                                 "file the array frames is read if there is one)")
    e.add_argument("--plot", metavar="PNG", help="save the diagnostic figure as PNG")
    return p, s, e


def main(argv=None) -> int:
    p, s, e = _parser()
    a = p.parse_args(argv)
    sp = s if a.cmd == "simulate" else e
    if a.L is not None and a.L < 1:
        sp.error("--L must be a positive integer")
    if a.exposures is not None and a.exposures < MIN_EXPOSURES:
        sp.error(f"--exposures must be at least {MIN_EXPOSURES}")
    if a.cmd == "simulate":
        if a.L is not None and a.exposures is not None and a.exposures != 2 * a.L + 1:
            s.error(f"--L {a.L} means {2 * a.L + 1} flashes but --exposures is {a.exposures}; "
                    f"give only one of them")
        if a.n < MIN_SIZE:
            s.error(f"--n must be at least {MIN_SIZE}")
        if not (np.isfinite(a.noise) and a.noise >= 0):
            s.error("--noise must be a finite number >= 0")
    else:
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
        if a.L is not None:
            L = a.L
        else:
            L = DEFAULT_L if a.exposures is None else max(1, (a.exposures - 1) // 2)
        K = 2 * L + 1 if a.exposures is None else a.exposures
        bp = simulate_pulsed_image(a.n, (a.vx, a.vy), L, exposures=K, noise=a.noise, rng=a.seed)
        true_angle = float(np.rad2deg(np.arctan2(-a.vy, a.vx)) % 180)
        _note(f"true: angle {true_angle:.2f} deg, speed {np.hypot(a.vx, a.vy):.3f} px per flash")
        est = estimate_velocity(bp, L=L, exposures=K, diagnostics=bool(a.plot))
    else:
        bp, K = _pulsed_image(_read(a.file, a.var, a.stack), a.stack, a.L, a.exposures, a.frames_axis)
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
