"""Command line: ``python -m fftvel simulate`` or ``python -m fftvel estimate FILE``."""

from __future__ import annotations

import argparse
import json
import sys

import numpy as np

from .core import estimate_velocity
from .simulate import simulate_pulsed_image


def _load(path: str, var: str | None) -> np.ndarray:
    """Load a 2-D image or a 3-D array (.npy, .mat) as float."""
    low = path.lower()
    if low.endswith(".npy"):
        return np.asarray(np.load(path), dtype=float)
    if low.endswith(".mat"):
        try:
            from scipy.io import loadmat
        except ImportError:
            raise SystemExit('reading .mat files needs scipy: pip install -e "python[data]"') from None

        data = loadmat(path)
        names = [k for k in data if not k.startswith("__")]
        return np.asarray(np.squeeze(data[var or names[0]]), dtype=float)
    import matplotlib.image as mpimg

    a = np.asarray(mpimg.imread(path), dtype=float)
    if a.ndim == 3:                # RGB(A) -> grey, MATLAB rgb2gray weights
        a = a[..., :3] @ np.array([0.2989, 0.5870, 0.1140])
    return a


def _pulsed_image(a: np.ndarray, stack: bool, L: int | None):
    """Return (bp, L). A stack holds frames along its last axis."""
    if not stack:
        if a.ndim != 2:
            raise SystemExit(f"expected a 2-D image, got shape {a.shape}; use --stack for a frame stack")
        return a, (4 if L is None else L)
    if a.ndim != 3:
        raise SystemExit(f"--stack expects a 3-D array (rows, cols, frames), got shape {a.shape}")
    k = a.shape[2]
    if L is None:
        if k % 2 == 0:
            raise SystemExit(f"the stack has {k} frames, an even number; pass --L to choose 2L+1 of them")
        L = (k - 1) // 2
    if 2 * L + 1 > k:
        raise SystemExit(f"--L {L} needs {2 * L + 1} frames, the stack has {k}")
    return a[:, :, :2 * L + 1].mean(axis=2), L


def _summary(est) -> dict:
    return {k: (round(v, 4) if isinstance(v, float) else v)
            for k, v in vars(est).items() if k != "diag"}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="python -m fftvel", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("simulate", help="simulate a pulsed image with known v and estimate it")
    s.add_argument("--vx", type=float, default=3.0)
    s.add_argument("--vy", type=float, default=-2.0)
    s.add_argument("--L", type=int, default=4)
    s.add_argument("--n", type=int, default=256)
    s.add_argument("--noise", type=float, default=0.0)
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--plot", metavar="PNG", help="save the diagnostic figure")

    e = sub.add_parser("estimate", help="estimate v from an image (.png/.jpg/.npy/.mat)")
    e.add_argument("file")
    e.add_argument("--L", type=int, help="half number of exposures, 2L+1 in total "
                                         "(default 4 for an image; for --stack the stack depth decides)")
    e.add_argument("--stack", action="store_true",
                   help="the file holds frames along its last axis; average 2L+1 of them")
    e.add_argument("--var", help="variable name inside a .mat file")
    e.add_argument("--plot", metavar="PNG", help="save the diagnostic figure")

    a = p.parse_args(argv)
    if a.cmd == "simulate":
        L = a.L
        bp = simulate_pulsed_image(a.n, (a.vx, a.vy), L, noise=a.noise, rng=a.seed)
        true_angle = float(np.rad2deg(np.arctan2(-a.vy, a.vx)) % 180)
        print(f"true: angle {true_angle:.2f} deg, speed {np.hypot(a.vx, a.vy):.3f} px/step")
    else:
        bp, L = _pulsed_image(_load(a.file, a.var), a.stack, a.L)
    est = estimate_velocity(bp, L=L, diagnostics=bool(a.plot))
    print(json.dumps(_summary(est), indent=2))
    if a.plot:
        import matplotlib

        matplotlib.use("Agg")
        from .plotting import plot_estimate

        plot_estimate(bp, est).savefig(a.plot, dpi=120)
        print(f"figure written to {a.plot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
