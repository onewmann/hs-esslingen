"""Command line: ``python -m fftvel simulate`` or ``python -m fftvel estimate FILE``."""

from __future__ import annotations

import argparse
import json
import sys

import numpy as np

from .core import estimate_velocity
from .simulate import simulate_pulsed_image


def _load_image(path: str, var: str | None) -> np.ndarray:
    if path.endswith(".npy"):
        return np.load(path)
    if path.endswith(".mat"):
        from scipy.io import loadmat

        data = loadmat(path)
        names = [k for k in data if not k.startswith("__")]
        a = np.asarray(data[var or names[0]], dtype=float)
        if a.ndim == 3:            # a frame stack: average it into a pulsed image
            a = a.mean(axis=2)
        return a
    import matplotlib.image as mpimg

    a = np.asarray(mpimg.imread(path), dtype=float)
    if a.ndim == 3:                # RGB(A) -> grey, MATLAB rgb2gray weights
        a = a[..., :3] @ np.array([0.2989, 0.5870, 0.1140])
    return a


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
    e.add_argument("--L", type=int, default=4, help="half number of exposures (2L+1 in total)")
    e.add_argument("--var", help="variable name inside a .mat file")
    e.add_argument("--plot", metavar="PNG", help="save the diagnostic figure")

    a = p.parse_args(argv)
    if a.cmd == "simulate":
        bp = simulate_pulsed_image(a.n, (a.vx, a.vy), a.L, noise=a.noise, rng=a.seed)
        true_angle = float(np.rad2deg(np.arctan2(-a.vy, a.vx)) % 180)
        print(f"true: angle {true_angle:.2f} deg, speed {np.hypot(a.vx, a.vy):.3f} px/step")
    else:
        bp = _load_image(a.file, a.var)
    est = estimate_velocity(bp, L=a.L, diagnostics=bool(a.plot))
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
