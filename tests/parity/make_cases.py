"""Write the shared parity inputs and the Python reference results.

    python tests/parity/make_cases.py

creates ``cases.mat`` (input images, float32) and ``expected.csv`` (results of
the Python implementation). ``tests/matlab/test_parity.m`` and
``tests/python/test_parity.py`` check that MATLAB/Octave and Python both
reproduce ``expected.csv`` from ``cases.mat``.
"""

from pathlib import Path

import numpy as np
from scipy.io import savemat

from fftvel import estimate_velocity, simulate_pulsed_image

HERE = Path(__file__).resolve().parent

# (n_rows, n_cols, vx, vy, L, pad_factor, noise)
CASES = [
    (128, 128, 3.0, -2.0, 4, 2, 0.0),
    (128, 128, -1.2, 0.7, 4, 2, 0.0),
    (128, 128, 6.0, 6.0, 4, 2, 0.3),
    (128, 128, 0.4, 2.5, 10, 2, 0.0),
    (150, 190, 4.0, 1.0, 4, 2, 0.0),     # non-square: cropped to 150x150
    (101, 101, 2.0, -3.0, 4, 3, 0.0),    # odd n and odd FFT length
    (128, 128, 0.0, 0.0, 4, 2, 0.0),     # no motion: must come out invalid
    (96, 96, 12.0, -1.0, 3, 2, 0.1),
]


def main():
    rng = np.random.default_rng(20251008)
    mats = {}
    rows = []
    for i, (h, w, vx, vy, L, pad, noise) in enumerate(CASES, start=1):
        n = max(h, w)
        bp = simulate_pulsed_image(n, (vx, vy), L, noise=noise, rng=rng)[:h, :w]
        bp = bp.astype(np.float32)
        mats[f"bp{i:02d}"] = bp
        est = estimate_velocity(bp.astype(float), L=L, pad_factor=pad)
        rows.append([i, L, pad, est.angle_deg, est.speed, est.geba, est.quality, int(est.valid)])
    mats["params"] = np.array([[r[1], r[2]] for r in rows], dtype=float)
    savemat(HERE / "cases.mat", mats, do_compression=True)
    header = "case,L,pad,angle_deg,speed,geba,quality,valid"
    np.savetxt(HERE / "expected.csv", np.array(rows, dtype=float), delimiter=",",
               header=header, comments="", fmt=["%d", "%d", "%d", "%.17g", "%.17g", "%.17g", "%.17g", "%d"])
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
