"""Same input as tests/matlab/test_parity.m, same expected results."""

import numpy as np
from scipy.io import loadmat

from fftvel import angle_difference, estimate_velocity, pulsed_from_stack


def test_parity(repo):
    pdir = repo / "tests" / "parity"
    cases = loadmat(pdir / "cases.mat")
    ex = np.loadtxt(pdir / "expected.csv", delimiter=",", skiprows=1)
    assert ex.shape == (15, 16)
    for row in ex:
        i, L, pad = int(row[0]), int(row[1]), int(row[2])
        exposures = int(row[8]) or None
        roi = int(row[10]) or None
        x = cases[f"bp{i:02d}"]
        if row[12]:
            bp, K = pulsed_from_stack(x, L=L, exposures=exposures)
            assert K == (exposures or 2 * L + 1), i
        else:
            bp = x.astype(float)
        e = estimate_velocity(bp, L=L, pad_factor=pad, exposures=exposures, duty=row[9], roi=roi,
                              clip_outliers=bool(row[11]))
        assert abs(angle_difference(e.angle_deg, row[3])) < 1e-9, i
        assert abs(e.speed / row[4] - 1) < 1e-9, i
        assert abs(e.geba / row[5] - 1) < 1e-9, i
        assert abs(e.quality - row[6]) < 1e-9, i
        assert e.valid == bool(row[7]), i
        assert abs(e.threshold - row[13]) < 1e-12, i
        assert abs(e.quality_perp - row[14]) < 1e-9, i
        assert e.isotropic == bool(row[15]), i
