"""Same input as tests/matlab/test_parity.m, same expected results."""

import numpy as np
from scipy.io import loadmat

from fftvel import estimate_velocity

from conftest import angle_error


def test_parity(repo):
    pdir = repo / "tests" / "parity"
    cases = loadmat(pdir / "cases.mat")
    ex = np.loadtxt(pdir / "expected.csv", delimiter=",", skiprows=1)
    for row in ex:
        i, L, pad = int(row[0]), int(row[1]), int(row[2])
        e = estimate_velocity(cases[f"bp{i:02d}"].astype(float), L=L, pad_factor=pad)
        assert abs(angle_error(e.angle_deg, row[3])) < 1e-9, i
        assert abs(e.speed / row[4] - 1) < 1e-9, i
        assert abs(e.geba / row[5] - 1) < 1e-9, i
        assert abs(e.quality - row[6]) < 1e-9, i
        assert e.valid == bool(row[7]), i
