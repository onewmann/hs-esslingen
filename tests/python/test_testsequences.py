"""Ground truth from data/testsequences.mat (see tests/matlab/test_testsequences.m)."""

import numpy as np
import pytest
from scipy.io import loadmat

from fftvel import estimate_velocity

from conftest import angle_error


@pytest.mark.parametrize("name, step", [("model", 5), ("tyre", 4)])
@pytest.mark.parametrize("start", [0, 5, 10])
def test_constant_motion(repo, name, step, start):
    s = loadmat(repo / "data" / "testsequences.mat")[name].astype(float) / 255
    last = start + 8
    b = step * last                      # zero-filled border of the last frame
    bp = s[b:, b:, start:last + 1].mean(axis=2)
    e = estimate_velocity(bp, L=4)
    assert abs(angle_error(e.angle_deg, 135)) < 0.1
    assert abs(e.speed / (step * np.sqrt(2)) - 1) < 0.005
    assert e.valid
