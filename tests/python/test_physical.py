"""fftvel.to_physical and the --gsd/--period options of the command line.

tests/matlab/test_physical.m checks fftvel.to_physical in MATLAB with the
same numbers.
"""

import json
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.io import loadmat

from fftvel import Estimate, estimate_velocity, simulate_pulsed_image, to_physical
from fftvel.__main__ import main

# input and result shared with tests/matlab/test_physical.m
EST = dict(speed=3.6053444133199184, vx=2.997450239305243, vy=-2.0034471297107115,
           angle_deg=33.758125216267615)
GSD, PERIOD = 0.0123, 1 / 750          # mm per px, s
EXPECTED = dict(speed=33.25930221287625, vx=27.65147845759087, vy=-18.481799771581315)


def test_known_values():
    p = to_physical(SimpleNamespace(**EST), GSD, PERIOD)
    for k, v in EXPECTED.items():
        assert getattr(p, k) == pytest.approx(v, rel=1e-12), k
    assert p.angle_deg == EST["angle_deg"]
    # 2.5 px per flash step, 0.04 mm per px, 1 ms between flashes: 100 mm/s
    e = Estimate(90.0, 2.5, 0.0, -2.5, 3.0, 0.99, 0.5, True, False, 0.7, True, 4, 9, 256, 512)
    p = to_physical(e, gsd=0.04, period=1e-3)
    assert p.speed == pytest.approx(100, rel=1e-12) and p.vy == pytest.approx(-100, rel=1e-12)
    assert p.vx == 0 and p.angle_deg == 90
    # proportional to gsd, inversely to period; also for an estimate that is not valid
    q = to_physical(replace(e, valid=False), gsd=0.08, period=0.5e-3)
    assert q.speed == pytest.approx(4 * p.speed, rel=1e-12)


def test_estimate_end_to_end(repo):
    # parity cases: estimate as in test_parity.py, convert, compare with the expected speed and angle
    cases = loadmat(repo / "tests" / "parity" / "cases.mat")
    ex = np.loadtxt(repo / "tests" / "parity" / "expected.csv", delimiter=",", skiprows=1)
    for row in ex[[0, 5, 9]]:            # cases 1, 6 (odd n, pad 3) and 10 (duty)
        i = int(row[0])
        e = estimate_velocity(cases[f"bp{i:02d}"].astype(float), L=int(row[1]), pad_factor=int(row[2]),
                              duty=row[9])
        p = to_physical(e, GSD, PERIOD)
        speed = row[4] * GSD / PERIOD
        a = np.deg2rad(row[3])
        assert abs(p.speed / speed - 1) < 1e-9, i
        assert abs(p.vx - speed * np.cos(a)) < 1e-9 * speed, i
        assert abs(p.vy + speed * np.sin(a)) < 1e-9 * speed, i
        assert abs(p.angle_deg - row[3]) < 1e-9, i


@pytest.mark.parametrize("name, value", [("gsd", 0), ("gsd", -0.1), ("gsd", np.nan), ("gsd", np.inf),
                                         ("gsd", "0.1"), ("gsd", None), ("gsd", True), ("gsd", [0.1]),
                                         ("period", 0), ("period", -1e-3), ("period", np.inf),
                                         ("period", "1e-3")])
def test_invalid_arguments(name, value):
    args = {"gsd": GSD, "period": PERIOD, name: value}
    with pytest.raises(ValueError, match=name):
        to_physical(SimpleNamespace(**EST), **args)


def test_needs_an_estimate():
    with pytest.raises(TypeError, match="vx, vy"):
        to_physical(SimpleNamespace(speed=1.0, angle_deg=0.0), GSD, PERIOD)


# command line

@pytest.fixture(scope="module")
def npy(tmp_path_factory):
    path = tmp_path_factory.mktemp("phys") / "bp.npy"
    np.save(path, simulate_pulsed_image(128, (3.0, -2.0), 4, rng=1))
    return path


def run(argv, capsys):
    assert main([str(x) for x in argv]) == 0
    cap = capsys.readouterr()
    return json.loads(cap.out), cap.err


def fails(argv, capsys):
    with pytest.raises(SystemExit):
        main([str(x) for x in argv])
    return capsys.readouterr().err


def test_cli_gsd_and_period(npy, capsys):
    r, err = run(["estimate", npy, "--L", "4", "--gsd", "0.05", "--period", "0.002"], capsys)
    p = to_physical(estimate_velocity(np.load(npy), L=4), 0.05, 0.002)
    assert r["gsd"] == 0.05 and r["period"] == 0.002
    for key, value in (("speed_per_s", p.speed), ("vx_per_s", p.vx), ("vy_per_s", p.vy)):
        assert r[key] == float(f"{value:.6g}"), key
    assert abs(r["speed_per_s"] / (r["speed"] * 0.05 / 0.002) - 1) < 1e-4
    assert "speed_per_s = speed * gsd / period" in err
    r, _ = run(["estimate", npy, "--L", "4"], capsys)
    assert not {"speed_per_s", "vx_per_s", "vy_per_s", "gsd", "period"} & set(r)


@pytest.mark.parametrize("extra, message", [
    (["--gsd", "0.05"], "go together"),
    (["--period", "0.002"], "go together"),
    (["--gsd", "0", "--period", "0.002"], "--gsd must be"),
    (["--gsd", "0.05", "--period", "-1"], "--period must be"),
    (["--gsd", "nan", "--period", "0.002"], "--gsd must be"),
])
def test_cli_rejects_bad_units(npy, capsys, extra, message):
    assert message in fails(["estimate", npy, "--L", "4", *extra], capsys)
