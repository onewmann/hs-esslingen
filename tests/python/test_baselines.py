"""experiments/baselines.py: the cases, the baselines on a few easy images, the sub-pixel peak, the table,
figure and command line."""

import csv
import math
import multiprocessing
import sys

import numpy as np
import pytest

MODULES = ("baselines", "common", "noise")


@pytest.fixture
def bl(repo, monkeypatch):
    """The module baselines, imported from experiments/ and removed again afterwards."""
    monkeypatch.setattr(sys, "path", [str(repo / "experiments"), *sys.path])
    for name in MODULES:
        monkeypatch.delitem(sys.modules, name, raising=False)
    import baselines
    yield baselines
    for name in MODULES:
        sys.modules.pop(name, None)


def noise_rows(repo):
    with open(repo / "experiments" / "results" / "noise.csv", newline="") as f:
        return list(csv.DictReader(f))


def test_cases_are_those_of_the_noise_study(repo, bl):
    cases = bl.cases()
    rows = noise_rows(repo)
    assert len(cases) == len(rows) == 840
    for case, r in zip(cases, rows):
        # speed and angle to rounding only: np.exp can differ in the last bit between CPUs (NumPy takes
        # another SIMD path with AVX-512) and between NumPy versions
        assert (case[0], case[3], case[4], case[5]) == (r["texture"], int(r["L"]), float(r["noise"]),
                                                        int(r["seed"]))
        assert case[1:3] == pytest.approx((float(r["speed"]), float(r["angle"])), rel=1e-12, abs=0)
    quick = bl.cases(quick=True)
    assert len(quick) == 96 and sorted({c[4] for c in quick}) == [0, 0.5, 1, 2]


def test_cases_are_generated_as_in_noise_py(bl, monkeypatch):
    # noise.run builds its cases inside run(); catch them where they go to the pool
    import common
    import noise
    caught = []
    monkeypatch.setattr(common, "run_cases", lambda cases, pool, worker=None: caught.append(list(cases)) or [])
    monkeypatch.setattr(common, "write_csv", lambda name, rows: None)
    monkeypatch.setattr(noise, "plot", lambda rows, srows: None)
    for quick in (False, True):
        caught.clear()
        noise.run(None, quick)
        assert caught[0] == bl.cases(quick)


# Easy images, no noise, L = 4: every method must find v, or the comparison would be unfair.
# The velocities have vy <= 0 so that they equal true_velocity of their case (angle in [0, 180)).
EASY = [(5.0, -3.0), (-5.0, -3.0), (2.5, -7.0), (6.0, 0.0), (0.75, -6.4)]


@pytest.mark.parametrize("texture", ["random", "surface"])
def test_every_method_finds_an_easy_motion(bl, monkeypatch, texture):
    import common
    common.init_worker()
    # 128 px images and a 256 px piece of the camera texture instead of 256 and 1024 px, for run time
    monkeypatch.setattr(common, "N", 128)
    monkeypatch.setattr(common, "_texture", common._texture[384:640, 384:640])
    for v in EASY:
        case = (texture, math.hypot(*v), math.degrees(math.atan2(-v[1], v[0])) % 180, 4, 0.0, 3)
        assert bl.true_velocity(case) == pytest.approx(v, abs=1e-12)
        rows = bl.evaluate(case)
        assert [r["method"] for r in rows] == list(bl.METHODS)
        for r in rows:
            assert r["correct"] == 1, (v, r)                # within 2 deg and 5 %, as in the study
        pc = rows[-1]
        assert pc["vector_err"] < 0.02, (v, pc)             # the shift with its sign: v itself, not -v
        bp, frames = bl.simulate(case)
        a, b = bl.single_exposures(case, frames)
        assert bl.phase_correlation_shift(a, b) == pytest.approx(v, abs=0.02)
        assert pc["vector_err_parabolic"] < 0.02, (v, pc)  # fftvel.phase_correlation refines by upsampling too


def test_upsampled_peak_of_a_pure_shift(bl):
    # the lag map of exp(-2 pi i f.tau) is a Dirichlet kernel centred on tau
    m = 64
    k = np.fft.fftfreq(m) * m
    for tau in [(3.3, -7.45), (-0.5, 0.25), (12.0, 5.9)]:
        F = np.exp(-2j * np.pi * np.add.outer(k * tau[1], k * tau[0]) / m)
        got = bl._upsampled_peak(F, round(tau[0]), round(tau[1]))
        assert got == pytest.approx(tau, abs=2e-3)


def test_no_peak_gives_nan_and_not_correct(bl):
    assert all(math.isnan(x) for x in bl._strongest_peak(np.ones((16, 16)), 8))
    row = bl._row(("random", 3.0, 10.0, 4, 0.0, 1), "cepstrum", *bl._angle_speed(float("nan"), float("nan")))
    assert row["correct"] == 0 and math.isnan(row["speed_err"])


def test_evaluate_in_spawned_workers_and_replot(repo, bl, tmp_path, capsys):
    # spawn as on Windows: evaluate and the worker set-up must be top-level and importable
    import common
    cases = [bl.cases()[1], bl.cases()[420]]          # random and camera texture, no noise
    assert [c[0] for c in cases] == ["random", "surface"] and cases[1][4] == 0
    with multiprocessing.get_context("spawn").Pool(2, initializer=common.init_worker) as pool:
        rows = [r for rs in common.run_cases(cases, pool, bl.evaluate) for r in rs]
    assert [r["method"] for r in rows] == list(bl.METHODS) * 2
    noise = {int(r["seed"]): r for r in noise_rows(repo)}
    for r in rows:
        assert r["correct"] == 1, r                    # easy cases: every method is right
        if r["method"] == "phasecorr":
            assert r["vector_err"] < 0.02 and r["correct_parabolic"] == 1
        if r["method"] == "fftvel":                    # the estimate of the noise study, to rounding
            ref = noise[r["seed"]]                     # (the last bits depend on the CPU and NumPy)
            assert r["valid"] == int(ref["valid"])
            assert (r["est_angle"], r["est_speed"]) == \
                pytest.approx((float(ref["est_angle"]), float(ref["est_speed"])), rel=1e-9, abs=0)

    common.set_output(tmp_path)
    common.write_csv("baselines.csv", rows)
    t = bl.table(rows)
    assert {(x["texture"], x["noise"]) for x in t} == {("random", 0), ("random", "all"), ("surface", 0),
                                                       ("surface", "all")}
    assert all(x["correct"] == 100 for x in t)
    assert all(x["valid"] == 100 and x["wrong_valid"] == 0 for x in t if x["method"] == "fftvel")
    capsys.readouterr()
    bl.main(["--replot", "--out", str(tmp_path)])
    out = capsys.readouterr().out
    assert "Baselines: share correct" in out and "(2 images;" in out
    assert "Phase correlation without noise" in out and "100.0 % with fftvel.phase_correlation" in out
    assert (tmp_path / "figures" / "baselines.png").stat().st_size > 10000


@pytest.mark.parametrize("args, code, message", [
    (["--quick", "--out", "{committed}"], 2, "--quick must not write to experiments/results"),
    (["--quick", "--out", "{docs}"], 2, "--quick must not write to docs/figures"),
    (["--jobs", "0"], 2, "--jobs must be at least 1"),
    (["--replot", "--out", "{tmp}"], 1, "baselines.csv not found"),
])
def test_argument_errors(repo, bl, tmp_path, capsys, args, code, message):
    args = [a.format(committed=repo / "experiments" / "results", docs=repo / "docs", tmp=tmp_path) for a in args]
    with pytest.raises(SystemExit) as e:
        bl.main(args)
    if code == 1:
        assert message in str(e.value.code)
    else:
        assert e.value.code == code and message in capsys.readouterr().err
