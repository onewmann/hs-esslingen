"""experiments/calibrate_threshold.py: the design, one image of each kind, the statistics, the CSV
round trip and the argument checks. Runs a handful of small images, not the study."""

import math
import subprocess
import sys

import numpy as np
import pytest

MODULES = ["calibrate_threshold", "common"]


@pytest.fixture
def ct(repo, monkeypatch, tmp_path):
    """experiments/calibrate_threshold.py imported with its output folder set to tmp_path."""
    monkeypatch.setattr(sys, "path", [str(repo / "experiments"), *sys.path])
    for name in MODULES:
        monkeypatch.delitem(sys.modules, name, raising=False)
    import matplotlib
    matplotlib.use("Agg")
    import calibrate_threshold
    import common
    common.set_output(tmp_path)
    (tmp_path / "figures").mkdir()
    yield calibrate_threshold
    for name in MODULES:
        sys.modules.pop(name, None)


def tiny_cases(ct, recordings):
    """One small image of every kind in every part."""
    cases = [c for c in ct.calibration_cases(quick=True) if c["n"] == 32 and c["L"] == 1][::2]
    cases += [c for c in ct.held_out_cases(recordings, quick=True) if c["n"] == 32 and c["L"] == 1]
    cases += [c for c in ct.power_cases(quick=True) if c["n"] == 32][:2]
    cases += [c for c in ct.uneven_cases(quick=True)][::3]
    cases += [c for c in ct.real_cases(recordings, quick=True) if c["n"] == 256][:1]
    cases += [c for c in ct.blur_cases(quick=True)][-1:]
    cases += [c for c in ct.raw_cases(recordings, quick=True) if c["n"] == 32][:1]
    return cases


def test_design_of_the_full_run(ct, repo):
    rec = repo / "data" / "recordings"
    cal, held = ct.calibration_cases(), ct.held_out_cases(rec)
    real, power, uneven, blur = ct.real_cases(rec), ct.power_cases(), ct.uneven_cases(), ct.blur_cases()
    raw = ct.raw_cases(rec)
    assert len(cal) == 7 * 4 * 3 * 12 == 1008
    assert len(held) == 7 * 4 * 6 * 9 and len(real) == 4 * 10 * (28 + 6)
    assert len(power) == 7 * 4 * 48 and len(uneven) == 2 * 3 * 2 * 40 and len(blur) == 2 * 7 * 4 * 3
    assert len(raw) == 7 * 4 * (12 + 9) == 588
    seeds = [c["seed"] for c in cal + held + power + uneven + blur]
    assert len(seeds) == len(set(seeds))                       # every image has its own seed
    # raw: the surface crops of parts 1 and 2 again, same seeds (so the same positions)
    assert [c["seed"] for c in raw] == [c["seed"] for c in cal + held if c["kind"] == "surface"]
    assert {c["kind"] for c in cal} == set(ct.NULL_KINDS) and {c["L"] for c in cal} == {1, 2, 4, 10}
    for c in power + uneven:                                   # motion inside the measurable range
        K = 2 * c["L"] + 1
        assert 4 / K <= c["speed"] <= c["n"] / 4 and 0 <= c["angle"] < 180
    cells = {(c["kind"], c["frame"], c["n"], c["cell"]) for c in real}
    assert len(cells) == len(real)                             # non-overlapping crops, each once


def test_rows_of_one_image_per_kind(ct, repo):
    rec = repo / "data" / "recordings"
    if not (rec / "vid1.mat").is_file():
        pytest.skip("data/recordings not available")
    cases = tiny_cases(ct, rec)
    rows = [ct.run_case(c) for c in cases]
    assert {r["part"] for r in rows} == {"calibration", "held_out", "power", "uneven", "real", "blur", "raw"}
    for r in rows:
        assert r["k"] == pytest.approx(math.atanh(r["quality"]) * math.sqrt(r["m_eff"] - 5), rel=1e-12)
        assert ct.passes(r, 9.0) == bool(r["valid"])             # same rule as estimate_velocity
        assert r["threshold"] == pytest.approx(max(0.6, math.tanh(9 / math.sqrt(r["m_eff"] - 5))), rel=1e-12)
        assert ct.passes(r, 0.0) >= ct.passes(r, 9.0) >= ct.passes(r, 100.0)
    for r in rows:
        if r["part"] in ("power", "uneven"):
            assert r["correct"] == int(r["angle_err"] <= 2 and r["speed_err"] <= 0.05)
            assert (r["failure"] == "") == bool(r["correct"])
        if r["part"] == "blur":                                # fooled, a blur of d px reads as K flashes d/K apart
            assert r["blur_px"] == 6.0 and r["k_speed"] == pytest.approx(r["est_speed"] * (2 * r["L"] + 1))
        camera = r["part"] == "real" or r["kind"] == "vid1"    # share of saturated pixels of camera crops only
        assert ("saturated" in r) == camera and (not camera or 0.0 <= r["saturated"] <= 1.0)
    assert ct.run_case(cases[0]) == rows[0]                    # fixed seeds: the same row again

    # raw: the same crop as the surface image of the same seed, without the added noise
    rc = next(c for c in cases if c["part"] == "raw")
    sc = dict(rc, part="calibration", kind="surface")
    (raw_img, raw_where), (noisy_img, noisy_where) = ct.null_image(rc), ct.null_image(sc)
    assert raw_where == noisy_where
    assert np.std(noisy_img - raw_img, ddof=1) == pytest.approx(ct.STATIC_NOISE * np.std(raw_img, ddof=1), rel=0.2)

    # analysis, CSV round trip and figure on these rows
    by_part = {name: [r for r in rows if r["part"] == part] for name, part in
               zip(ct.CSV, ("calibration", "held_out", "real", "power", "uneven", "blur", "raw"), strict=True)}
    res = ct.analyse(*by_part.values())
    nroc = 2 * len(ct.K_THRESHOLDS) + 2 * len(ct.GEBA_THRESHOLDS) + len(ct.HALF_CUTS)
    assert len(res["roc"]) == nroc
    cur = next(r for r in res["roc"] if r["rule"] == "k_floor" and r["threshold"] == 9.0)
    assert cur["accepted_null"] == sum(r["valid"] for r in res["null"])
    import common
    for name, sel in by_part.items():
        common.write_csv(name, sel)
    ct.replot()
    assert (common.FIGURES / "threshold.png").is_file()
    back = ct.analyse(*(common.read_csv(n) for n in ct.CSV))
    assert back["cal_stats"]["k_max"] == res["cal_stats"]["k_max"]
    assert [r["accepted_correct"] for r in back["roc"]] == [r["accepted_correct"] for r in res["roc"]]


@pytest.mark.parametrize("n", [256, 512])
def test_independent_points_from_first_principles(ct, n):
    """m counted independently: integer radii 2 <= j <= 0.9*rmax of the padded spectrum, of which n/npad
    are independent. With min_quality 0 the threshold of estimate_velocity shows m even where the
    floor 0.6 hides it in the study (n >= 256)."""
    from fftvel import estimate_velocity
    img = np.random.default_rng(n).standard_normal((n, n))
    est = estimate_velocity(img, L=2, roi=n, min_quality=0.0, diagnostics=True)
    npad = 2 * n                                               # pad_factor 2
    rmax = npad // 2 - 1
    npts = sum(1 for j in range(rmax + 1) if 2 <= j <= 0.9 * rmax)
    m = min(npts, (2 * n * npts) // npad)
    assert est.npad == npad and ct.independent_points(est) == m
    assert est.threshold == pytest.approx(math.tanh(9 / math.sqrt(m - 5)), rel=1e-12)


def test_uneven_flash_times(ct):
    alt = dict(kind="alternating", L=2, param=0.2)
    assert np.allclose(np.diff(ct.flash_positions(alt)), [0.8, 1.2, 0.8, 1.2])
    drop = dict(kind="dropped", L=2, param=3)
    p = ct.flash_positions(drop)
    assert np.allclose(np.diff(p), [1, 1, 2, 1]) and p.size == 5 and p.min() == -p.max()


def test_statistics(ct):
    assert ct.upper_bound(0, 1008) == pytest.approx(1 - 0.05 ** (1 / 1008), rel=1e-9)
    assert ct.upper_bound(3, 3) == 1.0 and math.isnan(ct.upper_bound(0, 0))
    from scipy.optimize import brentq
    from scipy.stats import binom                              # the bound p solves P(X <= k | n, p) = 0.05
    p = brentq(lambda p: binom.cdf(59, 1360, p) - 0.05, 1e-9, 1 - 1e-9, xtol=1e-15)
    assert ct.upper_bound(59, 1360) == pytest.approx(p, rel=1e-9)
    assert ct.auc([2, 3], [1]) == 1.0 and ct.auc([1], [1]) == 0.5 and ct.auc([0], [1, 2]) == 0.0
    assert math.isnan(ct.auc([float("nan")], [1.0]))
    assert ct.null_k(0.0, 30) == 0.0 and ct.null_k(math.tanh(9 / 5), 30) == pytest.approx(9.0)


@pytest.mark.parametrize("args, code, message", [
    (["--quick", "--out", "{committed}"], 2, "--quick must not write to experiments/results"),
    (["--quick", "--out", "experiments/../experiments/results/"], 2, "--quick must not write to experiments/results"),
    (["--quick", "--out", "{docs}"], 2, "--quick must not write to docs/figures"),
    (["--recordings", "{tmp}", "--out", "{tmp}"], 2, "vid1.mat, vid2.mat, vid3.mat, vid4.mat not found"),
    (["--jobs", "0", "--out", "{tmp}"], 2, "--jobs must be at least 1"),
    (["--replot", "--out", "{tmp}"], 1, "threshold_calibration.csv"),
])
def test_arguments(repo, tmp_path, args, code, message):
    args = [a.format(tmp=tmp_path, committed=repo / "experiments" / "results", docs=repo / "docs") for a in args]
    r = subprocess.run([sys.executable, str(repo / "experiments" / "calibrate_threshold.py"), *args],
                       capture_output=True, text=True, timeout=120, cwd=repo)
    assert r.returncode == code and message in r.stderr, r.stderr
