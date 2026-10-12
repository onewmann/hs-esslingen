"""experiments/model_violations.py: the forward model against the simulator, its parts, a tiny run in
spawned workers, the table, figure and command line."""

import multiprocessing
import sys

import numpy as np
import pytest
from scipy import ndimage

from fftvel import estimate_velocity, simulate_pulsed_image

MODULES = ("model_violations", "common")
N = 64                       # small images keep the test fast


@pytest.fixture
def mv(repo, monkeypatch):
    """The module model_violations with the camera texture loaded, removed again afterwards."""
    monkeypatch.setattr(sys, "path", [str(repo / "experiments"), *sys.path])
    for name in MODULES:
        monkeypatch.delitem(sys.modules, name, raising=False)
    import common
    import model_violations
    common.init_worker()
    yield model_violations
    for name in MODULES:
        sys.modules.pop(name, None)


@pytest.mark.parametrize("i", [0, 1, 3])          # camera, random and camera texture
@pytest.mark.parametrize("K", [9, 8])
def test_fourier_forward_model_is_the_simulator_in_the_ideal_case(mv, i, K):
    img = mv.image_params(i)
    assert img.texture == ("surface" if i % 2 == 0 else "random")
    R = mv.scene(img, N)
    assert R.shape == (N + 2 * mv.REACH,) * 2 and R.mean() == pytest.approx(1, abs=1e-12)
    bp = mv.render(R, img.v, mv.flash_times(K), np.ones(K), N)
    ref = simulate_pulsed_image(N, img.v, L=4, exposures=K, texture=R)
    np.testing.assert_allclose(bp - R.mean(), ref, rtol=0, atol=1e-10 * ref.std())


def test_forward_model_equals_shifting_by_hand_for_integer_steps(mv):
    # v = (3, -2) px per step: right and up. Flash k shows the scene moved by t_k * v, i.e. the window at
    # (y, x) sees the scene at (y - t_k * vy, x - t_k * vx); plain array slicing, no FFT
    n = 128
    img = mv.image_params(0)
    R = mv.scene(img, n)
    r0 = (R.shape[0] - n) // 2
    v = (3.0, -2.0)
    hand = np.zeros((n, n))
    for tk in mv.flash_times(9):
        sy, sx = int(tk * v[1]), int(tk * v[0])
        hand += R[r0 - sy:r0 - sy + n, r0 - sx:r0 - sx + n] / 9
    np.testing.assert_allclose(mv.render(R, v, mv.flash_times(9), np.ones(9), n), hand, rtol=0, atol=1e-12)
    e = estimate_velocity(mv.sensor(hand)[0], L=4)
    assert e.valid and e.angle_deg == pytest.approx(np.degrees(np.arctan2(2, 3)), abs=0.2)
    assert e.speed == pytest.approx(np.hypot(3, 2), rel=0.005)


def _sinusoids(m, count=30, seed=1):
    """Scene made of sinusoids that are periodic on the ``m`` grid, as ``[(ky, kx, phase, amplitude)]``."""
    rng = np.random.default_rng(seed)
    return [(int(rng.integers(-50, 51)), int(rng.integers(-50, 51)), rng.uniform(0, 2 * np.pi),
             0.01 * rng.uniform(0.2, 1)) for _ in range(count)]


def _sinusoid_scene(comps, m, y, x, vy=0.0, vx=0.0, length=0.0):
    """Scene at (y, x), or its mean along the path of ``length`` steps of (vy, vx) centred there."""
    out = np.ones_like(x)
    for ky, kx, ph, a in comps:
        out = out + a * np.cos(2 * np.pi * (ky * y + kx * x) / m + ph) * np.sinc(length * (ky * vy + kx * vx) / m)
    return out


def test_forward_model_on_a_scene_known_in_closed_form(mv):
    # sub-pixel steps, timing and energy jitter, flash duration and continuous light, against the values of a
    # scene of sinusoids, which are known exactly at any shift and for any box average along v
    n = 48
    m = n + 2 * mv.REACH
    comps = _sinusoids(m)
    R = _sinusoid_scene(comps, m, *np.mgrid[0:m, 0:m].astype(float))
    r0 = (m - n) // 2
    y, x = np.mgrid[r0:r0 + n, r0:r0 + n].astype(float)
    v = (2.37, -5.81)
    t = mv.flash_times(9, 0.1, np.random.default_rng(3))
    w = 1 + 0.2 * np.random.default_rng(4).standard_normal(9)
    exact = sum(wk * _sinusoid_scene(comps, m, y - tk * v[1], x - tk * v[0]) for tk, wk in zip(t, w)) / 9
    np.testing.assert_allclose(mv.render(R, v, t, w, n), exact, rtol=0, atol=1e-12)
    # each flash lasting half a step: sub-steps against the exact box average
    t = mv.flash_times(9)
    exact = np.mean([_sinusoid_scene(comps, m, y - tk * v[1], x - tk * v[0], v[1], v[0], 0.5) for tk in t], axis=0)
    rel = mv.render(R, v, t, np.ones(9), n, duty=0.5)
    assert np.abs(rel - exact).max() < 1e-3 * exact.std()
    # continuous light over an exposure of 10 steps: the box average over 10 steps
    exact = _sinusoid_scene(comps, m, y, x, v[1], v[0], 10)
    rel = mv.render(R, v, t, np.ones(9), n, continuous=1.0, exposure_steps=10)
    assert np.abs(rel - exact).max() < 1e-3 * exact.std()


def test_seed_streams_are_separate(mv):
    # numpy pads seeds with zeros ([a, b] and [a, b, 0] are the same stream); every stream must differ
    for i in (0, 1, 7):
        keys = [(), (mv.SCENE_KEY,)] + [(f.code,) for f in mv.FACTORS]
        keys += [(f.code, li + 1) for f in mv.FACTORS for li in range(len(f.levels))]
        first = {tuple(mv.stream(i, *k).integers(0, 2 ** 62, 2)) for k in keys}
        assert len(first) == len(keys), i
    # the scene does not depend on the velocity draw: crop offsets against log speed over many images
    speeds, offsets = [], []
    for i in range(0, 400, 2):
        speeds.append(np.log(mv.image_params(i).speed))
        offsets.append(mv.stream(i, mv.SCENE_KEY).integers(0, 321, 2))
    offsets = np.array(offsets)
    assert max(abs(np.corrcoef(speeds, offsets[:, j])[0, 1]) for j in (0, 1)) < 0.3


def test_flash_duration_by_substeps_is_close_to_the_closed_form(mv):
    img = mv.image_params(2)                       # 9.4 px per step: 19 sub-steps per flash at duty 0.5
    R = mv.scene(img, N)
    bp = mv.render(R, img.v, mv.flash_times(9), np.ones(9), N, duty=0.5)
    ref = simulate_pulsed_image(N, img.v, L=4, duty=0.5, texture=R)
    assert np.std(bp - R.mean() - ref) < 1e-3 * np.std(ref)
    sharp = mv.render(R, img.v, mv.flash_times(9), np.ones(9), N)
    assert np.std(sharp - bp) > 0.1 * np.std(ref)      # the blur itself is large


def test_continuous_light_over_k_steps_is_the_flash_pattern_smeared(mv):
    img = mv.image_params(1)
    R = mv.scene(img, N)
    full = mv.render(R, img.v, mv.flash_times(9), np.ones(9), N, continuous=1.0)
    # continuous light over 9 steps is the same as 9 flashes, each lasting a whole step
    ref = mv.render(R, img.v, mv.flash_times(9), np.ones(9), N, duty=1.0)
    np.testing.assert_allclose(full, ref, rtol=0, atol=1e-12)
    longer = mv.render(R, img.v, mv.flash_times(9), np.ones(9), N, continuous=1.0, exposure_steps=10)
    assert np.std(longer - full) > 0.01 * np.std(full)


def test_spline_render_equals_shifting_the_whole_scene(mv):
    img = mv.image_params(3)
    R = mv.scene(img, N)
    t = mv.flash_times(9)
    r0 = (R.shape[0] - N) // 2
    full = np.mean([ndimage.shift(R, (tk * img.v[1], tk * img.v[0]), order=3, mode="mirror")[r0:r0 + N, r0:r0 + N]
                    for tk in t], axis=0)
    np.testing.assert_allclose(mv.render_spline(R, img.v, t, N), full, rtol=0, atol=1e-12)


def test_flash_times_and_substeps(mv):
    np.testing.assert_array_equal(mv.flash_times(9), np.arange(-4, 5))
    np.testing.assert_array_equal(mv.flash_times(8), np.arange(8) - 3.5)
    t = mv.flash_times(9, 0.1, np.random.default_rng(1))
    assert t.mean() == pytest.approx(0, abs=1e-12) and np.all(np.diff(t) > 0) and np.std(np.diff(t)) > 0.02
    # timing jitter: each flash off its own nominal time, independently (not a random walk of the steps)
    K = 4001
    d = mv.flash_times(K, 0.1, np.random.default_rng(2)) - mv.flash_times(K)
    assert np.std(d) == pytest.approx(0.1, rel=0.05)
    assert abs(np.corrcoef(d[:-1], d[1:])[0, 1]) < 0.1
    tau, count = mv.substeps(2.0)                  # 2 px: 8 sub-steps 0.25 px apart, centred
    assert count == 8 and tau.mean() == pytest.approx(0, abs=1e-15) and np.allclose(np.diff(tau), 1 / 8)
    assert mv.substeps(0.1)[1] == 1


def test_periodic_tile_is_seamless_and_periodic(mv):
    h, w = 12, 16
    y, x = np.mgrid[0:h, 0:w]
    u = 10 * x / w + 0.1 * np.random.default_rng(2).standard_normal((h, w))     # ramp: big jump at the border
    p = mv.periodic_component(u)
    assert p.mean() == pytest.approx(u.mean())
    assert np.abs(p[:, 0] - p[:, -1]).mean() < 0.1 * np.abs(u[:, 0] - u[:, -1]).mean()
    R = 1 + 0.07 * np.random.default_rng(3).standard_normal((200, 200))
    T = mv.tiled(R, 16)
    assert T.shape == R.shape
    np.testing.assert_array_equal(T[:, 16:32], T[:, :16])
    np.testing.assert_array_equal(T[12:24, :], T[:12, :])


def test_disc_otf(mv):
    otf = mv.disc_otf(3, 128)
    assert otf[0, 0] == pytest.approx(1)
    assert np.abs(otf.imag).max() < 1e-12            # symmetric point spread function
    assert otf.real.min() < -0.05                    # a disc, not a Gaussian: negative lobes


def test_sensor(mv):
    rng = np.random.default_rng(4)
    rel = 1 + 0.1 * rng.standard_normal((200, 200))
    out, nsr = mv.sensor(rel)
    np.testing.assert_array_equal(out, mv.EXPOSURE * rel)
    assert nsr == 0
    dn, nsr = mv.sensor(rel, bits=8)
    assert np.all(dn == np.rint(dn)) and nsr > 0
    assert mv.sensor(10 * rel, bits=8)[0].max() == 255                             # clipped at full scale
    noisy, nsr = mv.sensor(np.ones((300, 300)) + 1e-3 * rng.standard_normal((300, 300)), full_well=10000,
                           rng=np.random.default_rng(5))
    expected = np.sqrt(mv.EXPOSURE * 10000 + mv.READ_NOISE ** 2) / 10000
    assert np.std(noisy) == pytest.approx(expected, rel=0.03)


def test_constant_offset_without_noise_changes_nothing(mv):
    img = mv.image_params(1)
    rel = mv.render(mv.scene(img, N), img.v, mv.flash_times(9), np.ones(9), N)
    a = estimate_velocity(mv.sensor(rel)[0], L=4)
    b = estimate_velocity(mv.sensor(0.1 * rel + 0.9)[0], L=4)
    assert b.speed == pytest.approx(a.speed, rel=1e-9) and b.angle_deg == pytest.approx(a.angle_deg, abs=1e-7)


def test_make_image_is_reproducible_and_paired_over_levels(mv):
    for f in mv.FACTORS:
        last = len(f.levels) - 1
        out, img, nsr = mv.make_image(f.name, last, 5, N)
        again = mv.make_image(f.name, last, 5, N)[0]
        assert out.shape == (N, N) and np.all(np.isfinite(out)) and img.index == 5, f.name
        np.testing.assert_array_equal(out, again)
        assert (nsr > 0) == (f.name in ("shot_noise", "quantisation", "offset", "defocus_noisy")), f.name
    # the same glint positions at every level: the difference of two levels is the glint map
    lo, hi = mv.make_image("glints", 1, 6, N)[0], mv.make_image("glints", 2, 6, N)[0]
    ideal = mv.make_image("glints", 0, 6, N)[0]
    g1, g2 = (lo - ideal) / mv.EXPOSURE, (hi - ideal) / mv.EXPOSURE
    np.testing.assert_allclose(g2 * mv.BY_NAME["glints"].levels[1], g1 * mv.BY_NAME["glints"].levels[2], atol=1e-9)
    assert g1.max() == pytest.approx(mv.BY_NAME["glints"].levels[1], rel=0.2)


def test_variants(mv):
    assert mv.variants("defocus", 3) == {"default": {"L": 4}}
    assert mv.variants("pulse_width", 0.5) == {"default": {"L": 4}, "matched": {"L": 4, "duty": 0.5}}
    assert mv.variants("glints", 10)["matched"] == {"L": 4, "clip_outliers": True}
    assert mv.variants("true_K", 8) == {"default": {"L": 4}, "matched": {"L": 4, "exposures": 8}}
    assert mv.variants("assumed_K", 10) == {"default": {"L": 4, "exposures": 10}}
    # the option reaches the estimator: duty 0.5 on a fast image (9.4 px per step) changes the estimate
    default, matched = mv.evaluate(("pulse_width", 3, 2, 128))
    assert (default["variant"], matched["variant"]) == ("default", "matched")
    assert matched["est_speed"] != default["est_speed"] and matched["correct"] == 1


def test_cases(mv):
    checks, sweeps = mv.cases()
    assert len(checks) == len(mv.CHECKS) * mv.PER_LEVEL >= 4 * 30
    assert len(sweeps) == sum(len(f.levels) for f in mv.FACTORS) * mv.PER_LEVEL
    assert all(4 <= len(f.levels) <= 6 for f in mv.FACTORS)
    quick = mv.cases(quick=True)
    assert len(quick[1]) == sum(len(f.levels) for f in mv.FACTORS) * mv.PER_LEVEL_QUICK


def test_tiny_run_in_spawned_workers_table_and_replot(mv, tmp_path, capsys):
    # spawn as on Windows: the workers and their set-up must be top-level and importable
    import common
    n = 128
    checks = [("fourier", 0, n), ("K_8", 1, n)]
    sweeps = [("glints", 0, 0, n), ("glints", 5, 0, n), ("true_K", 1, 1, n), ("assumed_K", 2, 1, n)]
    with multiprocessing.get_context("spawn").Pool(2, initializer=common.init_worker) as pool:
        rows = [r for rs in common.run_cases(checks, pool, mv.check) for r in rs]
        rows += [r for rs in common.run_cases(sweeps, pool, mv.evaluate) for r in rs]
    assert [r["factor"] for r in rows] == ["check"] * 2 + ["glints"] * 4 + ["true_K"] * 2 + ["assumed_K"]
    assert [r["variant"] for r in rows[2:]] == ["default", "matched"] * 3 + ["default"]
    for r in rows[:2]:
        assert abs(r["dev_speed"]) < 1e-9 and abs(r["dev_angle"]) < 1e-9 and r["valid"] == r["sim_valid"]
    ideal = rows[2]
    assert ideal["level"] == 0 and ideal["valid"] == 1 and ideal["correct"] == 1, ideal
    assert rows[3] == {**ideal, "variant": "matched"}                      # no glints: clipping changes nothing
    in_process = mv.evaluate(("glints", 0, 0, n))
    assert in_process[0] == ideal                                          # the same in this process

    common.set_output(tmp_path)
    common.write_csv("violations.csv", rows)
    t = mv.table(common.read_csv("violations.csv"))
    assert [(x["factor"], x["level"], x["variant"]) for x in t] == [
        ("glints", 0, "default"), ("glints", 0, "matched"), ("glints", 100, "default"), ("glints", 100, "matched"),
        ("assumed_K", 9, "default"), ("true_K", 8, "default"), ("true_K", 8, "matched")]      # order of FACTORS
    assert all(x["correct_valid"] + x["wrong_valid"] + x["not_valid"] == pytest.approx(100) for x in t)
    assert mv.check_table(rows)[0]["agree"]
    capsys.readouterr()
    mv.main(["--replot", "--out", str(tmp_path)])
    out = capsys.readouterr().out
    assert "Check: forward model against simulate_pulsed_image" in out and "Model violations" in out
    assert "Wrong but valid results" in out and "(4 images in 7 rows of the table)" in out
    assert (tmp_path / "figures" / "violations.png").stat().st_size > 10000


@pytest.mark.parametrize("args, code, message", [
    (["--quick", "--out", "{committed}"], 2, "--quick must not write to experiments/results"),
    (["--quick", "--out", "{docs}"], 2, "--quick must not write to docs/figures"),
    (["--jobs", "0"], 2, "--jobs must be at least 1"),
    (["--replot", "--out", "{tmp}"], 1, "violations.csv not found"),
])
def test_argument_errors(repo, mv, tmp_path, capsys, args, code, message):
    args = [a.format(committed=repo / "experiments" / "results", docs=repo / "docs", tmp=tmp_path) for a in args]
    with pytest.raises(SystemExit) as e:
        mv.main(args)
    if code == 1:
        assert message in str(e.value.code)
    else:
        assert e.value.code == code and message in capsys.readouterr().err
