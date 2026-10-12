"""Write the shared parity inputs and the Python reference results.

    python tests/parity/make_cases.py

creates ``cases.mat`` (input images, float32; uint8 for the frame stack) and
``expected.csv`` (results of the Python implementation).
``tests/matlab/test_parity.m`` and ``tests/python/test_parity.py`` check that
MATLAB/Octave and Python both reproduce ``expected.csv`` from ``cases.mat``.

It also writes ``fixed_inputs.mat``: images without motion (simulated static
scenes and white noise, 32 and 64 px), a static repeating texture, and
inputs of the MATLAB tests that would otherwise be drawn at random and fail
for some draws: the 8-bit image of the hot-pixel check and the image of 8
exposures at 0.6 px per step (``tests/matlab/test_model_options.m``) and the
texture of ``tests/matlab/test_phase_correlation.m``. These three are the
images that ``tests/python/test_model_options.py`` and
``tests/python/test_phase_correlation.py`` draw. All are saved once so that
both test suites check the same images whatever their random generator
(MATLAB and Octave draw different numbers from the same seed).

Columns of ``expected.csv`` after ``valid`` give the options of each case
(0 means the default) and further results:
``exposures`` (0: 2L+1), ``duty``, ``roi`` (0: largest square), ``clip``
(clip_outliers), ``stack`` (1: the input is a rows x cols x 1 x frames uint8
stack, averaged with pulsed_from_stack first), ``threshold``,
``quality_perp``, ``isotropic``, ``reject_blur`` (option), ``uint8`` (1: the
image is stored and passed as uint8, not converted to float),
``quality_blur``, ``blur_length``, ``blur_like`` and ``saturated``.
New cases go at the end of CASES: all cases draw from one random generator
in turn, so the inputs of the earlier cases stay the same.
"""

from pathlib import Path

import numpy as np
from scipy.io import savemat

from fftvel import (angle_difference, estimate_velocity, phase_correlation, pulsed_from_stack,
                    simulate_pulsed_image)
from fftvel.core import crop_square
from fftvel.simulate import make_texture

HERE = Path(__file__).resolve().parent

# (n_rows, n_cols, vx, vy, L, pad_factor, noise, options)
# options: exposures, duty, roi, clip, reject_blur (estimator); sim_exposures, sim_duty (simulator);
# hot (two saturated pixels on an 8-bit image), stack (frames in a uint8 stack), tile;
# blur (one continuous exposure over this many px, no flashes), sat (uint8 image with a saturated
# block), glints (two bright spots clipped at a float maximum)
CASES = [
    (128, 128, 3.0, -2.0, 4, 2, 0.0, {}),
    (128, 128, -1.2, 0.7, 4, 2, 0.0, {}),
    (128, 128, 6.0, 6.0, 4, 2, 0.3, {}),
    (128, 128, 0.4, 2.5, 10, 2, 0.0, {}),
    (150, 190, 4.0, 1.0, 4, 2, 0.0, {}),     # non-square: cropped to 150x150
    (101, 101, 2.0, -3.0, 4, 3, 0.0, {}),    # odd n and odd FFT length
    (128, 128, 0.0, 0.0, 4, 2, 0.0, {}),     # no motion: must come out invalid
    (96, 96, 12.0, -1.0, 3, 2, 0.1, {}),
    (96, 96, 1.5, 1.0, 4, 2, 0.0, {"sim_exposures": 8, "exposures": 8}),    # even flash count
    (96, 96, 9.0, 4.0, 4, 2, 0.0, {"sim_duty": 0.3, "duty": 0.3}),          # blur in each exposure
    (120, 150, 3.0, -1.0, 4, 2, 0.0, {"roi": 64}),                          # central 64 px only
    (96, 96, 1.5, -0.5, 4, 2, 0.0, {"hot": True, "clip": True}),            # two hot pixels clipped
    (96, 96, 2.0, -1.5, 5, 2, 0.0, {"stack": 10, "exposures": 10}),         # 10-frame camera stack
    (64, 64, 2.0, 1.0, 4, 4, 0.0, {}),                                      # threshold with pad 4
    (96, 96, 0.0, 0.0, 4, 2, 0.0, {"tile": True}),                          # repeating texture: isotropic
    (128, 128, 0.0, 0.0, 4, 2, 0.05, {"blur": 6.0}),                         # single exposure, blur 6 px
    (128, 128, 0.0, 0.0, 4, 2, 0.05, {"blur": 6.0, "reject_blur": True}),    # the same with reject_blur
    (128, 128, 7.0, 3.0, 4, 2, 0.0, {"reject_blur": True}),                 # fast flashes: not blur-like
    (96, 96, 5.0, 2.0, 4, 2, 0.0, {"sat": True}),                           # uint8 with 255 pixels
    (96, 96, 2.5, -1.0, 4, 2, 0.0, {"glints": True, "clip": True}),         # float maximum, clipped glints
]


def tile_texture(m, period, overlay, seed):
    """Random tile repeated every ``period`` px plus an independent random overlay."""
    rng = np.random.default_rng(seed)
    tile = make_texture(period, 1.0, rng)
    reps = m // period + 1
    return np.tile(tile, (reps, reps))[:m, :m] + overlay * make_texture(m, 1.5, rng)


def eight_bit(bp):
    """Texture with mean 100 and standard deviation 15 DN, rounded to 0..255."""
    return np.clip(np.round(100 + 15 * bp / bp.std()), 0, 255)


def blurred(n, length, angle_deg, rng):
    """One continuous exposure of a random texture moving by ``length`` px: 10 exposures of duty 0.99."""
    a = np.deg2rad(angle_deg)
    v = (length / 10 * np.cos(a), -length / 10 * np.sin(a))
    return simulate_pulsed_image(n, v, 4, exposures=10, duty=0.99, rng=rng)


def with_glints(x, top):
    """``x`` plus two Gaussian spots (sigma 2.5 px) of peak 20 standard deviations, clipped at ``top``."""
    y, c = np.mgrid[0:x.shape[0], 0:x.shape[1]]
    spots = sum(np.exp(-((c - cx) ** 2 + (y - cy) ** 2) / (2 * 2.5 ** 2)) for cy, cx in ((30, 40), (60, 70)))
    return np.minimum(x + 20 * x.std() * spots, top)


def shifted(tex, dx, dy):
    """``tex`` moved by ``(dx, dy)`` px (x along columns, y down), periodic band-limited shift."""
    fy = np.fft.fftfreq(tex.shape[0])
    fx = np.fft.fftfreq(tex.shape[1])
    return np.real(np.fft.ifft2(np.fft.fft2(tex) * np.exp(-2j * np.pi * np.add.outer(fy * dy, fx * dx))))


def check_drawn_inputs(out):
    """The images drawn for the MATLAB tests pass their checks with half the tolerance or less."""
    clean = out["hot96"].astype(float)
    x = clean.copy()
    x[30, 30] = 255
    x[50, 60] = 255
    off = estimate_velocity(x)
    on = estimate_velocity(x, clip_outliers=True)
    print(f"hot96: without clipping {off.speed:.4f} px, clipped {on.speed:.4f} px at {on.angle_deg:.3f} deg")
    assert off.valid and abs(off.speed / np.hypot(20, 30) - 1) < 0.005
    assert on.valid and abs(angle_difference(on.angle_deg, np.rad2deg(np.arctan2(0.5, 1.5)))) < 0.25
    assert abs(on.speed / np.hypot(1.5, 0.5) - 1) < 0.005
    assert estimate_velocity(clean, clip_outliers=True) == estimate_velocity(clean)
    e = estimate_velocity(out["exposures8"].astype(float), exposures=8)
    print(f"exposures8: {e.speed:.4f} px at {e.angle_deg:.3f} deg")
    assert e.valid and abs(angle_difference(e.angle_deg, 25)) < 0.25 and abs(e.speed / 0.6 - 1) < 0.005
    tex = out["texture256"].astype(float)
    a = tex[40:104, 20:116]
    b = shifted(tex, 5.3, -2.6)[40:104, 20:116]
    err = [np.subtract(phase_correlation(a, b)[:2], (5.3, -2.6)),
           np.subtract(phase_correlation(b, a)[:2], (-5.3, 2.6)),
           np.subtract(phase_correlation(a.T, b.T)[:2], (-2.6, 5.3))]
    sign = np.abs(err).max()
    a = tex[::4, ::2]
    half = max(np.abs(np.subtract(phase_correlation(a, shifted(a, *d))[:2], d)).max()
               for d in ((63.6, 0.3), (-63.6, -0.3), (0.3, 31.6), (-0.3, -31.6)))
    row = abs(phase_correlation(tex[:1, 40:168], shifted(tex[:1, :], -5.7, 0)[:, 40:168])[0] + 5.7)
    print(f"texture256: errors {sign:.4f} (sign), {half:.4f} (near half the size), {row:.4f} px (single row)")
    assert sign < 0.01 and half < 0.025 and row < 0.015


def fixed_inputs():
    """Images for the no-motion and isotropy tests, checked here for a clear margin, and images
    that the MATLAB tests would otherwise draw (the same calls as in the Python tests)."""
    out = {}
    for k in (1, 2, 3):
        out[f"static32_{k}"] = simulate_pulsed_image(32, (0, 0), 4, rng=k)
        out[f"noise32_{k}"] = np.random.default_rng(k).standard_normal((32, 32))
        out[f"noise64_{k}"] = np.random.default_rng(1000 + k).standard_normal((64, 64))
    out["tile16"] = crop_square(tile_texture(512, 16, 0.3, 7), 128)
    # tests/python/test_model_options.py: hot_pixel_image(0) before the hot pixels are set, and
    # 8 exposures at 0.6 px per step and 25 deg; tests/python/test_phase_correlation.py: the texture
    a = np.deg2rad(25)
    drawn = {"hot96": eight_bit(simulate_pulsed_image(96, (1.5, -0.5), 4, rng=0)),
             "exposures8": simulate_pulsed_image(256, (0.6 * np.cos(a), -0.6 * np.sin(a)), 4, exposures=8, rng=8),
             "texture256": make_texture(256, 1.5, np.random.default_rng(11))}
    out = {k: v.astype(np.float32) for k, v in out.items()}
    drawn = {k: v.astype(np.float32) for k, v in drawn.items()}
    check_drawn_inputs(drawn)
    for name, img in out.items():
        pads = (2, 6, 8) if name.startswith("noise") else (2,)
        for pad in pads:
            e = estimate_velocity(img.astype(float), pad_factor=pad)
            print(f"{name} pad {pad}: quality {e.quality:.3f} threshold {e.threshold:.3f} "
                  f"isotropic {e.isotropic} valid {e.valid}")
            if name == "tile16":
                assert e.quality > e.threshold + 0.1 and e.in_range and e.isotropic, name
            else:
                assert e.quality < e.threshold - 0.05 and not e.valid, (name, pad)
    return out | drawn


def main():
    rng = np.random.default_rng(20251008)
    mats = {}
    rows = []
    for i, (h, w, vx, vy, L, pad, noise, opt) in enumerate(CASES, start=1):
        n = max(h, w)
        sim = dict(noise=noise, rng=rng, exposures=opt.get("sim_exposures"), duty=opt.get("sim_duty", 0.0))
        exposures = opt.get("exposures")
        if "stack" in opt:
            _, fr = simulate_pulsed_image(n, (vx, vy), L, return_frames=True, **sim)
            fr = np.round(255 * (fr - fr.min()) / np.ptp(fr)).astype(np.uint8)
            x = fr[:, :, None, :opt["stack"]]                # rows x cols x 1 x frames, as from getdata
            bp, K = pulsed_from_stack(x, L=L, exposures=exposures)
            assert K == (exposures or 2 * L + 1)
        else:
            if opt.get("tile"):
                x = tile_texture(512, 16, 0.3, 8)[:h, :w]
            elif "blur" in opt:
                x = blurred(n, opt["blur"], 30.0, rng)
                x = x + noise * x.std() * rng.standard_normal(x.shape)
            else:
                x = simulate_pulsed_image(n, (vx, vy), L, **sim)[:h, :w]
            if opt.get("hot"):
                x = eight_bit(x)
                x[30, 30] = 255                               # static hot pixels 36 px apart
                x[50, 60] = 255
            if opt.get("glints"):
                x = with_glints(x, x.max() + 0.5 * x.std())
            if opt.get("sat"):
                x = eight_bit(x)
                x[10:20, 60:75] = 255                         # a saturated block, 150 px
                x = x.astype(np.uint8)                        # passed as uint8
                bp = x
            else:
                x = x.astype(np.float32)
                bp = x.astype(float)
        mats[f"bp{i:02d}"] = x
        est = estimate_velocity(bp, L=L, pad_factor=pad, exposures=exposures, duty=opt.get("duty", 0.0),
                                roi=opt.get("roi"), clip_outliers=opt.get("clip", False),
                                reject_blur=opt.get("reject_blur", False))
        rows.append([i, L, pad, est.angle_deg, est.speed, est.geba, est.quality, int(est.valid),
                     exposures or 0, opt.get("duty", 0.0), opt.get("roi", 0), int(opt.get("clip", False)),
                     int("stack" in opt), est.threshold, est.quality_perp, int(est.isotropic),
                     int(opt.get("reject_blur", False)), int(x.dtype == np.uint8 and "stack" not in opt),
                     est.quality_blur, est.blur_length, int(est.blur_like), est.saturated])
    mats["params"] = np.array([[r[1], r[2]] for r in rows], dtype=float)
    savemat(HERE / "cases.mat", mats, do_compression=True)
    header = ("case,L,pad,angle_deg,speed,geba,quality,valid,exposures,duty,roi,clip,stack,"
              "threshold,quality_perp,isotropic,reject_blur,uint8,quality_blur,blur_length,blur_like,saturated")
    np.savetxt(HERE / "expected.csv", np.array(rows, dtype=float), delimiter=",",
               header=header, comments="",
               fmt=["%d", "%d", "%d", "%.17g", "%.17g", "%.17g", "%.17g", "%d",
                    "%d", "%.17g", "%d", "%d", "%d", "%.17g", "%.17g", "%d",
                    "%d", "%d", "%.17g", "%.17g", "%d", "%.17g"])
    for r in rows:
        print(r)
    savemat(HERE / "fixed_inputs.mat", fixed_inputs(), do_compression=True)


if __name__ == "__main__":
    main()
