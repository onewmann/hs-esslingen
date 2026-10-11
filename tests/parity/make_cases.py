"""Write the shared parity inputs and the Python reference results.

    python tests/parity/make_cases.py

creates ``cases.mat`` (input images, float32; uint8 for the frame stack) and
``expected.csv`` (results of the Python implementation).
``tests/matlab/test_parity.m`` and ``tests/python/test_parity.py`` check that
MATLAB/Octave and Python both reproduce ``expected.csv`` from ``cases.mat``.

It also writes ``fixed_inputs.mat``: images without motion (simulated static
scenes and white noise, 32 and 64 px) and a static repeating texture, saved
once so that both test suites check the same images whatever their random
generator (MATLAB and Octave draw different numbers from the same seed).

Columns of ``expected.csv`` after ``valid`` give the options of each case
(0 means the default) and further results:
``exposures`` (0: 2L+1), ``duty``, ``roi`` (0: largest square), ``clip``
(clip_outliers), ``stack`` (1: the input is a rows x cols x 1 x frames uint8
stack, averaged with pulsed_from_stack first), ``threshold``,
``quality_perp`` and ``isotropic``.
"""

from pathlib import Path

import numpy as np
from scipy.io import savemat

from fftvel import estimate_velocity, pulsed_from_stack, simulate_pulsed_image
from fftvel.core import crop_square
from fftvel.simulate import make_texture

HERE = Path(__file__).resolve().parent

# (n_rows, n_cols, vx, vy, L, pad_factor, noise, options)
# options: exposures, duty, roi, clip (estimator); sim_exposures, sim_duty (simulator);
# hot (two saturated pixels on an 8-bit image), stack (frames in a uint8 stack), tile
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


def fixed_inputs():
    """Images for the no-motion and isotropy tests, checked here for a clear margin."""
    out = {}
    for k in (1, 2, 3):
        out[f"static32_{k}"] = simulate_pulsed_image(32, (0, 0), 4, rng=k)
        out[f"noise32_{k}"] = np.random.default_rng(k).standard_normal((32, 32))
        out[f"noise64_{k}"] = np.random.default_rng(1000 + k).standard_normal((64, 64))
    out["tile16"] = crop_square(tile_texture(512, 16, 0.3, 7), 128)
    out = {k: v.astype(np.float32) for k, v in out.items()}
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
    return out


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
            else:
                x = simulate_pulsed_image(n, (vx, vy), L, **sim)[:h, :w]
            if opt.get("hot"):
                x = eight_bit(x)
                x[30, 30] = 255                               # static hot pixels 36 px apart
                x[50, 60] = 255
            x = x.astype(np.float32)
            bp = x.astype(float)
        mats[f"bp{i:02d}"] = x
        est = estimate_velocity(bp, L=L, pad_factor=pad, exposures=exposures, duty=opt.get("duty", 0.0),
                                roi=opt.get("roi"), clip_outliers=opt.get("clip", False))
        rows.append([i, L, pad, est.angle_deg, est.speed, est.geba, est.quality, int(est.valid),
                     exposures or 0, opt.get("duty", 0.0), opt.get("roi", 0), int(opt.get("clip", False)),
                     int("stack" in opt), est.threshold, est.quality_perp, int(est.isotropic)])
    mats["params"] = np.array([[r[1], r[2]] for r in rows], dtype=float)
    savemat(HERE / "cases.mat", mats, do_compression=True)
    header = ("case,L,pad,angle_deg,speed,geba,quality,valid,exposures,duty,roi,clip,stack,"
              "threshold,quality_perp,isotropic")
    np.savetxt(HERE / "expected.csv", np.array(rows, dtype=float), delimiter=",",
               header=header, comments="",
               fmt=["%d", "%d", "%d", "%.17g", "%.17g", "%.17g", "%.17g", "%d",
                    "%d", "%.17g", "%d", "%d", "%d", "%.17g", "%.17g", "%d"])
    for r in rows:
        print(r)
    savemat(HERE / "fixed_inputs.mat", fixed_inputs(), do_compression=True)


if __name__ == "__main__":
    main()
