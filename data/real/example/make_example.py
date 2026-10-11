"""Write the synthetic example of data/real: three multi-flash images, two
single-flash reference frames and meta.csv.

    python data/real/example/make_example.py [OUT_DIR]

The images are fftvel.simulate_pulsed_image of data/surface_texture.png (a
real Basler texture) with exactly known motion, 256 x 256 px, stored as
16-bit PNG like the frames of a 12-bit camera. They let
experiments/real_measurements.py be tried before real recordings exist; the
errors it reports for them come from the estimator alone, not from a set
speed, a gsd or a flash period that is off.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "python" / "src"))

from fftvel import simulate_pulsed_image  # noqa: E402

N = 256
GAIN = 3 * 65536          # counts per unit of the texture (values 0.21 to 0.41)
COLUMNS = ["file", "exposures", "flash_period_s", "duty", "gsd_mm_per_px", "set_speed_mm_s",
           "set_angle_deg", "static", "roi", "exposure_s", "ref_frame_a", "ref_frame_b", "ref_dt_s",
           "note"]

# file, K, flash period (s), duty, gsd (mm/px), set speed (mm/s), set angle (deg), static, noise, seed
CASES = [
    ("flash_a.png", 9, 1e-3, 0.02, 0.05, 200.0, 30.0, 0, 0.1, 11),
    ("flash_b.png", 8, 5e-4, 0.2, 0.05, 560.0, 120.0, 0, 0.3, 12),
    ("static.png", 9, 1e-3, 0.02, 0.05, None, None, 1, 0.1, 13),
]


def to_uint16(img):
    """Texture units to 16-bit counts around mid-grey."""
    return np.clip(np.round(32768 + GAIN * img), 0, 65535).astype(np.uint16)


def write_png(path, img):
    from PIL import Image
    Image.fromarray(to_uint16(img)).save(path)


def fmt(x):
    return "" if x is None else f"{x:g}"


def main(out=HERE):
    import matplotlib.image as mpimg

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    texture = mpimg.imread(REPO / "data" / "surface_texture.png").astype(float)
    rows = []
    for name, K, T, duty, gsd, speed, angle, static, noise, seed in CASES:
        v_px = 0.0 if static else speed * T / gsd           # motion per flash step in px
        a = np.deg2rad(angle or 0.0)
        v = (v_px * np.cos(a), -v_px * np.sin(a))           # image coordinates, y down
        rng = np.random.default_rng(seed)
        bp, frames = simulate_pulsed_image(N, v, 1, texture=texture, noise=noise, rng=rng,
                                           return_frames=True, exposures=K, duty=duty)
        write_png(out / name, bp)
        row = dict(file=name, exposures=K, flash_period_s=fmt(T), duty=fmt(duty), gsd_mm_per_px=fmt(gsd),
                   set_speed_mm_s=fmt(speed), set_angle_deg=fmt(angle), static=static, roi="",
                   exposure_s=fmt((K - 1 + duty) * T * 1.05), ref_frame_a="", ref_frame_b="", ref_dt_s="",
                   note=("synthetic, stopped surface" if static else
                         f"synthetic, {v_px:g} px per flash, noise {100 * noise:g} % of the signal"))
        if name == "flash_a.png":
            # two single-flash frames, the first and the last flash: (K-1) flash periods apart
            for key, q in (("ref_frame_a", 0), ("ref_frame_b", K - 1)):
                single = frames[:, :, q]
                single = single + noise * single.std(ddof=1) * rng.standard_normal(single.shape)
                ref = f"ref_{key[-1]}.png"
                write_png(out / ref, single)
                row[key] = ref
            row["ref_dt_s"] = fmt((K - 1) * T)
        rows.append(row)
    with open(out / "meta.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(CASES)} multi-flash images, 2 reference frames and meta.csv to {out}")


if __name__ == "__main__":
    main(*sys.argv[1:2])
