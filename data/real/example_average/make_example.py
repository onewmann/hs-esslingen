"""Write the synthetic example of average mode in data/real: three stacks of
short-exposure frames and meta.csv.

    python data/real/example_average/make_example.py [OUT_DIR]

Each stack holds 9 frames of 192 x 192 px of data/surface_texture.png (a
real Basler texture), shifted by exact band-limited (Fourier) shifts and
blurred over the exposure of each frame (0.5 ms of a 5 ms frame period,
duty 0.1), with noise of 2 counts added to every frame and stored as 8-bit
counts like Mono8 frames. They are saved as .npz files in the format that
python/examples/live_camera.py --save writes: frames (rows x cols x frames),
t (s), exposure_s, frame_rate and clock.

* moving_regular.npz: 6 px per frame at 30 deg, regular timestamps.
* moving_dropped.npz: the same motion, but the frame after the fourth is
  missing (as after a dropped frame): one step and one frame interval are
  twice as long. The step check of experiments/real_measurements.py must
  flag it.
* static.npz: the stopped surface; it must come out not valid.

The errors that experiments/real_measurements.py reports for
moving_regular.npz come from the estimator and the noise alone, not from a
set speed, a gsd or a frame rate that is off.
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

N = 192                   # side of the frames in px
K = 9                     # frames per stack
PERIOD = 0.005            # frame period in s (200 frames per second)
EXPOSURE = 0.0005         # exposure of one frame in s, duty 0.1
GSD = 0.05                # mm per px
SPEED, ANGLE = 60.0, 30.0     # mm/s and deg: 6 px per frame
OFFSET, GAIN = 70, 600    # counts: mean level and counts per unit of the texture (mean removed: -0.09 to 0.08)
NOISE = 2.0               # standard deviation of the noise in counts
COLUMNS = ["file", "mode", "frames", "first_frame", "frame_period_s", "exposure_s", "duty", "gsd_mm_per_px",
           "set_speed_mm_s", "set_angle_deg", "static", "roi", "note"]

# file, moving, frame positions (in frame periods), frame_period_s in meta.csv, seed, note
CASES = [
    ("moving_regular.npz", True, list(range(9)), PERIOD, 21,
     "synthetic, 6 px per frame, regular timestamps"),
    ("moving_dropped.npz", True, [0, 1, 2, 3, 5, 6, 7, 8, 9], None, 22,
     "synthetic, 6 px per frame, one frame missing (one step twice as long)"),
    ("static.npz", False, list(range(9)), None, 23, "synthetic, stopped surface"),
]


def frames_at(texture, positions, moving):
    """Frames of the texture at the given positions (in frame periods), each blurred over its exposure."""
    v_px = SPEED * PERIOD / GSD if moving else 0.0
    a = np.deg2rad(ANGLE)
    v = (v_px * np.cos(a), -v_px * np.sin(a))           # image coordinates, y down
    n = positions[-1] - positions[0] + 1                 # frames simulated, at steps of one period
    _, frames = simulate_pulsed_image(N, v, 1, exposures=n, duty=EXPOSURE / PERIOD, texture=texture,
                                      return_frames=True)
    return frames[:, :, [p - positions[0] for p in positions]]


def to_uint8(frames, rng):
    """Texture units to 8-bit counts with noise of NOISE counts in every frame."""
    counts = OFFSET + GAIN * frames + NOISE * rng.standard_normal(frames.shape)
    return np.clip(np.round(counts), 0, 255).astype(np.uint8)


def main(out=HERE):
    import matplotlib.image as mpimg

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    texture = mpimg.imread(REPO / "data" / "surface_texture.png").astype(float)
    rows = []
    for name, moving, positions, period, seed, note in CASES:
        rng = np.random.default_rng(seed)
        frames = to_uint8(frames_at(texture, positions, moving), rng)
        t = PERIOD * np.asarray(positions, dtype=float)
        np.savez_compressed(out / name, frames=frames, t=t, exposure_s=EXPOSURE, frame_rate=1 / PERIOD,
                            clock="camera")
        rows.append(dict(file=name, mode="average", frames=K, first_frame="",
                         frame_period_s="" if period is None else f"{period:g}", exposure_s=f"{EXPOSURE:g}",
                         duty="", gsd_mm_per_px=f"{GSD:g}", set_speed_mm_s=f"{SPEED:g}" if moving else "",
                         set_angle_deg=f"{ANGLE:g}" if moving else "", static=0 if moving else 1, roi="",
                         note=note))
    with open(out / "meta.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(CASES)} stacks of {K} frames and meta.csv to {out}")


if __name__ == "__main__":
    main(*sys.argv[1:2])
