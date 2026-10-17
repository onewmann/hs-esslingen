"""Turn the phone video of this session into the frame stack and the meta.csv that
experiments/real_measurements.py reads.

    python data/measurements/2026-10-16_phone_120fps/extract_frames.py
    python experiments/real_measurements.py data/measurements/2026-10-16_phone_120fps/meta.csv --figure data/measurements/2026-10-16_phone_120fps/results.png

The video is read with OpenCV (pip install -e "./python[video]"), which gives
the timestamp of every frame as the phone stored it. The phone stored the
frames lying on their side (1080 x 1920 px, rows x columns, with a display
rotation of -90 deg); the script switches OpenCV's own rotation off, because
OpenCV versions differ in it, and turns each frame upright itself (1920 x
1080 px), a lossless rotation. Each frame is converted to grey and cut to
a 512 x 512 px square below the bright spot of the phone's lamp, which stays at
the same place in the image while the surface moves (rows 1300 to 1811,
columns 284 to 795). The stack goes to frames.npz in the format that
python/examples/live_camera.py --save writes: frames (rows x columns x frames,
8-bit), t (s from the first frame) and clock ('video': the times stored in the
file, without the transfer jitter of host times), plus opencv, the OpenCV
version. OpenCV versions decode the video slightly differently (grey values up
to 2 levels apart); the committed results.csv was made with OpenCV 5.0.0, and
with 4.12 the largest deviations are 0.88 % and 0.58 deg instead of 0.81 % and
0.52 deg. frames.npz is about 75 MB and not kept in git.

meta.csv lists consecutive, non-overlapping windows of 5 frames (mode
average). The phone moved by hand at an unknown speed, so gsd_mm_per_px,
set_speed_mm_s and set_angle_deg are empty: the reference is the step check of
real_measurements.py, the mean step between consecutive frames by phase
correlation.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROW, COL, SIZE = 1300, 284, 512          # the analysed square in the upright frame
STORED = (1080, 1920)                    # rows x columns of a frame as stored in the video
WINDOW = 5                               # frames per window in meta.csv
NOTE = "phone video of the desk top, hand-held, 120 fps"


def read_video(path: Path, row=ROW, col=COL, size=SIZE, max_frames=None):
    """Grey crops (frames x size x size, uint8) and timestamps (s from the first frame) of a video."""
    try:
        import cv2
    except ImportError:
        raise SystemExit('reading the video needs OpenCV: python -m pip install -e "./python[video]"') from None
    cap = cv2.VideoCapture(str(path))
    if not cap.isOpened():
        raise SystemExit(f"cannot open {path}")
    if hasattr(cv2, "CAP_PROP_ORIENTATION_AUTO"):      # versions differ; the frames are turned below
        cap.set(cv2.CAP_PROP_ORIENTATION_AUTO, 0)
    frames, times = [], []
    try:
        while max_frames is None or len(frames) < max_frames:
            ok, img = cap.read()
            if not ok:
                break
            times.append(cap.get(cv2.CAP_PROP_POS_MSEC) / 1000)
            grey = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
            if grey.shape != STORED:
                raise SystemExit(f"frame of {grey.shape[0]} x {grey.shape[1]} px, expected {STORED[0]} x {STORED[1]} "
                                 f"px as stored in the video (this OpenCV may have turned it already)")
            grey = np.rot90(grey, -1)                   # the display rotation of -90 deg: upright, 1920 x 1080 px
            frames.append(grey[row:row + size, col:col + size].copy())
    finally:
        cap.release()
    if not frames:
        raise SystemExit(f"no frame read from {path}")
    t = np.asarray(times, dtype=float)
    return np.stack(frames), t - t[0]


def write_meta(path: Path, n_frames: int, stack: str = "frames.npz", window: int = WINDOW) -> int:
    """meta.csv with consecutive windows of ``window`` frames; returns the number of rows."""
    starts = range(1, n_frames - window + 2, window)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["file", "mode", "frames", "first_frame", "static", "roi", "note"])
        for s in starts:
            w.writerow([stack, "average", window, s, 0, SIZE, NOTE])
    return len(starts)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--video", type=Path, default=HERE / "desk_120fps.mp4")
    p.add_argument("--out", type=Path, default=HERE, help="folder for frames.npz and meta.csv (default: here)")
    a = p.parse_args(argv)
    frames, t = read_video(a.video)
    a.out.mkdir(parents=True, exist_ok=True)
    import cv2
    np.savez_compressed(a.out / "frames.npz", frames=np.moveaxis(frames, 0, 2), t=t, clock="video",
                        opencv=cv2.__version__)
    rows = write_meta(a.out / "meta.csv", len(frames))
    dt = np.diff(t)
    print(f"{len(frames)} frames, median interval {1e3 * np.median(dt):.2f} ms, {int((dt > 1.5 * np.median(dt)).sum())} "
          f"longer gaps; {a.out / 'frames.npz'} and meta.csv with {rows} windows of {WINDOW} frames")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
