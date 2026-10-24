"""Turn the two graph-paper videos of this session into frame stacks and the meta.csv that
experiments/real_measurements.py reads.

    python data/measurements/2026-10-17_phone_grid/extract_frames.py
    python experiments/real_measurements.py data/measurements/2026-10-17_phone_grid/meta.csv --figure data/measurements/2026-10-17_phone_grid/results.png

The frames are read with the reader of the session 2026-10-16_phone_120fps
(OpenCV, turned upright, grey) and cut to the central 512 x 512 px square
(rows 704 to 1215, columns 284 to 795): frames_moving.npz and frames_still.npz,
not kept in git. The 1024 px field of every frame of a video is held in memory
for the grid measurement, about 1.7 GB at the peak for these videos; a longer
video needs proportionally more.

meta.csv lists consecutive windows of 5 frames of both videos. The still video
is listed with static 1: the phone was held still over the paper, and its
windows must come out not valid. For the moving video, the ground sample
distance of each window comes from the paper itself: the squares are 5 mm, and
their size in pixels is the period of the grid lines, measured by
autocorrelation of the column and row means of every frame of the window, on
the 1024 x 1024 px field with the same centre as the analysed square (four
lines of the 512 px square alone are too few for a steady value). The phone was held by hand, so its height and with it the scale
changed during the recording; a window whose period differs by more than 3 %
from the median of the last five accepted windows (a misreading: blur, the
paper tilting, its edge entering the field) gets no gsd_mm_per_px. There is no set speed: the reference is the step check.

grid_lines.csv gives, for every window of both videos, the period of the grid
in px and the direction of the grid lines in the image (deg, counter-clockwise
from +x with y up, as angle_deg in results.csv), measured from the shift of the
row means between the left and right half of the field and of the column means
between its upper and lower half, and the ground sample distance where meta.csv
has one. The hand
turned the phone, so the lines do not stay at 0 deg; directions.py compares
the estimates with them.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SQUARE_MM = 5.0                          # side of a square of the graph paper
ROW, COL, SIZE = 704, 284, 512           # the analysed square in the upright frame
SCALE_FIELD = 1024                       # field around it, with the same centre, for the grid period
WINDOW = 5
VIDEOS = (("grid_moving_120fps.mp4", "frames_moving.npz", 0, "graph paper (5 mm squares), phone moved by hand, 120 fps"),
          ("grid_still_120fps.mp4", "frames_still.npz", 1, "graph paper (5 mm squares), phone held still, 120 fps"))


def _reader():
    path = HERE.parent / "2026-10-16_phone_120fps" / "extract_frames.py"
    spec = importlib.util.spec_from_file_location("phone_120fps_extract", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.read_video


def grid_period(img, lo=70, hi=300):
    """Period of the grid lines in px, from the column and the row means.

    Takes the first local maximum of the (biased) autocorrelation between
    ``lo`` and ``hi`` px that reaches 70 % of the largest one there, so that a
    multiple of the period is not taken.
    """
    from scipy.ndimage import gaussian_filter1d

    out = []
    for prof in (img.mean(axis=0), img.mean(axis=1)):
        p = prof - gaussian_filter1d(prof, 30)          # remove the shading, keep the thin lines
        p = p - p.mean()
        n = p.size
        a = np.correlate(p, p, "full")[n - 1:] / n     # biased: falls with the lag, so the period beats its multiples
        top = min(hi, n // 2)
        peaks = [k for k in range(lo, top) if a[k] >= a[k - 1] and a[k] >= a[k + 1]]
        if not peaks:
            raise ValueError("no grid period found")
        best = max(a[k] for k in peaks)
        k = next(k for k in peaks if a[k] >= 0.7 * best)
        y0, y1, y2 = a[k - 1], a[k], a[k + 1]
        out.append(k + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2))
    return float(np.mean(out))


def _shift(a, b, max_shift=60):
    """Shift in px with b(x) close to a(x - shift), by cross-correlation and a parabola through its peak."""
    from scipy.ndimage import gaussian_filter1d

    a = a - gaussian_filter1d(a, 30)
    b = b - gaussian_filter1d(b, 30)
    a, b = a - a.mean(), b - b.mean()
    n = a.size
    c = np.array([np.dot(a[max(0, -k):n - max(0, k)], b[max(0, k):n - max(0, -k)])
                  for k in range(-max_shift, max_shift + 1)])
    k = int(np.argmax(c))
    if not 0 < k < len(c) - 1:
        return float("nan")                         # the peak at the edge of the search range: shift too large
    return k - max_shift + 0.5 * (c[k - 1] - c[k + 1]) / (c[k - 1] - 2 * c[k] + c[k + 1])


def line_angle(img):
    """Direction of the grid lines in deg, counter-clockwise from +x with y up.

    The horizontal lines lie lower or higher in the right half than in the
    left half, the vertical lines further right or left in the lower half than
    in the upper half; both give the rotation, and their mean is returned.
    Valid while the shift between the halves, half the width times tan(angle),
    stays below half the period and below the 60 px that _shift searches
    (about 6.7 deg on the 1024 px field); beyond that the result is nan or wrong.
    """
    h, w = img.shape
    hh, hw = h // 2, w // 2
    s_rows = _shift(img[:, :hw].mean(axis=1), img[:, hw:].mean(axis=1))
    s_cols = _shift(img[:hh, :].mean(axis=0), img[hh:, :].mean(axis=0))
    return float((np.degrees(np.arctan2(-s_rows, hw)) + np.degrees(np.arctan2(s_cols, hh))) / 2)


def window_gsd(field, window=WINDOW, tol=0.03):
    """gsd in mm/px per window of ``window`` frames of ``field`` (frames x rows x columns), and the periods.

    The scale changes slowly as the hand moves up and down. A window whose
    period differs by more than ``tol`` from the median of the last five
    accepted windows gets None: its period is a misreading (blur, the paper
    tilting, its edge entering the field), not a change of height.
    """
    periods = np.array([np.mean([grid_period(f.astype(float)) for f in field[s:s + window]])
                        for s in range(0, len(field) - window + 1, window)])
    out, accepted = [], []
    for p in periods:
        ok = not accepted or abs(p / np.median(accepted[-5:]) - 1) <= tol
        out.append(SQUARE_MM / p if ok else None)
        if ok:
            accepted.append(p)
    return out, periods


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--out", type=Path, default=HERE, help="folder for the frame stacks and meta.csv (default: here)")
    a = p.parse_args(argv)
    read_video = _reader()
    a.out.mkdir(parents=True, exist_ok=True)
    rows, lines = [], []
    import cv2
    for video, stack, static, note in VIDEOS:
        m = (SCALE_FIELD - SIZE) // 2
        field, t = read_video(HERE / video, ROW - m, COL - m, SCALE_FIELD)
        frames = field[:, m:m + SIZE, m:m + SIZE]
        np.savez_compressed(a.out / stack, frames=np.moveaxis(frames, 0, 2), t=t, clock="video", opencv=cv2.__version__)
        gsds, periods = window_gsd(field)
        for i, s in enumerate(range(1, len(frames) - WINDOW + 2, WINDOW)):
            gsd = "" if static or gsds[i] is None else f"{gsds[i]:.6f}"
            rows.append([stack, "average", WINDOW, s, static, SIZE, gsd, note])
            angle = np.mean([line_angle(f.astype(float)) for f in field[s - 1:s - 1 + WINDOW]])
            lines.append([f"{stack}#{s}-{s + WINDOW - 1}", f"{periods[i]:.2f}", gsd, f"{angle:.3f}"])
        print(f"{video}: {len(frames)} frames -> {stack}")
    with open(a.out / "meta.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["file", "mode", "frames", "first_frame", "static", "roi", "gsd_mm_per_px", "note"])
        w.writerows(rows)
    with open(a.out / "grid_lines.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["file", "period_px", "gsd_mm_per_px", "line_angle_deg"])
        w.writerows(lines)
    print(f"meta.csv and grid_lines.csv: {len(rows)} windows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
