"""Test sequences with known motion (data/testsequences.mat): every window of
9 frames, averaged into one pulsed image, against phase correlation.

    python experiments/test_sequences.py [--quick] [--out DIR] [--replot]

does the same as ``python experiments/run.py --only testsequences ...`` (see
run.py). Writes testsequences.csv and testsequences.png. ``--quick`` only
selects the output folder; this study has no smaller variant. This is a
script, not a test module: no name in it starts with "test".
"""

from __future__ import annotations

import numpy as np

import common
from common import BLUE, INK2, ORANGE, SURFACE
from fftvel import estimate_velocity, phase_correlation

CSV = ["testsequences.csv"]

SEQUENCES = (("model", "model: (5, 5) px/frame"),
             ("tyre", "tyre: (4, 4) px/frame"),
             ("accelerating", "accelerating: vx 2-8, vy 3"),
             ("vehicle", "vehicle: irregular steps 3.5-9.9 px"))


def _zero_border(frame):
    """Rows/columns of zero fill at the top and left (the sequences move down and to the right)."""
    rows = int(np.argmax(frame.any(axis=1)))
    cols = int(np.argmax(frame.any(axis=0)))
    return max(rows, cols)


def run():
    from scipy.io import loadmat

    data = loadmat(common.REPO / "data" / "testsequences.mat")
    rows = []
    for name, _ in SEQUENCES:
        s = data[name].astype(float) / 255
        nf = s.shape[2]
        shift = np.array([phase_correlation(s[:, :, f], s[:, :, f + 1])[:2] for f in range(nf - 1)])
        for st in range(nf - 8):
            last = st + 8
            border = max(_zero_border(s[:, :, k]) for k in range(st, last + 1))
            e = estimate_velocity(s[border:, border:, st:last + 1].mean(axis=2), L=4)
            r = float(np.hypot(*shift[st:last].mean(axis=0)))
            rows.append(dict(sequence=name, first_frame=st + 1, est_speed=e.speed, ref_speed=r,
                             est_angle=e.angle_deg, quality=e.quality, quality_perp=e.quality_perp,
                             isotropic=int(e.isotropic), valid=int(e.valid), quality_blur=e.quality_blur,
                             blur_length=e.blur_length, blur_like=int(e.blur_like), saturated=e.saturated))
    common.write_csv("testsequences.csv", rows)
    plot(rows)
    return rows


def plot(rows):
    fig, axes = common.figure(4, 14.5)
    for ax, (name, title) in zip(axes, SEQUENCES):
        sel = [r for r in rows if r["sequence"] == name]
        first = np.array([r["first_frame"] for r in sel])
        est = np.array([r["est_speed"] for r in sel])
        ref = np.array([r["ref_speed"] for r in sel])
        valid = np.array([bool(r["valid"]) for r in sel])
        ax.plot(first, ref, color=INK2, lw=2, label="phase correlation (reference)")
        ax.plot(first[valid], est[valid], "o", color=BLUE, ms=6, mec=SURFACE, mew=1.2,
                label="single-image estimate")
        if (~valid).any():
            ax.plot(first[~valid], est[~valid], "o", color=ORANGE, ms=6, mfc="none", label="not valid")
        ax.set_ylim(0, 16.5)
        ax.set_title(title, color=INK2, fontsize=10, loc="left")
        common.style(ax, "first frame of the 9-frame window", "speed (px/frame)" if ax is axes[0] else "")
    axes[0].legend(frameon=False, fontsize=9, labelcolor=INK2, loc="lower left")
    common.title(fig, "Test sequences with known motion",
                 "every window of 9 frames, cropped to the area without zero-filled borders; "
                 "the method assumes equal steps, which the last two sequences break")
    common.save(fig, "testsequences.png")


def summary(rows):
    for name, _ in SEQUENCES:
        sel = [r for r in rows if r["sequence"] == name]
        val = [r for r in sel if r["valid"]]
        err = [abs(r["est_speed"] / r["ref_speed"] - 1) for r in val]
        ang = [r["est_angle"] for r in val]
        blur = f", {int(sum(r['blur_like'] for r in val))} of them blur-like" if sel and "blur_like" in sel[0] else ""
        print(f"{name:13s} {len(val)}/{len(sel)} windows valid{blur}, max relative speed error of valid "
              f"{100 * max(err) if err else float('nan'):.1f} %, angles "
              + (f"{min(ang):.2f}..{max(ang):.2f}" if ang else "-"))


def replot():
    rows = common.read_csv("testsequences.csv")
    plot(rows)
    summary(rows)


def main(argv=None):
    import run as cli
    cli.main(argv, only=["testsequences"], description=__doc__)


if __name__ == "__main__":
    main()
