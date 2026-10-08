"""Sliding 9-frame windows over the test sequences in data/testsequences.mat.

Each frame is the previous one shifted with zero fill, so every window is
cropped to the area that holds image content in all of its frames. The
reference is the mean frame-to-frame shift from phase correlation.

    python python/examples/demo_testsequences.py [--save fig.png]
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

from fftvel import estimate_velocity, phase_correlation

REPO = Path(__file__).resolve().parents[2]
L = 4

p = argparse.ArgumentParser()
p.add_argument("--save")
a = p.parse_args()

data = loadmat(REPO / "data" / "testsequences.mat")
names = ["model", "tyre", "accelerating"]
fig, axes = plt.subplots(1, len(names), figsize=(13, 3.8), constrained_layout=True)
for ax, name in zip(axes, names):
    s = data[name].astype(float) / 255
    nf = s.shape[2]
    shift = np.array([phase_correlation(s[:, :, f], s[:, :, f + 1])[:2] for f in range(nf - 1)])
    starts = np.arange(nf - 2 * L)
    est, ref, valid = [], [], []
    for st in starts:
        last = st + 2 * L
        border = int(np.ceil(np.abs(shift[:last].sum(axis=0)).max()))   # zero fill so far
        bp = s[border:, border:, st:last + 1].mean(axis=2)
        e = estimate_velocity(bp, L=L)
        est.append(e.speed)
        valid.append(e.valid)
        ref.append(np.hypot(*shift[st:last].mean(axis=0)))
    est, ref, valid = np.array(est), np.array(ref), np.array(valid)
    ax.plot(starts + 1, ref, "k-", lw=1.5, label="phase correlation")
    ax.plot(starts[valid] + 1, est[valid], "o", color="#2a78d6", label="estimate (valid)")
    ax.plot(starts[~valid] + 1, est[~valid], "o", mfc="none", color="#eb6834", label="estimate (not valid)")
    ax.set_ylim(0, 1.5 * ref.max())
    ax.set_title(name)
    ax.set_xlabel("first frame of window")
    ax.set_ylabel("speed (px/frame)")
    err = np.median(np.abs(est[valid] - ref[valid])) if valid.any() else float("nan")
    print(f"{name:13s} {len(starts)} windows, {valid.sum()} valid, median |error| of valid {err:.3f} px/frame")
axes[0].legend(loc="lower left", frameon=False)
if a.save:
    fig.savefig(a.save, dpi=120)
else:
    plt.show()
