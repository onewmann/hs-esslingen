"""Figure of the processing chain for one simulated image of the camera texture,
and the short version of it for the top of the README.

    python experiments/pipeline_figure.py [--quick] [--out DIR] [--replot]

does the same as ``python experiments/run.py --only pipeline ...`` (see
run.py). Writes pipeline.png (every stage of the estimate) and hero.png (one
frame, the mean of 9 frames of the moving surface and its log spectrum with
the direction and the line spacing). The simulated image is the mean of 9
frames with short exposures, which is equally one exposure with 9 short
flashes. There is no CSV file, so ``--replot`` computes both figures again;
the image is fixed by its seed, so every run draws the same figures.
``--quick`` only selects the output folder.
"""

from __future__ import annotations

import math

import numpy as np

import common
from common import INK, INK2, ORANGE, SURFACE
from fftvel import estimate_velocity, simulate_pulsed_image

CSV = []

V_TRUE = (5, -3)        # px per frame, x right and y down: 5.83 px at 30.96 deg
L = 4                   # 2L + 1 = 9 frames
SEED = 3
SPECTRUM_HALF = 0.3     # hero.png shows |f_x|, |f_y| <= 0.3 cycles/px of the log spectrum
BOX = dict(boxstyle="round,pad=0.3", fc=INK, ec="none", alpha=0.62)    # label boxes on the images


def run():
    import matplotlib.image as mpimg
    from fftvel.plotting import plot_estimate

    tex = mpimg.imread(common.REPO / "data" / "surface_texture.png").astype(float)
    bp, frames = simulate_pulsed_image(common.N, V_TRUE, L, texture=tex, rng=SEED, return_frames=True)
    e = estimate_velocity(bp, L=L, diagnostics=True)
    fig = plot_estimate(bp, e, frame=frames[:, :, L])
    for ax in fig.axes:                     # panel titles in the terms of averaged frames
        if ax.get_title() == "Single exposure":
            ax.set_title("One frame")
        elif ax.get_title().startswith("Pulsed image"):
            ax.set_title(f"{e.exposures} frames averaged")
    fig.patch.set_facecolor("white")
    common.save(fig, "pipeline.png", dpi=110, facecolor="white")
    hero(frames[:, :, L], bp, e)


def _spectrum_extent(npad):
    """Extent of the fftshift-ed spectrum in cycles/px, f_y pointing up (as in fftvel.plotting)."""
    half = (npad // 2) / npad
    edge = (npad - 1 - npad // 2) / npad
    return [-half, edge, -edge, half]


def hero(frame, bp, e):
    """hero.png: one frame, the mean of 2L+1 frames of the moving surface and its log spectrum, in one row."""
    import matplotlib.pyplot as plt
    from matplotlib import patheffects

    true_angle = math.degrees(math.atan2(-V_TRUE[1], V_TRUE[0])) % 180
    true_speed = math.hypot(*V_TRUE)
    k = 2 * L + 1

    fig, axes = plt.subplots(1, 3, figsize=(13, 5.0), layout="constrained")
    fig.get_layout_engine().set(rect=(0, 0.07, 1, 0.93), w_pad=0.15)    # room for the caption line
    fig.patch.set_facecolor(SURFACE)
    titles = ("(a) One frame, real texture",
              f"(b) Moving, {k} frames averaged",
              "(c) Log spectrum of (b)")
    for ax, img, t in zip(axes[:2], (frame, bp), titles):
        ax.imshow(img, cmap="gray", interpolation="none")
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(t, color=INK, fontsize=15, loc="left")

    # motion of the surface during the burst: from the first to the last frame, 2L steps of v
    n = bp.shape[0]
    ax = axes[1]
    x0, y0 = 0.10 * n, 0.88 * n
    halo = [patheffects.withStroke(linewidth=4.5, foreground=INK)]
    ax.annotate("", xy=(x0 + 2 * L * V_TRUE[0], y0 + 2 * L * V_TRUE[1]), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color="white", lw=2.2, mutation_scale=18,
                                path_effects=halo))
    ax.text(x0 + 2 * L * V_TRUE[0] + 8, y0 + 2 * L * V_TRUE[1] + 6, f"motion over\nthe {k} frames",
            color="white", fontsize=13, va="center", bbox=BOX)

    ax = axes[2]
    S = e.diag["spectrum"]
    ext = _spectrum_extent(S.shape[0])
    fx = np.linspace(ext[0], ext[1], S.shape[1])
    shown = S[np.ix_(np.abs(fx) <= SPECTRUM_HALF, np.abs(fx) <= SPECTRUM_HALF)]
    ax.imshow(S, cmap="gray", extent=ext, origin="upper", interpolation="bilinear",
              vmin=np.percentile(shown, 3), vmax=np.percentile(shown, 99.8))
    ax.set_xlim(-SPECTRUM_HALF, SPECTRUM_HALF)
    ax.set_ylim(-SPECTRUM_HALF, SPECTRUM_HALF)
    ax.set_xticks([-0.2, 0, 0.2])
    ax.set_yticks([-0.2, 0, 0.2])
    ax.tick_params(colors=INK2, labelsize=12)
    for s in ax.spines.values():
        s.set_color(INK2)
    ax.set_xlabel("$f_x$ (cycles/px)", color=INK2, fontsize=13)
    ax.set_ylabel("$f_y$ (cycles/px)", color=INK2, fontsize=13, labelpad=0)
    ax.set_title(titles[2], color=INK, fontsize=15, loc="left")

    # direction found: an arrow from the centre, perpendicular to the lines
    t = math.radians(e.angle_deg)
    u = np.array([math.cos(t), math.sin(t)])          # along v
    p = np.array([-math.sin(t), math.cos(t)])         # along the lines
    tip = 0.27 * u
    ax.annotate("", xy=tip, xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=2.4, mutation_scale=20))
    ax.text(*(tip + 0.07 * p - 0.05 * u), f"direction {e.angle_deg:.2f}°", color="white", fontsize=13,
            ha="center", va="center", bbox=dict(BOX, fc=ORANGE, alpha=0.9))

    # line spacing: the bright lines lie 1/d apart along v; mark the first gap beside the arrow
    gap = 1 / e.speed
    off = -0.075 * p
    a0, a1 = off, off + gap * u
    for q in (np.zeros(2), gap * u):                  # short marks on the centre line and the next one
        seg = np.array([q - 0.045 * p, q + 0.045 * p]) + off
        ax.plot(seg[:, 0], seg[:, 1], color="white", lw=1.6)
    ax.annotate("", xy=a1, xytext=a0,
                arrowprops=dict(arrowstyle="<|-|>", color="white", lw=1.6, mutation_scale=12,
                                shrinkA=0, shrinkB=0))
    mid = 0.5 * (a0 + a1) - 0.1 * p
    ax.text(*mid, f"line spacing\n1 / {e.speed:.2f} px", color="white", fontsize=13, ha="center", va="center",
            bbox=BOX)

    fig.text(0.5, 0.03, f"In simulation, real camera texture, {n} × {n} px: estimate {e.angle_deg:.2f}°, {e.speed:.3f} px per frame; "
             f"true {true_angle:.2f}°, {true_speed:.3f} px per frame",
             ha="center", va="center", color=INK, fontsize=14)
    common.save(fig, "hero.png", dpi=130)


def replot():
    run()


def main(argv=None):
    import run as cli
    cli.main(argv, only=["pipeline"], description=__doc__)


if __name__ == "__main__":
    main()
