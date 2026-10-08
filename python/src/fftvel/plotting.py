"""Figure with every stage of the estimate (needs matplotlib)."""

from __future__ import annotations

import numpy as np

from .core import Estimate, crop_square

__all__ = ["plot_estimate"]


def plot_estimate(bp, est: Estimate, frame=None, fig=None):
    """Show pulsed image, log spectrum, ripple profile and the speed fit.

    ``est`` must come from ``estimate_velocity(..., diagnostics=True)``.
    Returns the matplotlib figure.
    """
    import matplotlib.pyplot as plt

    if est.diag is None:
        raise ValueError("call estimate_velocity(..., diagnostics=True)")
    d = est.diag
    if fig is None:
        fig = plt.figure(figsize=(12, 7), constrained_layout=True)
    gs = fig.add_gridspec(2, 3)

    if frame is not None:
        ax = fig.add_subplot(gs[0, 0])
        ax.imshow(crop_square(np.asarray(frame, float)), cmap="gray")
        ax.set_title("Single exposure")
        ax.axis("off")
        ax_bp, ax_spec, ax_txt = fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[1, 2])
        ax_wr, ax_prof = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    else:
        ax_bp, ax_spec, ax_txt = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])
        ax_wr, ax_prof = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1:])

    ax_bp.imshow(crop_square(np.asarray(bp, float)), cmap="gray")
    ax_bp.set_title(f"Pulsed image ({2 * est.L + 1} exposures)")
    ax_bp.axis("off")

    S = d["spectrum"]
    npad = S.shape[0]
    half = (npad // 2) / npad
    ext = [-half, (npad - 1 - npad // 2) / npad, -(npad - 1 - npad // 2) / npad, half]
    ax_spec.imshow(S, cmap="gray", extent=[ext[0], ext[1], ext[2], ext[3]], origin="upper")
    t = np.deg2rad(est.angle_deg)
    ax_spec.plot([-0.5 * np.cos(t), 0.5 * np.cos(t)], [-0.5 * np.sin(t), 0.5 * np.sin(t)], "r-", lw=1)
    ax_spec.set_xlim(-0.5, 0.5)
    ax_spec.set_ylim(-0.5, 0.5)
    ax_spec.set_xlabel("$f_x$ (cycles/px)")
    ax_spec.set_ylabel("$f_y$ (cycles/px)")
    ax_spec.set_title("Log spectrum, direction of v")

    ax_wr.plot(d["theta"], d["wr"], "k-", lw=1)
    ax_wr.axvline(est.angle_deg, color="r", ls="--", lw=1)
    ax_wr.set_xlim(0, 180)
    ax_wr.set_xlabel(r"$\theta$ (deg)")
    ax_wr.set_ylabel(r"WR($\theta$)")
    ax_wr.set_title(f"Ripple profile, GEBA = {est.geba:.2f}")

    ax_prof.plot(d["profile_j"], d["profile"], "k-", lw=1, label="measured")
    ax_prof.plot(d["profile_j"], d["model"], "r-", lw=1, label="Dirichlet model")
    ax_prof.set_xlabel("distance from centre (bins)")
    ax_prof.set_ylabel("mean log spectrum")
    ax_prof.legend(loc="best", frameon=False)
    ax_prof.set_title(f"Projection along v, fit quality {est.quality:.2f}")

    state = "valid" if est.valid else "not valid"
    lines = [
        f"direction  {est.angle_deg:6.1f} deg",
        f"speed      {est.speed:6.2f} px/step",
        f"v (image)  [{est.vx:.2f}, {est.vy:.2f}]",
        f"GEBA       {est.geba:6.2f}",
        f"quality    {est.quality:6.2f} ({state})",
    ]
    ax_txt.axis("off")
    ax_txt.text(0.0, 0.5, "\n".join(lines), family="monospace", fontsize=11, va="center")
    return fig
