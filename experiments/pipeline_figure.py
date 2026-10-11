"""Figure of the processing chain for one simulated image of the camera texture.

    python experiments/pipeline_figure.py [--quick] [--out DIR] [--replot]

does the same as ``python experiments/run.py --only pipeline ...`` (see
run.py). Writes pipeline.png; there is no CSV file, so ``--replot`` computes
the figure again. ``--quick`` only selects the output folder.
"""

from __future__ import annotations

import common
from fftvel import estimate_velocity, simulate_pulsed_image

CSV = []


def run():
    import matplotlib.image as mpimg
    from fftvel.plotting import plot_estimate

    tex = mpimg.imread(common.REPO / "data" / "surface_texture.png").astype(float)
    bp, frames = simulate_pulsed_image(common.N, (5, -3), 4, texture=tex, rng=3, return_frames=True)
    e = estimate_velocity(bp, L=4, diagnostics=True)
    fig = plot_estimate(bp, e, frame=frames[:, :, 4])
    fig.patch.set_facecolor("white")
    common.save(fig, "pipeline.png", dpi=110, facecolor="white")


def replot():
    run()


def main(argv=None):
    import run as cli
    cli.main(argv, only=["pipeline"], description=__doc__)


if __name__ == "__main__":
    main()
