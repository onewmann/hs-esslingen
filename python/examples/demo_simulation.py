"""Simulate a pulsed-exposure image with known motion and estimate it.

    python python/examples/demo_simulation.py [--random] [--save fig.png]
"""

import argparse
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

from fftvel import estimate_velocity, simulate_pulsed_image
from fftvel.plotting import plot_estimate

REPO = Path(__file__).resolve().parents[2]

p = argparse.ArgumentParser()
p.add_argument("--random", action="store_true", help="random texture instead of the camera texture")
p.add_argument("--save", help="write the figure to this file instead of showing it")
a = p.parse_args()

L = 4                       # 2L+1 = 9 exposures
v_true = (5.0, -3.0)        # px per exposure step, x right, y down
texture = None if a.random else mpimg.imread(REPO / "data" / "surface_texture.png").astype(float)

bp, frames = simulate_pulsed_image(256, v_true, L, texture=texture, noise=0.1, rng=1,
                                   return_frames=True)
est = estimate_velocity(bp, L=L, diagnostics=True)

true_angle = np.rad2deg(np.arctan2(-v_true[1], v_true[0])) % 180
print(f"true:      direction {true_angle:6.2f} deg, speed {np.hypot(*v_true):6.3f} px/step")
print(f"estimate:  direction {est.angle_deg:6.2f} deg, speed {est.speed:6.3f} px/step "
      f"(quality {est.quality:.2f}, GEBA {est.geba:.2f})")

fig = plot_estimate(bp, est, frame=frames[:, :, L])
if a.save:
    fig.savefig(a.save, dpi=120)
else:
    plt.show()
