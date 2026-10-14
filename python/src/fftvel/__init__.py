"""Surface velocity from averaged camera frames (2-D FFT and Radon transform).

A camera takes ``K`` frames at a fixed frame rate while a textured surface
moves; their mean (:func:`pulsed_from_stack`) shows the surface ``K`` times,
shifted by equal steps. A strobe flashing ``K`` times during one exposure
gives the same image in a single frame, a multi-flash image, and the same
code measures it. The method comes from J. Horn's doctoral thesis (2006, see
README.md), and like the thesis the function names call such an image a
pulsed-exposure image. The option ``exposures`` gives ``K``: the number of
averaged frames, or of flashes. ``L`` is the half-count of the thesis,
``K = 2L+1``.

Public API (the names in ``__all__``):

* :func:`pulsed_from_stack` - the mean of a burst of frames as one pulsed
  image.
* :func:`estimate_velocity` - direction and speed from one pulsed image (the
  mean of ``K`` frames, or a multi-flash image); returns an
  :class:`Estimate`. All options are keyword-only.
* :func:`to_physical` - an estimate in pixels per step converted to length
  units per second, given the ground sample distance and the time step
  (frame interval or flash period); returns a :class:`PhysicalVelocity`.
* :func:`angle_difference` - difference of two directions modulo 180 deg.
* :func:`simulate_pulsed_image`, :func:`make_texture` - synthetic images
  with known velocity.
* :func:`phase_correlation` - frame-to-frame shift, a check that is
  independent of the estimator.
* :func:`plot_estimate` - figure with every stage of an estimate (needs
  matplotlib, which is imported only when it is called).

:mod:`fftvel.io` reads images and frame stacks from files (``load``,
``frames_axis``, ``to_grey``). The steps of the estimator (``log_spectrum``,
``mean_projection``, ``ripple_profile``, ``fit_speed``,
``dirichlet_template``, ...) and its constants can be imported from
:mod:`fftvel.core`; they are not part of the public API and may change.
:mod:`fftvel.stream` estimates from cameras and video files.

>>> import fftvel
>>> bp = fftvel.simulate_pulsed_image(256, (5, -3), 4, rng=0)   # L = 4, so K = 2L+1 = 9 exposures
>>> est = fftvel.estimate_velocity(bp, exposures=9)
>>> round(est.speed, 1), est.valid
(5.8, True)
"""

from importlib import metadata as _metadata

from .core import Estimate, PhysicalVelocity, angle_difference, estimate_velocity, pulsed_from_stack, to_physical
from .phasecorr import phase_correlation
from .plotting import plot_estimate
from .simulate import make_texture, simulate_pulsed_image

try:
    __version__ = _metadata.version("fftvel")
except _metadata.PackageNotFoundError:     # used from the source tree without installing
    __version__ = "0.1.0"

__all__ = [
    "Estimate",
    "PhysicalVelocity",
    "angle_difference",
    "estimate_velocity",
    "make_texture",
    "phase_correlation",
    "plot_estimate",
    "pulsed_from_stack",
    "simulate_pulsed_image",
    "to_physical",
]
