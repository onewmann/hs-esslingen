"""Velocity of a moving texture from a single pulsed-exposure image (2-D FFT + Radon).

Public API (the names in ``__all__``):

* :func:`estimate_velocity` - direction and speed from one pulsed image;
  returns an :class:`Estimate`. All options are keyword-only.
* :func:`to_physical` - an estimate in pixels per step converted to length
  units per second, given the ground sample distance and the time step.
* :func:`pulsed_from_stack` - pulsed image from a stack of single frames.
* :func:`angle_difference` - difference of two directions modulo 180 deg.
* :func:`simulate_pulsed_image`, :func:`make_texture` - synthetic images
  with known velocity.
* :func:`phase_correlation` - frame-to-frame shift, an independent reference.

The steps of the estimator (``log_spectrum``, ``mean_projection``,
``ripple_profile``, ``fit_speed``, ``dirichlet_template``, ...) and its
constants can be imported from :mod:`fftvel.core`; they are not part of the
public API and may change. :mod:`fftvel.stream` estimates from cameras and
video files, :mod:`fftvel.plotting` draws the diagnostic figure.
"""

from importlib import metadata as _metadata

from .core import Estimate, angle_difference, estimate_velocity, pulsed_from_stack, to_physical
from .reference import phase_correlation
from .simulate import make_texture, simulate_pulsed_image

try:
    __version__ = _metadata.version("fftvel")
except _metadata.PackageNotFoundError:     # used from the source tree without installing
    __version__ = "0.1.0"

__all__ = [
    "Estimate",
    "angle_difference",
    "estimate_velocity",
    "make_texture",
    "phase_correlation",
    "pulsed_from_stack",
    "simulate_pulsed_image",
    "to_physical",
]
