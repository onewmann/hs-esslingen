"""Velocity of a moving texture from a single pulsed-exposure image (2-D FFT + Radon)."""

from importlib import metadata as _metadata

from .core import (
    Estimate,
    angle_difference,
    crop_square,
    estimate_velocity,
    fit_speed,
    hann_window,
    log_spectrum,
    mean_projection,
    pulsed_from_stack,
    ripple_profile,
)
from .reference import phase_correlation
from .simulate import make_texture, simulate_pulsed_image

try:
    __version__ = _metadata.version("fftvel")
except _metadata.PackageNotFoundError:     # used from the source tree without installing
    __version__ = "0.1.0"

__all__ = [
    "Estimate",
    "angle_difference",
    "crop_square",
    "estimate_velocity",
    "fit_speed",
    "hann_window",
    "log_spectrum",
    "make_texture",
    "mean_projection",
    "phase_correlation",
    "pulsed_from_stack",
    "ripple_profile",
    "simulate_pulsed_image",
]
