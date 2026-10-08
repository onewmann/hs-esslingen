"""Velocity of a moving texture from a single pulsed-exposure image (2-D FFT + Radon)."""

from .core import (
    Estimate,
    crop_square,
    estimate_velocity,
    fit_speed,
    hann_window,
    log_spectrum,
    mean_projection,
    ripple_profile,
)
from .reference import phase_correlation
from .simulate import make_texture, simulate_pulsed_image

__version__ = "1.0.0"

__all__ = [
    "Estimate",
    "crop_square",
    "estimate_velocity",
    "fit_speed",
    "hann_window",
    "log_spectrum",
    "make_texture",
    "mean_projection",
    "phase_correlation",
    "ripple_profile",
    "simulate_pulsed_image",
]
