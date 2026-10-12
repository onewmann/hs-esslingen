"""Estimates from a stream of frames (camera or video file).

Two ways to get a pulsed image from a camera:

* Flash mode: a strobe fires ``K = 2L+1`` times (or ``exposures``) during
  one exposure, so every frame is already a pulsed image.
  :func:`estimate_frame` estimates one such frame. The speed is in pixels
  per flash step.
* Average mode: under continuous light, ``K`` consecutive frames are
  averaged (:class:`FrameWindow`, :func:`estimate_window`,
  :func:`estimate_stream`). The model needs equal time steps between them,
  so every window is checked against its timestamps (:func:`frame_step`)
  and flagged when the spread of the frame intervals exceeds ``tol`` times
  their median. The speed is in pixels per frame.
"""

from __future__ import annotations

import warnings
from collections import deque
from dataclasses import dataclass
from typing import Iterable, Iterator, Optional, Tuple

import numpy as np

from .core import (DEFAULT_L, DEFAULT_MIN_QUALITY, GREY_WEIGHTS, SAT_WARN, Estimate, _check_option,
                   crop_square, estimate_velocity)

__all__ = ["FrameWindow", "StreamEstimate", "estimate_frame", "estimate_stream", "estimate_window",
           "frame_step", "to_grey"]


def to_grey(frame: np.ndarray) -> np.ndarray:
    """Grey-scale float image; colour input uses the MATLAB ``rgb2gray`` weights (RGB order)."""
    a = np.asarray(frame, dtype=float)
    if a.ndim == 3 and a.shape[2] == 1:
        a = a[..., 0]
    elif a.ndim == 3:
        a = a[..., :3] @ np.array(GREY_WEIGHTS)
    return a


def _check_number(x, name: str, lowest: float, inclusive: bool):
    """Return ``x`` as float if it is a finite real number >= ``lowest`` (``inclusive``) or > ``lowest``."""
    ok = (not isinstance(x, (str, bool)) and np.isscalar(x) and np.isreal(x) and np.isfinite(x)
          and (x >= lowest if inclusive else x > lowest))
    if not ok:
        what = f">= {lowest:g}" if inclusive else f"> {lowest:g}"
        raise ValueError(f"{name} must be a finite number {what}, got {x!r}")
    return float(x)


def _check_seconds(x: Optional[float], name: str) -> Optional[float]:
    return None if x is None else _check_number(x, name, 0, False)


def frame_step(times, tol: float = 0.1, period: Optional[float] = None) -> Tuple[float, bool]:
    """Time step of a window of frames and whether the frames are equally spaced.

    Returns ``(step, regular)``. Without ``period`` the timestamps are taken
    as exact (camera clock): ``step`` is the median frame interval, and the
    window is regular if ``step > 0`` and the spread of the intervals,
    ``max - min``, is at most ``tol * step``.

    ``period`` is the nominal frame period reported by the camera, for
    timestamps taken on the host, which carry transfer jitter. If the mean
    interval agrees with ``period`` within ``tol``, ``step`` is ``period``
    and the window is regular unless an interval is off by half a period or
    more (a dropped frame). Otherwise the spread test above applies.
    """
    tol = _check_number(tol, "tol", 0, True)
    period = _check_seconds(period, "period")
    dt = np.diff(np.asarray(times, dtype=float))
    if dt.size == 0:
        return float("nan"), False
    if period is not None and abs(np.mean(dt) / period - 1) <= tol:
        return period, bool(np.all(np.abs(dt / period - 1) < 0.5))
    step = float(np.median(dt))
    return step, bool(step > 0 and np.ptp(dt) <= tol * step)


def _window_length(L: int, exposures: Optional[int]) -> int:
    L = _check_option(L, "L", "positive")
    if exposures is not None:
        return _check_option(exposures, "exposures", "exposures")
    return 2 * L + 1


class FrameWindow:
    """Ring buffer of the last ``K`` frames with their timestamps.

    ``K`` is ``exposures`` if given, else ``2L+1``. Each frame is converted
    to grey and cropped to its central ``roi``-by-``roi`` square (default:
    the largest square) when it is pushed. ``tol`` and ``period`` are
    passed to :func:`frame_step`.
    """

    def __init__(self, L: int = DEFAULT_L, roi: Optional[int] = None, tol: float = 0.1, *,
                 exposures: Optional[int] = None, period: Optional[float] = None):
        n = _window_length(L, exposures)
        self.L = L
        self.exposures = exposures
        self.roi = None if roi is None else _check_option(roi, "roi", "roi")
        self.tol = _check_number(tol, "tol", 0, True)
        self.period = _check_seconds(period, "period")
        self.size = n
        self._frames: deque = deque(maxlen=n)
        self._times: deque = deque(maxlen=n)

    def push(self, frame: np.ndarray, t: float) -> None:
        img = to_grey(frame)
        img = crop_square(img, None if self.roi is None else min(self.roi, *img.shape))
        self._frames.append(img)
        self._times.append(float(t))

    @property
    def ready(self) -> bool:
        return len(self._frames) == self._frames.maxlen

    def pulsed_image(self) -> Tuple[np.ndarray, float, bool]:
        """Return ``(bp, step, regular)``: mean image, frame step, equal-steps flag (see :func:`frame_step`)."""
        if not self.ready:
            raise RuntimeError("window not full yet")
        bp = np.mean(np.stack(self._frames, axis=2), axis=2)
        step, regular = frame_step(self._times, self.tol, self.period)
        return bp, step, regular


@dataclass
class StreamEstimate:
    t: float               # timestamp of the (last) frame
    angle_deg: float
    speed: float           # px per step: per frame (average mode) or per flash (flash mode)
    speed_per_s: float     # px per second (nan if the step is unknown)
    quality: float
    valid: bool            # estimate valid and frame steps equal
    regular: bool          # frame steps equal within the tolerance (always True in flash mode)
    step: float = float("nan")              # time step in s: frame step or flash period (nan if unknown)
    estimate: Optional[Estimate] = None     # full result of estimate_velocity
    warning: str = ""                       # saturated pixels (see _saturation_warning); "" if none


def _saturation_warning(e: Estimate) -> str:
    """Text of the warning for more than SAT_WARN saturated pixels, or ''; also issued once as a UserWarning.

    Python shows a UserWarning with the same text only once, so a live loop
    is not flooded; the share of every estimate is in ``warning`` and in
    ``estimate.saturated``.
    """
    if e.saturated <= SAT_WARN:
        return ""
    warnings.warn(f"more than {100 * SAT_WARN:g} % of the analysed pixels are saturated; saturated areas and "
                  f"glints can give wrong results that are marked valid", stacklevel=3)
    return f"{100 * e.saturated:.1f} % of the analysed pixels are saturated"


def estimate_frame(frame: np.ndarray, t: float = float("nan"), L: int = DEFAULT_L,
                   roi: Optional[int] = None, flash_period: Optional[float] = None,
                   min_quality: float = DEFAULT_MIN_QUALITY, *, exposures: Optional[int] = None,
                   duty: float = 0.0, clip_outliers: bool = False) -> StreamEstimate:
    """Estimate from one multi-flash frame (flash mode).

    ``frame`` holds ``K = 2L+1`` (or ``exposures``) strobe flashes within one
    exposure. It is converted to grey and cropped to its central
    ``roi``-by-``roi`` square (default: the largest square). ``flash_period``
    is the time between flashes in seconds; without it ``speed_per_s`` is
    nan. ``duty`` is the flash pulse length divided by the flash period. The
    flash timing cannot be checked from frame timestamps, so ``regular`` is
    always True.
    """
    flash_period = _check_seconds(flash_period, "flash_period")
    img = to_grey(frame)
    img = crop_square(img, None if roi is None else min(_check_option(roi, "roi", "roi"), *img.shape))
    e = estimate_velocity(img, L=L, min_quality=min_quality, exposures=exposures, duty=duty,
                          clip_outliers=clip_outliers)
    step = float("nan") if flash_period is None else flash_period
    return StreamEstimate(
        t=float(t), angle_deg=e.angle_deg, speed=e.speed, speed_per_s=e.speed / step,
        quality=e.quality, valid=e.valid, regular=True, step=step, estimate=e, warning=_saturation_warning(e))


def estimate_window(frames: Iterable[Tuple[np.ndarray, float]], L: int = DEFAULT_L,
                    roi: Optional[int] = None, tol: float = 0.1,
                    min_quality: float = DEFAULT_MIN_QUALITY, *, exposures: Optional[int] = None,
                    duty: float = 0.0, clip_outliers: bool = False,
                    period: Optional[float] = None) -> StreamEstimate:
    """Estimate from exactly ``K = 2L+1`` (or ``exposures``) consecutive ``(frame, t)`` pairs.

    Average mode, e.g. one camera burst. ``duty`` is the exposure time
    divided by the frame step; ``tol`` and ``period`` are passed to
    :func:`frame_step`.
    """
    win = FrameWindow(L, roi, tol, exposures=exposures, period=period)
    t = float("nan")
    for frame, t in frames:
        win.push(frame, t)
    return _estimate(win, t, min_quality, duty, clip_outliers)


def _estimate(win: FrameWindow, t: float, min_quality: float, duty: float,
              clip_outliers: bool) -> StreamEstimate:
    bp, step, regular = win.pulsed_image()
    e = estimate_velocity(bp, L=win.L, min_quality=min_quality, exposures=win.exposures,
                          duty=duty, clip_outliers=clip_outliers)
    return StreamEstimate(
        t=float(t), angle_deg=e.angle_deg, speed=e.speed,
        speed_per_s=e.speed / step if step > 0 else float("nan"),
        quality=e.quality, valid=e.valid and regular, regular=regular, step=step, estimate=e,
        warning=_saturation_warning(e))


def estimate_stream(frames: Iterable[Tuple[np.ndarray, float]], L: int = DEFAULT_L,
                    roi: Optional[int] = None, every: Optional[int] = None,
                    tol: float = 0.1, min_quality: float = DEFAULT_MIN_QUALITY, *,
                    exposures: Optional[int] = None, duty: float = 0.0,
                    clip_outliers: bool = False,
                    period: Optional[float] = None) -> Iterator[StreamEstimate]:
    """Yield an estimate for every ``every`` new frames (default ``K``, no overlap).

    Average mode, meant for sources that deliver every frame, such as video
    files. ``frames`` yields ``(frame, t)`` pairs with ``t`` in seconds.
    ``K`` consecutive frames (``2L+1``, or ``exposures``) are averaged.
    ``roi`` is the side of the central square that is analysed (default:
    largest square). Windows with unequal frame intervals (see
    :func:`frame_step`) are still estimated but reported with
    ``regular=False`` and ``valid=False``. ``duty`` is the exposure time
    divided by the frame step. For live cameras, grab a burst of ``K``
    frames per estimate and call :func:`estimate_window`, so that the time
    spent estimating does not mix old buffered frames into a window.
    """
    win = FrameWindow(L, roi, tol, exposures=exposures, period=period)
    every = every or win.size
    since = 0
    for frame, t in frames:
        win.push(frame, t)
        since += 1
        if win.ready and since >= every:
            since = 0
            yield _estimate(win, t, min_quality, duty, clip_outliers)
