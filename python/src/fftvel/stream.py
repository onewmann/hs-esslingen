"""Estimates from a stream of frames (camera or video file).

A pulsed image is the mean of ``2L+1`` consecutive frames. The model needs
equal time steps between them, so every window is checked against its
timestamps and flagged when the frame intervals differ by more than ``tol``.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Iterable, Iterator, Optional, Tuple

import numpy as np

from .core import crop_square, estimate_velocity

__all__ = ["FrameWindow", "StreamEstimate", "estimate_stream", "estimate_window", "to_grey"]


def to_grey(frame: np.ndarray) -> np.ndarray:
    """Grey-scale float image; colour input uses the MATLAB ``rgb2gray`` weights (RGB order)."""
    a = np.asarray(frame, dtype=float)
    if a.ndim == 3 and a.shape[2] == 1:
        a = a[..., 0]
    elif a.ndim == 3:
        a = a[..., :3] @ np.array([0.2989, 0.5870, 0.1140])
    return a


class FrameWindow:
    """Ring buffer of the last ``2L+1`` frames with their timestamps."""

    def __init__(self, L: int = 4, roi: Optional[int] = None, tol: float = 0.2):
        self.L = L
        self.roi = roi
        self.tol = tol
        self._frames: deque = deque(maxlen=2 * L + 1)
        self._times: deque = deque(maxlen=2 * L + 1)

    def push(self, frame: np.ndarray, t: float) -> None:
        img = to_grey(frame)
        img = crop_square(img, None if self.roi is None else min(self.roi, *img.shape))
        self._frames.append(img)
        self._times.append(float(t))

    @property
    def ready(self) -> bool:
        return len(self._frames) == self._frames.maxlen

    def pulsed_image(self) -> Tuple[np.ndarray, float, bool]:
        """Return ``(bp, step, regular)``: mean image, median frame interval, equal-steps flag."""
        if not self.ready:
            raise RuntimeError("window not full yet")
        bp = np.mean(np.stack(self._frames, axis=2), axis=2)
        dt = np.diff(np.asarray(self._times))
        step = float(np.median(dt))
        regular = bool(step > 0 and np.max(np.abs(dt - step)) <= self.tol * step)
        return bp, step, regular


@dataclass
class StreamEstimate:
    t: float               # timestamp of the last frame in the window
    angle_deg: float
    speed: float           # px per frame
    speed_per_s: float     # px per second (nan if the step is unknown)
    quality: float
    valid: bool            # fit quality reached and frame steps equal
    regular: bool          # frame steps equal within the tolerance


def estimate_window(frames: Iterable[Tuple[np.ndarray, float]], L: int = 4,
                    roi: Optional[int] = None, tol: float = 0.2,
                    min_quality: float = 0.6) -> StreamEstimate:
    """Estimate from exactly ``2L+1`` consecutive ``(frame, t)`` pairs (e.g. one camera burst)."""
    win = FrameWindow(L, roi, tol)
    t = float("nan")
    for frame, t in frames:
        win.push(frame, t)
    return _estimate(win, t, min_quality)


def _estimate(win: FrameWindow, t: float, min_quality: float) -> StreamEstimate:
    bp, step, regular = win.pulsed_image()
    e = estimate_velocity(bp, L=win.L, min_quality=min_quality)
    return StreamEstimate(
        t=float(t), angle_deg=e.angle_deg, speed=e.speed,
        speed_per_s=e.speed / step if step > 0 else float("nan"),
        quality=e.quality, valid=e.valid and regular, regular=regular)


def estimate_stream(frames: Iterable[Tuple[np.ndarray, float]], L: int = 4,
                    roi: Optional[int] = None, every: Optional[int] = None,
                    tol: float = 0.2, min_quality: float = 0.6) -> Iterator[StreamEstimate]:
    """Yield an estimate for every ``every`` new frames (default ``2L+1``, no overlap).

    Meant for sources that deliver every frame, such as video files. ``frames``
    yields ``(frame, t)`` pairs with ``t`` in seconds. ``roi`` is the side of
    the central square that is analysed (default: largest square). Windows
    with unequal frame intervals are still estimated but reported with
    ``regular=False`` and ``valid=False``. For live cameras, grab a burst of
    ``2L+1`` frames per estimate and call :func:`estimate_window`, so that the
    time spent estimating does not mix old buffered frames into a window.
    """
    win = FrameWindow(L, roi, tol)
    every = every or (2 * L + 1)
    since = 0
    for frame, t in frames:
        win.push(frame, t)
        since += 1
        if win.ready and since >= every:
            since = 0
            yield _estimate(win, t, min_quality)
