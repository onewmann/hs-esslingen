"""Velocity from a live camera or a video file.

    python python/examples/live_camera.py --source pylon            # Basler camera (pip install pypylon)
    python python/examples/live_camera.py --source 0                # any webcam through OpenCV
    python python/examples/live_camera.py --source recording.mp4    # video file

Needs the [video] extra for OpenCV (pip install -e "./python[examples]").

Every estimate uses 2L+1 consecutive frames. Live cameras deliver them as one
burst per estimate, so the time spent estimating never mixes old buffered
frames into a window. The pulsed-exposure model needs equal time steps; each
window is checked against its timestamps (camera clock for Basler, the
driver's buffer time or the host clock for webcams, the container time for
files) and reported as not valid when the intervals differ by more than 20 %.
``--show`` opens a live plot.
"""

from __future__ import annotations

import argparse
import time
from typing import Iterator, List, Tuple

import numpy as np

from fftvel.stream import estimate_stream, estimate_window

Frames = List[Tuple[np.ndarray, float]]


def _try_set(cam, name: str, value) -> bool:
    try:
        getattr(cam, name).SetValue(value)
        return True
    except Exception:  # noqa: BLE001 - node missing on this camera model
        return False


def _tick_seconds(cam) -> float:
    """Length of one timestamp tick: GigE models report their tick frequency, USB3 counts ns."""
    try:
        return 1.0 / float(cam.GevTimestampTickFrequency.GetValue())
    except Exception:  # noqa: BLE001
        return 1e-9


def pylon_bursts(n: int, frame_rate: float | None) -> Iterator[Frames]:
    """Bursts of ``n`` consecutive frames with camera timestamps from the first Basler camera."""
    from pypylon import pylon

    cam = pylon.InstantCamera(pylon.TlFactory.GetInstance().CreateFirstDevice())
    cam.Open()
    try:
        _try_set(cam, "PixelFormat", "Mono8")
        if frame_rate:
            ok = _try_set(cam, "AcquisitionFrameRateEnable", True)
            if not (_try_set(cam, "AcquisitionFrameRate", frame_rate)        # USB3, newer GigE
                    or _try_set(cam, "AcquisitionFrameRateAbs", frame_rate)):  # GigE ace classic
                print("warning: could not set the frame rate" + ("" if ok else " (no frame-rate node)"))
        tick = _tick_seconds(cam)
        cam.MaxNumBuffer.SetValue(max(10, n + 2))
        while True:
            cam.StartGrabbingMax(n, pylon.GrabStrategy_OneByOne)
            burst: Frames = []
            while cam.IsGrabbing():
                res = cam.RetrieveResult(5000, pylon.TimeoutHandling_ThrowException)
                try:
                    if res.GrabSucceeded():
                        burst.append((res.Array.copy(), res.TimeStamp * tick))
                finally:
                    res.Release()
            if len(burst) == n:
                yield burst
    finally:
        cam.StopGrabbing()
        cam.Close()


def webcam_bursts(index: int, n: int) -> Iterator[Frames]:
    """Bursts of ``n`` frames from a webcam; frames buffered during an estimate are dropped first."""
    import cv2

    cap = cv2.VideoCapture(index)
    if not cap.isOpened():
        raise SystemExit(f"cannot open camera {index}")
    cap.set(cv2.CAP_PROP_BUFFERSIZE, 1)            # not every backend honours this
    try:
        while True:
            for _ in range(5):                       # drain stale buffers
                cap.grab()
            frames, driver_t, host_t = [], [], []
            for _ in range(n):
                ok, frame = cap.read()
                if not ok:
                    return
                host_t.append(time.monotonic())
                driver_t.append(cap.get(cv2.CAP_PROP_POS_MSEC) * 1e-3)
                frames.append(frame[..., ::-1] if frame.ndim == 3 else frame)   # BGR -> RGB
            d = np.diff(driver_t)
            times = driver_t if min(driver_t) > 0 and np.all(d > 0) else host_t
            yield list(zip(frames, times))
    finally:
        cap.release()


def file_frames(path: str) -> Iterator[Tuple[np.ndarray, float]]:
    """Every frame of a video file with its presentation time."""
    import cv2

    cap = cv2.VideoCapture(path)
    if not cap.isOpened():
        raise SystemExit(f"cannot open {path!r}")
    try:
        while True:
            ok, frame = cap.read()
            if not ok:
                break
            yield (frame[..., ::-1] if frame.ndim == 3 else frame), cap.get(cv2.CAP_PROP_POS_MSEC) * 1e-3
    finally:
        cap.release()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--source", default="pylon", help="'pylon', a webcam index or a video file")
    p.add_argument("--L", type=int, default=4, help="half number of frames per estimate")
    p.add_argument("--roi", type=int, default=512, help="side of the analysed central square (px)")
    p.add_argument("--every", type=int, help="video files: new frames between estimates (default 2L+1)")
    p.add_argument("--frame-rate", type=float, help="set the Basler frame rate (Hz)")
    p.add_argument("--show", action="store_true", help="live plot of speed over time")
    a = p.parse_args()

    n = 2 * a.L + 1
    if a.source == "pylon":
        results = (estimate_window(b, a.L, a.roi) for b in pylon_bursts(n, a.frame_rate))
    elif a.source.isdigit():
        results = (estimate_window(b, a.L, a.roi) for b in webcam_bursts(int(a.source), n))
    else:
        results = estimate_stream(file_frames(a.source), L=a.L, roi=a.roi, every=a.every)

    if a.show:
        import matplotlib.pyplot as plt
        plt.ion()
        fig, ax = plt.subplots(figsize=(7, 3.5))
        ax.set_xlabel("time (s)")
        ax.set_ylabel("speed (px/frame)")
        good, = ax.plot([], [], "o", color="#2a78d6", label="valid")
        bad, = ax.plot([], [], "o", mfc="none", color="#eb6834", label="not valid")
        ax.legend(loc="upper left", frameon=False)
        hist = []

    print(f"{'t (s)':>8} {'dir (deg)':>9} {'px/frame':>9} {'px/s':>8} {'quality':>7}  state")
    for r in results:
        state = "valid" if r.valid else ("uneven frame steps" if not r.regular else "not valid")
        print(f"{r.t:8.2f} {r.angle_deg:9.1f} {r.speed:9.2f} {r.speed_per_s:8.1f} {r.quality:7.2f}  {state}")
        if a.show:
            hist.append(r)
            v = [h for h in hist if h.valid]
            nv = [h for h in hist if not h.valid]
            good.set_data([h.t for h in v], [h.speed for h in v])
            bad.set_data([h.t for h in nv], [h.speed for h in nv])
            ax.relim()
            ax.autoscale_view()
            plt.pause(0.001)


if __name__ == "__main__":
    main()
