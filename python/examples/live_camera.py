"""Velocity from a live camera or a video file.

    python python/examples/live_camera.py --source pylon            # Basler camera (pip install pypylon)
    python python/examples/live_camera.py --source 0                # any webcam through OpenCV
    python python/examples/live_camera.py --source recording.mp4    # video file

Every estimate uses 2L+1 consecutive frames. The pulsed-exposure model needs
equal time steps between them; windows with uneven frame intervals are
reported as not valid. ``--show`` opens a live plot.
"""

from __future__ import annotations

import argparse
import time
from typing import Iterator, Tuple

import numpy as np

from fftvel.stream import estimate_stream


def pylon_frames(frame_rate: float | None) -> Iterator[Tuple[np.ndarray, float]]:
    """Frames and camera timestamps from the first Basler camera (pypylon)."""
    from pypylon import pylon

    cam = pylon.InstantCamera(pylon.TlFactory.GetInstance().CreateFirstDevice())
    cam.Open()
    try:
        try:
            cam.PixelFormat.SetValue("Mono8")
        except Exception:  # noqa: BLE001 - not every model exposes it
            pass
        if frame_rate:
            cam.AcquisitionFrameRateEnable.SetValue(True)
            cam.AcquisitionFrameRate.SetValue(frame_rate)
        # USB3 cameras count timestamp ticks in ns; GigE models may use another
        # tick rate (see GevTimestampTickFrequency).
        cam.StartGrabbing(pylon.GrabStrategy_OneByOne)
        while cam.IsGrabbing():
            res = cam.RetrieveResult(5000, pylon.TimeoutHandling_ThrowException)
            try:
                if res.GrabSucceeded():
                    yield res.Array.copy(), res.TimeStamp * 1e-9
            finally:
                res.Release()
    finally:
        cam.StopGrabbing()
        cam.Close()


def opencv_frames(source) -> Iterator[Tuple[np.ndarray, float]]:
    """Frames from OpenCV; file timestamps from the container, webcam ones from the clock."""
    import cv2

    cap = cv2.VideoCapture(source)
    if not cap.isOpened():
        raise SystemExit(f"cannot open {source!r}")
    is_file = isinstance(source, str)
    try:
        while True:
            ok, frame = cap.read()
            if not ok:
                break
            t = cap.get(cv2.CAP_PROP_POS_MSEC) * 1e-3 if is_file else time.monotonic()
            yield frame[..., ::-1] if frame.ndim == 3 else frame, t    # BGR -> RGB
    finally:
        cap.release()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--source", default="pylon", help="'pylon', a webcam index or a video file")
    p.add_argument("--L", type=int, default=4, help="half number of frames per estimate")
    p.add_argument("--roi", type=int, default=512, help="side of the analysed central square (px)")
    p.add_argument("--every", type=int, help="new frames between estimates (default 2L+1)")
    p.add_argument("--frame-rate", type=float, help="set the Basler frame rate (Hz)")
    p.add_argument("--show", action="store_true", help="live plot of speed over time")
    a = p.parse_args()

    if a.source == "pylon":
        frames = pylon_frames(a.frame_rate)
    elif a.source.isdigit():
        frames = opencv_frames(int(a.source))
    else:
        frames = opencv_frames(a.source)

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
    for r in estimate_stream(frames, L=a.L, roi=a.roi, every=a.every):
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
