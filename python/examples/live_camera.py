"""Velocity from a live camera or a video file.

    python python/examples/live_camera.py --mode flash --source pylon --L 4 --flash-period 0.001
    python python/examples/live_camera.py --mode average --source pylon --exposure-us 500
    python python/examples/live_camera.py --mode average --source 0           # any webcam through OpenCV
    python python/examples/live_camera.py --mode flash --source recording.avi

Basler cameras need pypylon (pip install pypylon), webcams and video files
OpenCV (pip install -e "./python[examples]", or keep an opencv-python that is
already installed and use "./python[plot,data]"; see the README).

--mode must be given:

flash    Every frame is one pulsed image: a strobe fires K = 2L+1 times (or
         --exposures K) during one exposure. Each frame is estimated on its
         own; live cameras deliver the newest frame for every estimate. The
         speed is in px per flash step, and in px/s with --flash-period.
         --duty is the flash pulse length divided by the flash period. The
         flash timing cannot be checked from the frame times.
average  Continuous light: K consecutive frames are averaged into a pulsed
         image (do not use this for multi-flash frames). Live cameras
         deliver a fresh burst of K frames per estimate, so frames buffered
         while estimating never end up in a window. The model needs equal
         frame steps: a window whose frame intervals spread (max - min) by
         more than --tol (0.1) times their median is reported as not valid
         ("uneven frame steps"). The speed is in px per frame and px/s. For
         Basler cameras the exposure time is read from the camera and passed
         to the estimator as duty = exposure / frame step; keep it below
         about 0.3 of the step.

Frame times: Basler cameras use the camera clock (the host clock if the
camera sends no timestamps, as the pylon emulator; the window is then checked
against the nominal frame period), webcams the driver's buffer time or the
host clock. Video files only carry nominal times (frame number / frame rate
stored in the file): they cannot reveal dropped frames, and px/s is only
right if that frame rate is. Image sequences (img_%04d.png) need --fps.

--show opens a live plot. Ctrl+C stops.
"""

from __future__ import annotations

import argparse
import time
from typing import Iterator, List, Optional, Tuple

import numpy as np

from fftvel.stream import estimate_frame, estimate_stream, estimate_window, frame_step

Frames = List[Tuple[np.ndarray, float]]
MAX_FAILED = 20       # consecutive failed grabs, or bursts with a failed grab, before giving up
READ_RETRIES = 10     # webcam: attempts to read a frame before giving up
MONO = ("Mono8", "Mono10", "Mono12", "Mono16")    # unpacked formats that pypylon delivers as arrays


def _cv2():
    """The OpenCV module, or a message how to install it."""
    try:
        import cv2
    except ImportError:
        raise SystemExit('webcams and video files need OpenCV: pip install opencv-python '
                         '(or pip install -e "./python[video]")') from None
    return cv2


def _try_set(cam, name: str, value) -> bool:
    try:
        getattr(cam, name).SetValue(value)
        return True
    except Exception:  # noqa: BLE001 - node missing or not writable on this camera model
        return False


def _get(cam, name: str):
    try:
        return getattr(cam, name).GetValue()
    except Exception:  # noqa: BLE001
        return None


def _tick_seconds(cam) -> float:
    """Length of one timestamp tick: GigE models report their tick frequency, USB3 counts ns."""
    try:
        return 1.0 / float(cam.GevTimestampTickFrequency.GetValue())
    except Exception:  # noqa: BLE001
        return 1e-9


def _exposure_seconds(cam) -> float:
    us = _get(cam, "ExposureTime") or _get(cam, "ExposureTimeAbs")    # USB3 and newer GigE, GigE ace classic
    return float(us) * 1e-6 if us else float("nan")


def _mono_reader(pylon, fmt: Optional[str]):
    """Function that returns a grab result as a 2-D array, for camera pixel format ``fmt``."""
    if fmt in MONO:
        return lambda res: res.Array.copy()           # Mono10/12/16 arrive as uint16, e.g. 0..4095
    if fmt and fmt.startswith("Mono"):                 # packed (Mono12p, ...): convert to Mono16
        conv = pylon.ImageFormatConverter()
        conv.OutputPixelFormat = pylon.PixelType_Mono16
        return lambda res: conv.Convert(res).Array.copy()
    raise SystemExit(f"The camera pixel format is {fmt} and Mono8 could not be set. Colour and Bayer "
                     f"formats are not supported, because the raw colour pattern would be analysed as "
                     f"texture. Select a Mono format (e.g. in the pylon Viewer) and save it in the user "
                     f"set the camera starts with.")


def _open_pylon(frame_rate: Optional[float], exposure_us: Optional[float]):
    """Open and configure the first Basler camera; return ``(pylon, cam, read, period)``.

    ``read`` turns a grab result into a 2-D array, ``period`` is the nominal
    frame period in s that the camera reports (None if unknown).
    """
    try:
        from pypylon import pylon
    except ImportError:
        raise SystemExit("--source pylon needs pypylon: pip install pypylon") from None

    cam = pylon.InstantCamera(pylon.TlFactory.GetInstance().CreateFirstDevice())
    cam.Open()
    try:
        _try_set(cam, "PixelFormat", "Mono8")
        fmt = _get(cam, "PixelFormat")
        read = _mono_reader(pylon, fmt)
        if frame_rate:
            ok = _try_set(cam, "AcquisitionFrameRateEnable", True)
            if not (_try_set(cam, "AcquisitionFrameRate", frame_rate)          # USB3, newer GigE
                    or _try_set(cam, "AcquisitionFrameRateAbs", frame_rate)):  # GigE ace classic
                print("warning: could not set the frame rate" + ("" if ok else " (no frame-rate node)"))
        if exposure_us:
            _try_set(cam, "ExposureAuto", "Off")
            if not (_try_set(cam, "ExposureTime", exposure_us) or _try_set(cam, "ExposureTimeAbs", exposure_us)):
                print("warning: could not set the exposure time")
        elif _get(cam, "ExposureAuto") not in (None, "Off"):
            print("note: auto exposure is on, so the exposure time can change; set --exposure-us")
        rate = _get(cam, "ResultingFrameRate") or _get(cam, "ResultingFrameRateAbs")
        print(f"camera: pixel format {fmt}, exposure {_exposure_seconds(cam) * 1e6:.0f} us, "
              f"frame rate {rate or float('nan'):.2f} Hz")
    except BaseException:
        cam.Close()
        raise
    return pylon, cam, read, (1.0 / rate if rate else None)


def _retrieve(pylon, cam, failed: int):
    """Next grab result, or None after printing why the grab failed. Gives up after MAX_FAILED in a row."""
    try:
        res = cam.RetrieveResult(5000, pylon.TimeoutHandling_ThrowException)
    except pylon.TimeoutException as err:
        raise SystemExit(f"{err}\nIf the camera waits for a hardware trigger (TriggerMode On), "
                         f"check the trigger signal.") from None
    if res.GrabSucceeded():
        return res
    print(f"grab failed ({failed + 1} in a row): {res.ErrorCode:#x} {res.ErrorDescription}")
    res.Release()
    if failed + 1 >= MAX_FAILED:
        raise SystemExit(f"{MAX_FAILED} grabs in a row failed. On GigE enable jumbo frames or lower "
                         f"GevSCPSPacketSize; otherwise lower the frame rate, the ROI or "
                         f"DeviceLinkThroughputLimit.")
    return None


def pylon_bursts(n: int, frame_rate: Optional[float] = None,
                 exposure_us: Optional[float] = None) -> Iterator[Tuple[Frames, float, Optional[float]]]:
    """Bursts of ``n`` consecutive frames from the first Basler camera (average mode).

    Yields ``(burst, exposure, period)``: the ``(frame, t)`` pairs with ``t``
    in s since the first frame, the exposure time in s (nan if unknown) and
    ``period``. If the camera delivers usable timestamps (checked on the
    first burst), ``t`` is camera time and ``period`` is None. Otherwise
    ``t`` is the host clock at retrieval and ``period`` is the nominal frame
    period, to check the window against (None if unknown).
    """
    pylon, cam, read, nominal = _open_pylon(frame_rate, exposure_us)
    try:
        tick = _tick_seconds(cam)
        cam.MaxNumBuffer.SetValue(max(10, n + 2))
        failed = bad = 0
        host_clock = None
        t0 = None
        while True:
            cam.StartGrabbingMax(n, pylon.GrabStrategy_OneByOne)
            frames, cam_t, host_t = [], [], []
            while cam.IsGrabbing():
                res = _retrieve(pylon, cam, failed)
                if res is None:
                    failed += 1
                    continue
                failed = 0
                try:
                    host_t.append(time.perf_counter())
                    cam_t.append(res.TimeStamp * tick)
                    frames.append(read(res))
                finally:
                    res.Release()
            if len(frames) < n:
                bad += 1
                if bad >= MAX_FAILED:
                    raise SystemExit(f"{MAX_FAILED} bursts in a row lost frames, so no window of {n} "
                                     f"consecutive frames was complete. Lower the frame rate or the ROI, "
                                     f"or on GigE enable jumbo frames.")
                continue
            bad = 0
            if host_clock is None:
                host_clock = not np.all(np.diff(cam_t) > 0)
                if host_clock:
                    print("note: the camera delivers no usable timestamps (all equal or not increasing); "
                          "using the host clock" + (" and the nominal frame period" if nominal else ""))
            times = np.asarray(host_t if host_clock else cam_t)
            if t0 is None:
                t0 = times[0]
            yield list(zip(frames, times - t0)), _exposure_seconds(cam), (nominal if host_clock else None)
    finally:
        cam.StopGrabbing()
        cam.Close()


def pylon_frames(frame_rate: Optional[float] = None,
                 exposure_us: Optional[float] = None) -> Iterator[Tuple[np.ndarray, float, float]]:
    """The newest frame from the first Basler camera, each time one is requested (flash mode).

    Yields ``(frame, t, exposure)`` with ``t`` in s since the first frame
    (camera clock, or host clock if the camera sends no timestamps) and the
    exposure time in s (nan if unknown).
    """
    pylon, cam, read, _ = _open_pylon(frame_rate, exposure_us)
    try:
        tick = _tick_seconds(cam)
        cam.StartGrabbing(pylon.GrabStrategy_LatestImageOnly)
        failed = 0
        t0 = None
        host_clock = False
        while True:
            res = _retrieve(pylon, cam, failed)
            if res is None:
                failed += 1
                continue
            failed = 0
            try:
                t = res.TimeStamp * tick
                frame = read(res)
            finally:
                res.Release()
            if t0 is None:
                host_clock = not t > 0
                if host_clock:
                    print("note: the camera delivers no usable timestamps; using the host clock")
            if host_clock:
                t = time.perf_counter()
            if t0 is None:
                t0 = t
            yield frame, t - t0, _exposure_seconds(cam)
    finally:
        cam.StopGrabbing()
        cam.Close()


def _open_webcam(index: int):
    cv2 = _cv2()
    cap = cv2.VideoCapture(index)
    if not cap.isOpened():
        raise SystemExit(f"cannot open camera {index}")
    cap.set(cv2.CAP_PROP_BUFFERSIZE, 1)            # not every backend honours this
    return cv2, cap


def _read(cap, index: int) -> np.ndarray:
    """Next webcam frame (RGB or grey); retries a failed read before giving up."""
    for k in range(READ_RETRIES):
        ok, frame = cap.read()
        if ok and frame is not None:
            return frame[..., ::-1] if frame.ndim == 3 else frame     # BGR -> RGB
        print(f"camera {index}: no frame, retrying ({k + 1}/{READ_RETRIES})")
        time.sleep(0.2)
    raise SystemExit(f"camera {index}: no frame after {READ_RETRIES} attempts (OpenCV backend "
                     f"{cap.getBackendName()}). Check the connection, close other programs that use "
                     f"the camera, or try another camera index.")


def webcam_bursts(index: int, n: int) -> Iterator[Frames]:
    """Bursts of ``n`` frames from a webcam (average mode); frames buffered during an estimate are dropped first.

    The frame intervals come from the driver's buffer times if they are
    positive and increasing, else from the host clock. Times are in s since
    the first frame, on the host clock, so they stay comparable when the
    driver times are missing in some bursts only.
    """
    cv2, cap = _open_webcam(index)
    try:
        t0 = None
        while True:
            for _ in range(5):                       # drain stale buffers
                cap.grab()
            frames, driver_t, host_t = [], [], []
            for _ in range(n):
                frames.append(_read(cap, index))
                host_t.append(time.perf_counter())
                driver_t.append(cap.get(cv2.CAP_PROP_POS_MSEC) * 1e-3)
            driver_t = np.asarray(driver_t)
            if driver_t.min() > 0 and np.all(np.diff(driver_t) > 0):
                times = host_t[0] + (driver_t - driver_t[0])
            else:
                times = np.asarray(host_t)
            if t0 is None:
                t0 = times[0]
            yield list(zip(frames, times - t0))
    finally:
        cap.release()


def webcam_frames(index: int) -> Iterator[Tuple[np.ndarray, float]]:
    """The newest webcam frame each time one is requested (flash mode), with host time since the first."""
    cv2, cap = _open_webcam(index)
    try:
        t0 = time.perf_counter()
        while True:
            for _ in range(5):                       # drain stale buffers
                cap.grab()
            frame = _read(cap, index)
            yield frame, time.perf_counter() - t0
    finally:
        cap.release()


def file_frames(path: str, fps: Optional[float] = None) -> Iterator[Tuple[np.ndarray, float]]:
    """Every frame of a video file with its nominal time (container time, or frame number / ``fps``)."""
    cv2 = _cv2()
    cap = cv2.VideoCapture(path)
    if not cap.isOpened():
        raise SystemExit(f"cannot open {path!r}")
    try:
        k = 0
        while True:
            ok, frame = cap.read()
            if not ok:
                break
            t = k / fps if fps else cap.get(cv2.CAP_PROP_POS_MSEC) * 1e-3
            yield (frame[..., ::-1] if frame.ndim == 3 else frame), t
            k += 1
    finally:
        cap.release()


def exposure_duty(exposure: float, step: float, warned: set) -> float:
    """Duty for the estimator in average mode: exposure time / frame step, with a warning for risky values."""
    if not (exposure > 0 and step > 0):
        return 0.0
    duty = exposure / step
    band = "half" if 0.4 <= duty <= 0.6 else ("full" if duty > 0.9 else None)
    if band and band not in warned:
        warned.add(band)
        if band == "half":
            print(f"warning: the exposure is {duty:.2f} of the frame step. Between 0.4 and 0.6 the "
                  f"speed fit depends strongly on the exact exposure (at 0.5 every second spectral "
                  f"line vanishes); use an exposure below about 0.3 of the step.")
        else:
            print(f"warning: the exposure is {duty:.2f} of the frame step, so the frames are strongly "
                  f"blurred and most estimates will not be valid; shorten the exposure.")
    return min(duty, 0.99)          # the estimator needs duty < 1


def _average_pylon(a, K: int, est: dict):
    warned: set = set()
    for burst, exposure, period in pylon_bursts(K, a.frame_rate, a.exposure_us):
        duty = a.duty
        if duty is None:
            step, _ = frame_step([t for _, t in burst], a.tol, period)
            duty = exposure_duty(exposure, step, warned)
        yield estimate_window(burst, tol=a.tol, duty=duty, period=period, **est)


def _flash_pylon(a, K: int, est: dict):
    checked = False
    for frame, t, exposure in pylon_frames(a.frame_rate, a.exposure_us):
        if not checked and a.flash_period and exposure > 0:
            checked = True
            need = (K - 1 + (a.duty or 0)) * a.flash_period
            if exposure < need:
                print(f"warning: the exposure ({exposure * 1e6:.0f} us) is shorter than {K} flashes "
                      f"{a.flash_period * 1e6:.0f} us apart need ({need * 1e6:.0f} us)")
        yield estimate_frame(frame, t, flash_period=a.flash_period, duty=a.duty or 0.0, **est)


def _state(r) -> str:
    if r.valid:
        return "valid"
    if not r.regular:
        return "uneven frame steps" if r.step > 0 else "no frame times"
    e = r.estimate
    if e is not None and e.isotropic:
        return "not valid (isotropic: no motion, or defocus?)"
    if e is not None and not e.in_range:
        return "not valid (below the measurable range)"
    return "not valid (low quality)"


def _parse(argv):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mode", required=True, choices=["flash", "average"],
                   help="flash: every frame holds K strobe flashes; average: K frames are averaged")
    p.add_argument("--source", default="pylon", help="'pylon', a webcam index or a video file")
    p.add_argument("--L", type=int, default=4, help="half number of flashes or frames, K = 2L+1 (4)")
    p.add_argument("--exposures", type=int, help="exact number K >= 2 of flashes or frames, overrides --L")
    p.add_argument("--duty", type=float, help="flash mode: flash pulse / flash period (0); average mode: "
                   "exposure / frame step (default: read from a Basler camera, else 0)")
    p.add_argument("--flash-period", type=float, help="flash mode: time between flashes in s, for px/s")
    p.add_argument("--roi", type=int, default=512,
                   help="side of the analysed central square in px (512 takes several seconds per estimate)")
    p.add_argument("--tol", type=float, default=0.1,
                   help="average mode: allowed spread of the frame intervals, relative to their median (0.1)")
    p.add_argument("--clip-outliers", action="store_true", help="clip static bright spots before the FFT")
    p.add_argument("--every", type=int, help="video files: frames between estimates (default: K in average "
                   "mode, 1 in flash mode)")
    p.add_argument("--fps", type=float, help="video files: frame rate; replaces the times stored in the file "
                   "(needed for image sequences)")
    p.add_argument("--frame-rate", type=float, help="Basler: set the frame rate (Hz)")
    p.add_argument("--exposure-us", type=float, help="Basler: set the exposure time (us), auto exposure off")
    p.add_argument("--show", action="store_true", help="live plot of speed over time")
    a = p.parse_args(argv)
    if a.L < 1:
        p.error("--L must be at least 1")
    if a.exposures is not None and a.exposures < 2:
        p.error("--exposures must be at least 2")
    if a.duty is not None and not 0 <= a.duty < 1:
        p.error("--duty must be in [0, 1)")
    if a.roi < 16:
        p.error("--roi must be at least 16")
    if not 0 <= a.tol <= 1:
        p.error("--tol must be in [0, 1]")
    if a.every is not None and a.every < 1:
        p.error("--every must be at least 1")
    for name in ("flash_period", "fps", "frame_rate", "exposure_us"):
        v = getattr(a, name)
        if v is not None and not (np.isfinite(v) and v > 0):
            p.error(f"--{name.replace('_', '-')} must be a positive number")
    if a.mode == "average" and a.flash_period is not None:
        p.error("--flash-period is for --mode flash; in average mode the frame times give px/s")
    if a.source != "pylon" and (a.frame_rate or a.exposure_us):
        p.error("--frame-rate and --exposure-us need --source pylon")
    return a


def main(argv=None):
    a = _parse(argv)
    K = a.exposures or 2 * a.L + 1
    est = dict(L=a.L, roi=a.roi, exposures=a.exposures, clip_outliers=a.clip_outliers)
    is_file = a.source != "pylon" and not a.source.isdigit()
    if is_file and "%" in a.source and not a.fps:
        raise SystemExit("an image sequence has no frame times (OpenCV assumes 25 fps); give --fps")

    if a.mode == "flash":
        duty = a.duty or 0.0
        if a.source == "pylon":
            results = _flash_pylon(a, K, est)
        elif a.source.isdigit():
            results = (estimate_frame(f, t, flash_period=a.flash_period, duty=duty, **est)
                       for f, t in webcam_frames(int(a.source)))
        else:
            every = a.every or 1
            results = (estimate_frame(f, t, flash_period=a.flash_period, duty=duty, **est)
                       for k, (f, t) in enumerate(file_frames(a.source, a.fps)) if k % every == 0)
        unit = "px/flash"
    else:
        if a.source == "pylon":
            results = _average_pylon(a, K, est)
        elif a.source.isdigit():
            results = (estimate_window(b, tol=a.tol, duty=a.duty or 0.0, **est)
                       for b in webcam_bursts(int(a.source), K))
        else:
            results = estimate_stream(file_frames(a.source, a.fps), every=a.every, tol=a.tol,
                                      duty=a.duty or 0.0, **est)
        unit = "px/frame"
    if is_file:
        print("note: the times of a video file are nominal (frame number / frame rate"
              + (" given by --fps)" if a.fps else " stored in the file)")
              + "; dropped frames cannot be detected from them")
    if a.mode == "flash" and not a.flash_period:
        print("note: px/s needs --flash-period")

    if a.show:
        import matplotlib.pyplot as plt
        plt.ion()
        fig, ax = plt.subplots(figsize=(7, 3.5))
        ax.set_xlabel("nominal time (s)" if is_file else "time (s)")
        ax.set_ylabel(f"speed ({unit})")
        good, = ax.plot([], [], "o", color="#2a78d6", label="valid")
        bad, = ax.plot([], [], "o", mfc="none", color="#eb6834", label="not valid")
        ax.legend(loc="upper left", frameon=False)
        hist = []

    tlabel = "t nom (s)" if is_file else "t (s)"
    header = f"{tlabel:>9} {'dir (deg)':>9} {unit:>9} {'px/s':>8} {'quality':>7}  state"
    try:
        for r in results:
            if header:                  # printed after the camera has reported its settings
                print(header)
                header = ""
            print(f"{r.t:9.3f} {r.angle_deg:9.1f} {r.speed:9.2f} {r.speed_per_s:8.1f} {r.quality:7.2f}  {_state(r)}")
            if a.show:
                if not plt.fignum_exists(fig.number):
                    break
                hist.append(r)
                v = [h for h in hist if h.valid]
                nv = [h for h in hist if not h.valid]
                good.set_data([h.t for h in v], [h.speed for h in v])
                bad.set_data([h.t for h in nv], [h.speed for h in nv])
                ax.relim()
                ax.autoscale_view()
                plt.pause(0.001)
    except KeyboardInterrupt:
        pass
    finally:
        if hasattr(results, "close"):
            results.close()
    if header:
        print("no estimate" + (f" (average mode needs {K} frames)" if is_file and a.mode == "average" else ""))
    if a.show:
        plt.ioff()
        plt.show()          # keeps the plot open after a file ends or after Ctrl+C
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
