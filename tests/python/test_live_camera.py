"""python/examples/live_camera.py with fake OpenCV and pypylon modules (no camera needed)."""

import importlib.util
import sys
import types

import numpy as np
import pytest

from fftvel import simulate_pulsed_image
from fftvel.stream import StreamEstimate

from conftest import REPO


@pytest.fixture
def lc():
    spec = importlib.util.spec_from_file_location("live_camera", REPO / "python" / "examples" / "live_camera.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class _Clock:
    """Stands in for the time module: perf_counter advances by ``step`` per call, sleep does nothing."""

    def __init__(self, step):
        self.t = 100.0
        self.step = step

    def perf_counter(self):
        self.t += self.step
        return self.t

    def sleep(self, s):
        pass


# --- fake OpenCV ---------------------------------------------------------------

def _fake_cv2(frames, msec=lambda k: 40.0 * k, fail=lambda k: False):
    """cv2 stand-in whose VideoCapture returns ``frames`` (cycling) and POS_MSEC ``msec(k)``."""

    class VideoCapture:
        def __init__(self, source):
            self.k = -1
            self.reads = 0

        def isOpened(self):
            return True

        def set(self, *args):
            return False

        def grab(self):
            return True

        def read(self):
            self.reads += 1
            if fail(self.reads):
                return False, None
            self.k += 1
            if isinstance(frames, list) and self.k >= len(frames):
                return False, None
            f = frames[self.k % len(frames)]
            return True, np.repeat(f[:, :, None], 3, axis=2)       # BGR

        def get(self, prop):
            return msec(self.k)

        def getBackendName(self):
            return "FAKE"

        def release(self):
            pass

    return types.SimpleNamespace(VideoCapture=VideoCapture, CAP_PROP_BUFFERSIZE=38, CAP_PROP_POS_MSEC=0)


def _uint8(img):
    return np.clip(128 + 40 * img / img.std(), 0, 255).astype(np.uint8)


def _table(out):
    """Rows of the result table as (t, angle, speed, px/s, quality, state)."""
    rows = []
    for line in out.splitlines():
        parts = line.split(None, 5)
        try:
            rows.append(tuple(float(x) for x in parts[:5]) + (parts[5],))
        except (ValueError, IndexError):
            pass
    return rows


def test_flash_mode_on_a_video_file(lc, monkeypatch, capsys):
    flashes = [_uint8(simulate_pulsed_image(128, (3.0, -2.0), 4, rng=k)) for k in range(3)]
    monkeypatch.setitem(sys.modules, "cv2", _fake_cv2(flashes))
    assert lc.main(["--mode", "flash", "--source", "rec.avi", "--roi", "128", "--flash-period", "0.001"]) == 0
    out = capsys.readouterr().out
    assert "nominal" in out and "t nom (s)" in out and "px/flash" in out
    rows = _table(out)
    assert len(rows) == 3
    for t, angle, speed, pxs, q, state in rows:
        assert state == "valid"
        assert abs(speed / np.hypot(3, 2) - 1) < 0.01 and abs(angle - 33.69) < 0.5
        assert abs(pxs / (1000 * np.hypot(3, 2)) - 1) < 0.01
    assert [r[0] for r in rows] == [0.0, 0.04, 0.08]


def test_average_mode_on_a_video_file_and_fps(lc, monkeypatch, capsys):
    _, fr = simulate_pulsed_image(128, (2.0, 1.0), 4, rng=3, return_frames=True)
    frames = [_uint8(fr[:, :, k]) for k in range(9)]
    monkeypatch.setitem(sys.modules, "cv2", _fake_cv2(frames))
    assert lc.main(["--mode", "average", "--source", "rec.avi", "--roi", "128", "--fps", "200"]) == 0
    rows = _table(capsys.readouterr().out)
    assert len(rows) == 1
    t, angle, speed, pxs, q, state = rows[0]
    assert state == "valid" and abs(t - 8 / 200) < 1e-9
    assert abs(speed / np.hypot(2, 1) - 1) < 0.01 and abs(pxs / (200 * np.hypot(2, 1)) - 1) < 0.01
    with pytest.raises(SystemExit, match="--fps"):
        lc.main(["--mode", "average", "--source", "img_%04d.png"])     # image sequence without --fps


def test_options_are_checked(lc):
    for argv in (["--source", "0"],                                            # no --mode
                 ["--mode", "average", "--flash-period", "0.001"],
                 ["--mode", "flash", "--source", "rec.avi", "--exposure-us", "100"],
                 ["--mode", "flash", "--duty", "1"],
                 ["--mode", "flash", "--exposures", "1"],
                 ["--mode", "flash", "--flash-period", "0"]):
        with pytest.raises(SystemExit):
            lc.main(argv)


def test_webcam_retries_and_host_clock(lc, monkeypatch, capsys):
    _, fr = simulate_pulsed_image(96, (2.0, 1.0), 4, rng=4, return_frames=True)
    frames = np.stack([_uint8(fr[:, :, k]) for k in range(9)])
    monkeypatch.setattr(lc, "time", _Clock(1 / 30))
    # driver time -1 (as DirectShow): the host clock is used; the first two reads fail
    monkeypatch.setitem(sys.modules, "cv2", _fake_cv2(frames, msec=lambda k: -1.0, fail=lambda n: n <= 2))
    burst = next(lc.webcam_bursts(0, 9))
    assert "retrying (2/10)" in capsys.readouterr().out
    assert np.allclose(np.diff([t for _, t in burst]), 1 / 30) and burst[0][1] == 0
    r = lc.estimate_window(burst, L=4)
    assert r.regular and abs(r.speed_per_s / (30 * np.hypot(2, 1)) - 1) < 0.01
    monkeypatch.setitem(sys.modules, "cv2", _fake_cv2(frames, fail=lambda n: True))
    with pytest.raises(SystemExit, match="no frame after 10 attempts"):
        next(lc.webcam_frames(0))


def test_webcam_driver_times_missing_in_the_first_burst(lc, monkeypatch):
    # driver time 0 during the first burst (as MSMF at start), then 25 ms steps from the driver;
    # the times stay on one time base, and the intervals of later bursts come from the driver
    frames = np.zeros((9, 32, 32), dtype=np.uint8)
    clock = _Clock(0.03)
    monkeypatch.setattr(lc, "time", clock)
    monkeypatch.setitem(sys.modules, "cv2", _fake_cv2(frames, msec=lambda k: 0.0 if k < 9 else 25.0 * k))
    gen = lc.webcam_bursts(0, 9)
    t1 = [t for _, t in next(gen)]
    t2 = [t for _, t in next(gen)]
    gen.close()
    assert t1[0] == 0 and np.allclose(np.diff(t1), 0.03)                 # host clock
    assert np.allclose(np.diff(t2), 0.025) and t2[0] > t1[-1]           # driver intervals, later times


def test_state_column(lc):
    def r(valid=False, regular=True, step=0.01, isotropic=False, in_range=True):
        e = types.SimpleNamespace(isotropic=isotropic, in_range=in_range)
        return StreamEstimate(0.0, 0.0, 1.0, 100.0, 0.5, valid, regular, step, e)
    assert lc._state(r(valid=True)) == "valid"
    assert lc._state(r(regular=False)) == "uneven frame steps"
    assert lc._state(r(regular=False, step=0.0)) == "no frame times"
    assert lc._state(r(isotropic=True)).startswith("not valid (isotropic")
    assert lc._state(r(in_range=False)) == "not valid (below the measurable range)"
    assert lc._state(r()) == "not valid (low quality)"


def test_short_video_file_gives_a_message(lc, monkeypatch, capsys):
    frames = [np.zeros((32, 32), dtype=np.uint8)] * 5
    monkeypatch.setitem(sys.modules, "cv2", _fake_cv2(frames))
    assert lc.main(["--mode", "average", "--source", "rec.avi", "--roi", "32"]) == 0
    assert "no estimate (average mode needs 9 frames)" in capsys.readouterr().out


# --- fake pypylon ------------------------------------------------------------

class _Node:
    def __init__(self, value, writable=True):
        self.value = value
        self.writable = writable

    def SetValue(self, v):
        if not self.writable:
            raise RuntimeError("node not writable")
        self.value = v

    def GetValue(self):
        return self.value


class _Result:
    def __init__(self, ok, array, ts):
        self.ok = ok
        self.Array = array
        self.TimeStamp = ts
        self.ErrorCode = 0 if ok else 0xE1000014
        self.ErrorDescription = "" if ok else "The buffer was incompletely grabbed."

    def GrabSucceeded(self):
        return self.ok

    def Release(self):
        pass


class _Camera:
    """Basler camera stand-in. ``nodes`` holds the GenICam features it has; others are missing."""

    def __init__(self, frames, ts=lambda k: 0, ok=lambda k: True, **nodes):
        self.frames, self.ts, self.ok = frames, ts, ok
        self.nodes = dict(PixelFormat=_Node("Mono8"), MaxNumBuffer=_Node(10), ExposureAuto=_Node("Off"),
                          ExposureTime=_Node(1000.0), AcquisitionFrameRateEnable=_Node(False),
                          AcquisitionFrameRate=_Node(100.0), ResultingFrameRate=_Node(100.0))
        self.nodes.update({k: v for k, v in nodes.items() if v is not None})
        for k in [k for k, v in nodes.items() if v is None]:
            del self.nodes[k]
        self.k = 0
        self.left = 0
        self.open = False

    def __getattr__(self, name):
        if name in self.__dict__.get("nodes", {}):
            return self.nodes[name]
        raise AttributeError(name)

    def Open(self):
        self.open = True

    def Close(self):
        self.open = False

    def StartGrabbingMax(self, n, strategy):
        self.left = n

    def StartGrabbing(self, strategy):
        self.left = np.inf

    def IsGrabbing(self):
        return self.left > 0

    def StopGrabbing(self):
        self.left = 0

    def RetrieveResult(self, timeout, mode):
        self.left -= 1
        k = self.k
        self.k += 1
        return _Result(self.ok(k), self.frames[k % len(self.frames)], self.ts(k))


def _fake_pylon(monkeypatch, cam):
    class Converter:
        OutputPixelFormat = None

        def Convert(self, res):
            return types.SimpleNamespace(Array=res.Array.astype(np.uint16) * 16)

    pylon = types.SimpleNamespace(
        InstantCamera=lambda dev: cam, TlFactory=types.SimpleNamespace(GetInstance=lambda: types.SimpleNamespace(
            CreateFirstDevice=lambda: None)),
        GrabStrategy_OneByOne=0, GrabStrategy_LatestImageOnly=1, TimeoutHandling_ThrowException=0,
        TimeoutException=TimeoutError, ImageFormatConverter=Converter, PixelType_Mono16=16)
    monkeypatch.setitem(sys.modules, "pypylon", types.SimpleNamespace(pylon=pylon))


def test_pylon_without_timestamps_uses_host_clock_and_exposure(lc, monkeypatch, capsys):
    _, fr = simulate_pulsed_image(96, (2.0, 1.0), 4, rng=5, return_frames=True)
    frames = [_uint8(fr[:, :, k]) for k in range(9)]
    # GigE ace classic: only the ...Abs nodes exist; the camera sends timestamp 0 (as the pylon emulator)
    cam = _Camera(frames, ExposureTime=None, AcquisitionFrameRate=None, ResultingFrameRate=None,
                  ExposureTimeAbs=_Node(4000.0), AcquisitionFrameRateAbs=_Node(100.0),
                  ResultingFrameRateAbs=_Node(100.0))
    _fake_pylon(monkeypatch, cam)
    monkeypatch.setattr(lc, "time", types.SimpleNamespace(
        perf_counter=iter(0.01 * np.arange(100) + np.tile([0, 0.002], 50)).__next__))   # 2 ms host jitter
    gen = lc.pylon_bursts(9, frame_rate=50, exposure_us=2000)
    burst, exposure, period = next(gen)
    gen.close()
    out = capsys.readouterr().out
    assert "no usable timestamps" in out and not cam.open
    assert cam.nodes["AcquisitionFrameRateAbs"].value == 50 and cam.nodes["ExposureTimeAbs"].value == 2000
    assert exposure == pytest.approx(0.002) and period == pytest.approx(0.01)
    r = lc.estimate_window(burst, L=4, period=period)
    assert r.regular and r.valid and r.step == pytest.approx(0.01)
    assert lc.estimate_window(burst, L=4).regular is False          # 40 % spread without the nominal period
    warned = set()
    assert lc.exposure_duty(0.001, 0.01, warned) == pytest.approx(0.1) and not warned
    assert lc.exposure_duty(0.005, 0.01, warned) == pytest.approx(0.5) and warned == {"half"}
    assert lc.exposure_duty(0.02, 0.01, warned) == 0.99 and warned == {"half", "full"}
    assert lc.exposure_duty(float("nan"), 0.01, warned) == 0.0
    assert "0.50 of the frame step" in capsys.readouterr().out


def test_pylon_camera_timestamps_and_flash_frames(lc, monkeypatch):
    frames = [np.full((40, 40), k, dtype=np.uint8) for k in range(20)]
    cam = _Camera(frames, ts=lambda k: 5_000_000_000 + 10_000_000 * k)       # ns, 100 fps
    _fake_pylon(monkeypatch, cam)
    burst, exposure, period = next(lc.pylon_bursts(9))
    assert period is None and exposure == pytest.approx(0.001)
    assert np.allclose([t for _, t in burst], 0.01 * np.arange(9))
    gen = lc.pylon_frames()
    f0, t0, _ = next(gen)
    f1, t1, _ = next(gen)
    assert t0 == 0 and t1 == pytest.approx(0.01) and f1[0, 0] == f0[0, 0] + 1
    gen.close()


def test_pylon_flash_frames_without_timestamps(lc, monkeypatch, capsys):
    frames = [np.full((40, 40), k, dtype=np.uint8) for k in range(20)]
    _fake_pylon(monkeypatch, _Camera(frames))                                # timestamp 0, as the emulator
    monkeypatch.setattr(lc, "time", _Clock(0.2))
    gen = lc.pylon_frames()
    _, t0, _ = next(gen)
    _, t1, _ = next(gen)
    gen.close()
    assert t0 == 0 and t1 == pytest.approx(0.2)
    assert "no usable timestamps; using the host clock" in capsys.readouterr().out


def test_pylon_failed_grabs_stop_with_a_message(lc, monkeypatch, capsys):
    frames = [np.zeros((40, 40), dtype=np.uint8)]
    _fake_pylon(monkeypatch, _Camera(frames, ok=lambda k: False))
    with pytest.raises(SystemExit, match="20 grabs in a row failed"):
        next(lc.pylon_frames())
    assert "grab failed (20 in a row): 0xe1000014 The buffer was incompletely grabbed." in capsys.readouterr().out
    _fake_pylon(monkeypatch, _Camera(frames, ts=lambda k: k + 1, ok=lambda k: k % 4 != 3))   # every 4th fails
    with pytest.raises(SystemExit, match="20 bursts in a row lost frames"):
        next(lc.pylon_bursts(9))


def test_pylon_pixel_formats(lc, monkeypatch):
    frames = [np.arange(1600, dtype=np.uint16).reshape(40, 40)]
    for fmt, scale in (("Mono12", 1), ("Mono12p", 16)):
        _fake_pylon(monkeypatch, _Camera(frames, ts=lambda k: k + 1, PixelFormat=_Node(fmt, writable=False)))
        f, _, _ = next(lc.pylon_frames())
        assert f.dtype == np.uint16 and np.array_equal(f, frames[0] * scale)
    cam = _Camera(frames, PixelFormat=_Node("BayerRG8", writable=False))
    _fake_pylon(monkeypatch, cam)
    with pytest.raises(SystemExit, match="BayerRG8"):
        next(lc.pylon_frames())
    assert not cam.open


def test_video_file_through_opencv(lc, tmp_path, capsys):
    cv2 = pytest.importorskip("cv2")
    path = str(tmp_path / "flash.avi")
    w = cv2.VideoWriter(path, cv2.VideoWriter_fourcc(*"FFV1"), 50, (128, 128), isColor=False)
    if not w.isOpened():
        pytest.skip("OpenCV cannot write FFV1 video here")
    for k in range(3):
        g = _uint8(simulate_pulsed_image(128, (3.0, -2.0), 4, rng=20 + k))
        w.write(g)
        cv2.imwrite(str(tmp_path / f"img_{k:04d}.png"), g)
    w.release()
    for source, extra, times in ((path, [], [0, 0.02, 0.04]),
                                 (str(tmp_path / "img_%04d.png"), ["--fps", "100"], [0, 0.01, 0.02])):
        assert lc.main(["--mode", "flash", "--source", source, "--roi", "128", "--flash-period", "0.002"] + extra) == 0
        rows = _table(capsys.readouterr().out)
        assert [r[0] for r in rows] == times
        for t, angle, speed, pxs, q, state in rows:
            assert state == "valid" and abs(speed / np.hypot(3, 2) - 1) < 0.01
            assert abs(pxs / (500 * np.hypot(3, 2)) - 1) < 0.01
