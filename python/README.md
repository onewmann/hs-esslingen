# fftvel

Python implementation of the single-image velocity estimate (2-D FFT + Radon
transform of a pulsed-exposure image). It gives the same results as the
MATLAB/Octave package `+fftvel` in `../matlab` (checked by `tests/parity`).

Needs Python 3.10 or newer and, for the editable install, pip 21.3 or newer.
From this folder:

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[plot]"
python -m fftvel simulate --vx 5 --vy -3 --plot fig.png
python -m fftvel estimate image.png --L 4
```

The core needs only numpy. Extras: `data` (scipy, h5py and Pillow: `.mat` and
image files), `plot` (matplotlib: `--plot` and figures), `video`
(opencv-python-headless: webcams and video files in the live example),
`examples` (`plot`, `data` and `video`), `test` (pytest and the test
dependencies). If `opencv-python` is already installed, leave out `video` and
`examples` and install `".[plot,data]"`: both OpenCV packages provide `cv2`
and overwrite each other's files.

```python
from fftvel import estimate_velocity, pulsed_from_stack

est = estimate_velocity(bp, L=4)        # bp: mean of 2L+1 equally spaced exposures
print(est.angle_deg, est.speed, est.valid)

bp, K = pulsed_from_stack(frames, L=4)  # mean of the first 9 frames of a stack
est = estimate_velocity(bp, exposures=K, duty=0.2, roi=512, clip_outliers=True)
```

Tests, from the repository root:

```bash
python -m pip install -e "./python[test]"
python -m pytest tests/python
```

See the main README for the method, all options, results, limitations and the
live-camera example.
