# fftvel

The Python package of
[Surface velocity from a single multi-flash image](../README.md): direction
and speed of a moving textured surface from one multi-flash image, with the
2-D FFT and the Radon transform. It gives the same results as the
MATLAB/Octave package `+fftvel` in [`../matlab`](../matlab), to 1e-9 on the
shared cases of `tests/parity`, so far checked in GNU Octave; the code has not
run in MATLAB yet. The main README describes the method, the
results, where it fails and how this code was built.

## Install

Python 3.10 or newer. The editable install needs pip 21.3 or newer, which
older Python versions do not bring. From the repository root:

```bash
python -m pip install --upgrade pip
python -m pip install -e "./python[examples]"
```

The core (`python -m pip install -e ./python`) needs only NumPy and reads
`.npy` files. The extras add:

* `data`: SciPy, h5py and Pillow, for `.mat` and image files;
* `plot`: matplotlib, for `--plot` and the figures of the experiments;
* `video`: `opencv-python-headless`, for webcams and video files in the live
  script;
* `examples`: `plot`, `data` and `video` together;
* `test`: pytest and the test dependencies.

If `opencv-python` is already installed, for example for other scripts that
open windows, install `"./python[plot,data]"` instead and keep it: both
OpenCV packages provide the module `cv2` and overwrite each other's files.

## Command line

From the repository root:

```bash
python -m fftvel simulate --vx 5 --vy -3 --plot estimate.png
python -m fftvel estimate data/real/example/flash_b.png --exposures 8 --duty 0.2 --gsd 0.05 --period 0.0005
python -m fftvel estimate data/recordings/vid2.mat --stack --L 4
```

`simulate` estimates a simulated 256 x 256 px image with 9 flashes and known
motion. `estimate` reads PNG, TIFF, BMP, `.npy` and `.mat` files, writes the
result as JSON to stdout and notes on what it read and used to stderr. A
single image needs the flash count, `--L L` for 2L+1 flashes or
`--exposures K` for any K of at least 2. `--duty` is the flash duration
divided by the flash period. `--gsd` (length per pixel) and `--period`
(seconds between flashes) add the speed per second to the JSON: 559.68 mm/s
for the second command, an image made with 560 mm/s. `--stack` averages the
first 2L+1 (or K) frames of a stack; the 2025 recording in the third command
has uneven steps and comes out not valid. Further options are `--roi`
(default 512 px), `--clip-outliers`, `--frames-axis`, `--var` and `--plot`;
`python -m fftvel estimate --help` explains them all.

Files. Pillow reads PNG, TIFF, BMP and JPEG and keeps 16-bit, 32-bit integer
and float32 data; colour becomes grey, a multi-page TIFF a stack, and float64
TIFF needs tifffile. `.npy` files are memory-mapped, so only the frames used
are read. `.mat` files need SciPy, MATLAB v7.3 files also h5py. Octave's
default `save` writes text, which is not a MAT-file: use
`save('-v7', file, 'frames')`, with structs and cells saved as plain arrays.
Use lossless formats for measurements: JPEG adds an 8 x 8 px block grid that
can give wrong results marked valid, and the command line warns for JPEG
input.

Tests, from the repository root:

```bash
python -m pip install -e "./python[test]"
python -m pytest -q tests/python
```
