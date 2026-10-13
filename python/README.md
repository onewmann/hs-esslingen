# fftvel

The Python package of
[Surface velocity from a single multi-flash image](../README.md): direction
and speed of a moving textured surface from one multi-flash image, with the
2-D FFT and the Radon transform. It gives the same results as the
MATLAB/Octave package `+fftvel` in [`../matlab`](../matlab), to 1e-9 on the
shared cases of `tests/parity`, so far checked in GNU Octave; the code has not
run in MATLAB yet. The main README describes the method, the
results, where it fails and who did what.

## Install

Python 3.10 or newer. The editable install needs pip 21.3 or newer, which
older Python versions do not bring. From the repository root:

```bash
python -m pip install --upgrade pip
python -m pip install -e "./python[plot,data]"
```

Without a clone, pip installs the package straight from GitHub; the package
lies in the subfolder `python`, so the URL needs `#subdirectory=python`:

```bash
python -m pip install "fftvel @ git+https://github.com/onewmann/2D-FFT-Code#subdirectory=python"
```

This still downloads the files of the current version, about 70 MB, most of
it the recordings in `data/recordings`, although the package itself is about
100 KB.

The core (`python -m pip install -e ./python`) needs only NumPy and reads
`.npy` and `.npz` files. The extras add:

* `data`: SciPy, h5py and Pillow, for `.mat` and image files;
* `plot`: matplotlib, for `--plot` and the figures of the experiments;
* `video`: `opencv-python-headless`, for webcams and video files in the live
  script;
* `examples`: `plot`, `data` and `video` together;
* `test`: pytest and the test dependencies.

If `opencv-python` is already installed, for example for other scripts that
open windows, keep it and leave out `video` and `examples`: both OpenCV
packages provide the module `cv2` and overwrite each other's files.

## Command line

The install adds the command `fftvel`, a short form of `python -m fftvel`.
An editable install made before the command existed needs
`python -m pip install -e ./python` once more. From the repository root:

```bash
fftvel simulate --exposures 9 --vx 5 --vy -3 --plot estimate.png
fftvel estimate data/measurements/synthetic_flash/flash_b.png --exposures 8 --duty 0.2 --gsd 0.05 --period 0.0005
fftvel estimate data/recordings/vid2.mat --stack --exposures 9
```

`simulate` estimates a simulated 256 x 256 px image with 9 flashes and known
motion. `estimate` reads PNG, TIFF, BMP, `.npy`, `.npz` and `.mat` files,
writes the result as JSON to stdout and notes on what it read and used to
stderr. A single image needs the flash count, `--exposures K` for any K of
at least 2, or `--L L` for K = 2L+1 (L is the half-count of the thesis).
`--duty` is the flash duration divided by the flash period. `--gsd` (length
per pixel) and `--period` (seconds between flashes) add the speed per second
to the JSON: 559.68 mm/s for the second command, an image made with
560 mm/s. `--stack` averages the first K frames of a stack; the 2025
recording in the third command has uneven steps and comes out not valid.
Further options are `--roi` (default 512 px), `--clip-outliers`,
`--frames-axis`, `--var` and `--plot`; `fftvel estimate --help` explains them
all, and `fftvel --version` prints the version.

## Reading files

The module `fftvel.io` reads the files for the command line, and scripts can
use it as well:

```python
import fftvel
from fftvel.io import load

frames = load("data/measurements/synthetic_average/moving_regular.npz", stack=True)   # 192 x 192 x 9, uint8
bp, K = fftvel.pulsed_from_stack(frames, exposures=9)
est = fftvel.estimate_velocity(bp, exposures=K)
```

`load(path, var=None, stack=False)` returns the image or the frame stack in
its own data type. Pillow reads PNG, TIFF, BMP and JPEG and keeps 16-bit,
32-bit integer and float32 data; colour becomes grey (`to_grey`, with the
weights of MATLAB `rgb2gray`), a multi-page TIFF a stack, and float64 TIFF
and 16-bit colour TIFF need tifffile. Pillow reads only 8 bits of 16-bit
colour, so 16-bit colour in other formats (PNG, PPM) gives a `ValueError`;
save it as 16-bit grey, as TIFF or as `.npy`. `.npy` files are memory-mapped,
so only the frames used are read. In a `.npz` file, such as those that the
live script writes with `--save`, the array `frames` is read; `var` names
another one, in a `.npz` or a `.mat` file. `.mat` files need SciPy, MATLAB
v7.3 files also h5py. Octave's default `save` writes text, which is not a
MAT-file: use `save('-v7', file, 'frames')`, with structs and cells saved as
plain arrays.
`frames_axis(shape)` tells whether the first or the last axis of a stack
holds the frames. Errors are `FileNotFoundError`, `ImportError` (a reader
that the file needs is not installed) and `ValueError`; what was read is
logged at level INFO on the logger `fftvel.io`.

Use lossless formats for measurements: JPEG adds an 8 x 8 px block grid that
can give wrong results marked valid. `load` gives a `UserWarning` for JPEG
files, and the command line prints it.

## Tests

From the repository root:

```bash
python -m pip install -e "./python[test]"
python -m pytest -q tests/python
```

The test of video files in the live script needs OpenCV and is skipped
without it; `./python[test,video]` installs it too, as the CI does.
