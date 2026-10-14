# Contributing

The estimator exists twice: as the Python package
[`python/src/fftvel`](python/src/fftvel) and as the MATLAB/Octave package
[`matlab/+fftvel`](matlab/+fftvel). The two are written line by line in the
same order and must give the same numbers for the same input, to $10^{-9}$.
Python is the reference: the parity results are computed with it, and MATLAB
and Octave have to reproduce them. This page says how the two are kept
identical and what a change to the algorithm involves.

## Which function is where

The public functions have the same names in both packages. Python takes the
options as keyword arguments (`duty=0.2`), MATLAB as name-value pairs
(`'Duty', 0.2`, matched without regard to case).
[`matlab/+fftvel/Contents.m`](matlab/+fftvel/Contents.m) lists the public
MATLAB functions (`help fftvel.Contents`), and `tests/matlab/test_contents.m`
checks that it names every one with its H1 line and the version of
`python/pyproject.toml` and `CITATION.cff`.

| Python | MATLAB/Octave | Purpose |
|---|---|---|
| `fftvel.estimate_velocity` (`core.py`) | `fftvel.estimate_velocity` | the estimator; returns an `Estimate` (Python) or a struct with the same fields |
| `fftvel.pulsed_from_stack` (`core.py`) | `fftvel.pulsed_from_stack` | mean of the first K frames of a stack |
| `fftvel.to_physical` (`core.py`) | `fftvel.to_physical` | pixels per step to length units per second |
| `fftvel.angle_difference` (`core.py`) | `fftvel.angle_difference` | difference of two directions modulo 180 deg |
| `fftvel.simulate_pulsed_image` (`simulate.py`) | `fftvel.simulate_pulsed_image` | simulated mean of K exposures (averaged frames, or a multi-flash image) with known motion |
| `fftvel.make_texture` (`simulate.py`) | `fftvel.make_texture` | random texture for the simulator |
| `fftvel.phase_correlation` (`phasecorr.py`) | `fftvel.phase_correlation` | shift between two frames, the independent reference |
| `fftvel.plot_estimate` (`plotting.py`) | `fftvel.plot_estimate` | diagnostic figure; needs `diagnostics=True` (`'Diagnostics', true`) |

The steps of the estimator are module functions of `core.py` in Python and
functions of the subpackage `+internal` in MATLAB. Neither set is part of the
public API.

| Python (`fftvel.core`) | MATLAB/Octave (`fftvel.internal`) |
|---|---|
| module constants (`DEFAULT_L`, `LOG_FLOOR`, ...) | `constants` |
| `_check_option` | `check_option` |
| keyword arguments of `estimate_velocity` | `parse_options` |
| `crop_square` | `crop_square` |
| `hann_window` | `hann_window` |
| `dilate` | `dilate` |
| `saturated_share` | `saturated_share` |
| `log_spectrum` | `log_spectrum` |
| `mean_projection` | `mean_projection` |
| `ripple_profile` | `ripple_profile` |
| `parabolic_offset` | `parabolic_offset` |
| `_window_autocorr` | `window_autocorr` |
| `dirichlet_template` | `dirichlet_template` |
| `box_template` | `box_template` |
| `_template_scores` | `template_scores` |
| `fit_speed` | `fit_speed` |
| `fit_blur` | `fit_blur` |
| `null_quality` | local function `null_quality` in `estimate_velocity.m` |
| `_upsampled_peak` (`phasecorr.py`) | local function `upsampled_peak` in `phase_correlation.m` |

Only in Python: the command line (`__main__.py`), the file readers of
`fftvel.io`, `fftvel.stream` and `python/examples/live_camera.py`. Only in
MATLAB: `matlab/examples/live_basler.m`.
The demos `demo_simulation` and `demo_testsequences` exist in both
`python/examples/` and `matlab/examples/`.

## Constants

Every number of the algorithm is a named constant, and each lives in three
places: the module constants of `python/src/fftvel/core.py`,
`matlab/+fftvel/+internal/constants.m` and
[`tests/constants.json`](tests/constants.json), which also gives the meaning
of each constant and the functions that use it. `tests/python/test_constants.py`
and `tests/matlab/test_constants.m` check that both packages hold exactly the
names and values of the JSON file. A new constant therefore goes into all
three files, and a changed value as well. The table in
[docs/method.md](docs/method.md#constants) is a copy of the JSON file, and
`tests/python/test_docs_constants.py` checks that it holds the same names,
values and meanings in the same order, so it changes with the JSON file.

## Parity cases

[`tests/parity/make_cases.py`](tests/parity/make_cases.py) runs the Python
implementation on 20 input images and writes three files:

* `cases.mat`: the input images;
* `expected.csv`: the Python results (angle, speed, GEBA, the qualities, the
  threshold, the blur length, the flags and the saturated share) and the
  options of each case;
* `fixed_inputs.mat`: images without motion at 32 and 64 px, a static
  repeating tile, and the hot-pixel image, the image of 8 exposures at
  0.6 px per step and the texture that the Python option and
  phase-correlation tests draw. They are stored because MATLAB and Octave
  draw different random numbers from the same seed and some draws fail these
  checks; the MATLAB tests load them instead of drawing them.

`tests/python/test_parity.py` and `tests/matlab/test_parity.m` estimate the
inputs of `cases.mat` and compare with `expected.csv`: angle, speed, GEBA,
quality, perpendicular quality and blur fit to $10^{-9}$, threshold and
saturated share to $10^{-12}$, and the flags `valid`, `isotropic` and
`blur_like` exactly. The cases cover odd sizes and FFT lengths, a non-square
image, padding factors 3 and 4, 8 and 10 exposures, a duty, a small `roi`,
clipped hot pixels and glints, a `uint8` frame stack, a `uint8` image with
saturated pixels, a static repeating tile, an image without motion and single
blurred exposures with and without `reject_blur`.
[`tests/parity/make_phase_cases.py`](tests/parity/make_phase_cases.py) does
the same for `phase_correlation` with 10 image pairs in `phase_cases.mat`.

Regenerate the files from the repository root, with the extras of
`./python[test]` installed:

```bash
python tests/parity/make_cases.py
python tests/parity/make_phase_cases.py
```

All cases draw from one random generator in turn, so a new case goes at the
end of `CASES`; the inputs of the earlier cases then stay the same. On
unchanged code `expected.csv` comes out byte-identical, and the `.mat` files
differ only in the creation date in their header, so commit them only when
their content has changed.

## Changing the algorithm

1. Change `python/src/fftvel` and `matlab/+fftvel` in the same commit, in the
   same order of operations, so that the two files can still be read side by
   side.
2. Add or change constants in `core.py`, `constants.m` and
   `tests/constants.json` together.
3. Regenerate the parity cases with `make_cases.py` (and `make_phase_cases.py`
   if `phase_correlation` changed). Never edit `expected.csv` or the `.mat`
   files by hand, and check that the differences in `expected.csv` are the
   ones the change should cause.
4. Run both suites. `tests/matlab/run_tests.m` runs every
   `tests/matlab/test_*.m` it finds.

   ```bash
   python -m pytest -q tests/python
   cd tests/matlab && octave --no-gui --eval run_tests
   ```

   With MATLAB, from the repository root:
   `matlab -batch "cd tests/matlab; run_tests"`.

5. Rerun the studies. `python experiments/run.py --all` rewrites the CSV files
   in `experiments/results/` and the figures in `docs/figures/` (its docstring
   gives about 70 min on 4 cores); `--quick` checks the scripts with smaller
   grids and writes to the ignored `experiments/results_quick/`.

   ```bash
   python experiments/run.py --quick --all
   python experiments/run.py --all
   python experiments/run.py --replot --all
   ```

   The last command redraws the figures from the CSV files and prints a
   summary of each study. Not every number in the documentation is in these
   summaries; the others are computed from the same CSV files.
6. Update every number in `README.md`, `docs/*.md` and
   `data/measurements/README.md` that the new CSV files change.

The MATLAB code has to run in GNU Octave 8.4 as well as in MATLAB, without
toolboxes. Five constructs do not work in Octave 8.4 (checked):

* a function in `+fftvel/private/` is undefined when called from a function
  of the package, hence `+internal` with fully qualified calls such as
  `fftvel.internal.fit_speed`;
* `import fftvel.*` is not implemented;
* `exist('fftvel.estimate_velocity')` returns 0;
* `arguments` blocks are a parse error;
* the call syntax `f(Name=Value)` of MATLAB R2021a and newer is read as an
  assignment, so only the value is passed; examples and tests give options as
  pairs (`'Duty', 0.02`).

GitHub Actions ([`.github/workflows/ci.yml`](.github/workflows/ci.yml)) runs
the Python suite on Python 3.10, 3.13 and 3.14 together with a small subset of
the studies, and the MATLAB/Octave suite and the demos in GNU Octave and in
MATLAB. MATLAB needs no licence in CI on a public repository.
