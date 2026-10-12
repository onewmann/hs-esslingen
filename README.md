# Surface velocity from a single multi-flash image

2-D FFT and Radon transform, Python and MATLAB/Octave

[![tests](https://github.com/onewmann/2D-FFT-Code/actions/workflows/ci.yml/badge.svg)](https://github.com/onewmann/2D-FFT-Code/actions/workflows/ci.yml)

A strobe flashes K = 2L+1 times (9 by default) during one camera exposure, so a moving textured surface appears K times in one image, shifted by the same step each time.
In the 2-D spectrum of that image the shifts become evenly spaced straight lines: their angle gives the direction of motion, their spacing the distance moved per flash.
The intended use is contactless measurement of a vehicle's speed over ground and sideslip angle, as in the thesis the method comes from, or of a conveyor belt or a paper web, with one image per measurement.

## Status

| | |
|---|---|
| Simulation and synthetic sequences | Works. In simulation (256 px, 9 flashes, no noise) the direction is within 0.70 deg and the speed within 0.47 % from 0.4 to 64 px per flash. On the two synthetic test sequences with constant motion all 24 windows are valid. |
| Real multi-flash images | Recordings are being prepared and are not in the repository yet. The 2025 camera recordings are frame sequences with uneven steps; the estimator rejects all 6 of their windows, as it should. |
| Camera scripts | Experimental. Neither live script has run on camera hardware; the Python one was checked with the pylon camera emulator. |
| Tests | CI with Python 3.10 and 3.13 and GNU Octave. A MATLAB job exists and runs for the first time once the repository is on GitHub; until then the code has not run in MATLAB. |

## Who did what

The method comes from Jan Horn's doctoral thesis ([reference](#reference)): the pulsed-exposure model, the log spectrum and Radon transform, the ripple profile WR for the direction, the speed from the spacing of the repeated structures and the quality measure GEBA.
In a student project in 2025 I implemented it in MATLAB; those scripts are kept unchanged at tag [`student-project-2025`](https://github.com/onewmann/2D-FFT-Code/tree/student-project-2025/2D-FFT).
Everything after commit `8f0aec4` was done with the assistance of AI (Claude by Anthropic), as [How this was built](#how-this-was-built) describes.

![One exposure of a real surface texture, the same surface moving under 9 flashes in one image, and the log spectrum of that image](docs/figures/hero.png)

*In simulation: the 9 copies in (b) turn into the lines in (c), perpendicular to the motion and 1/5.83 cycles per px apart.*

## Try it

From a clone of this repository, with Python 3.10 or newer (the editable install needs pip 21.3 or newer, hence the first line):

```bash
python -m pip install --upgrade pip
python -m pip install -e "./python[examples]"
python -m fftvel simulate --vx 5 --vy -3 --plot estimate.png
```

This simulates a 256 x 256 px image with 9 flashes, prints the true velocity and the estimate and saves the diagnostic figure. In MATLAB or GNU Octave, from the repository root and without any toolbox:

```matlab
addpath('matlab');
bp = fftvel.simulate_pulsed_image(256, [5 -3], 4);
est = fftvel.estimate_velocity(bp, 'L', 4)
```

## Why one image

A camera that takes separate frames measures motion against the frame period, 10 ms at 100 frames per second. A multi-flash image measures it against the flash period, which a controller can set in microseconds. Take a vehicle at 30 m/s, a ground sample distance (the length on the surface that one pixel covers) of 0.2 mm/px and a 256 px field, 51.2 mm wide. At 100 frames per second two consecutive frames are 300 mm apart and share no content; frames 20 px apart would need 7500 frames per second. Flashes 133 µs apart give 20 px per flash inside one exposure, well within the 0.4 to 64 px per flash tested for a 256 px square.

If you have separate frames that overlap, use phase correlation (the shift between two images from the phase of their cross-power spectrum). It was correct more often in simulation and also gives the sign of the motion, which one multi-flash image cannot. On the 840 noisy images of the noise study, in simulation:

| Method | Input | Correct, random texture | Correct, camera texture |
|---|---|---:|---:|
| this method (`fftvel`) | one multi-flash image | 80.2 % | 63.3 % |
| cepstrum / autocorrelation | one multi-flash image | 69.5 % / 72.1 % | 59.8 % / 55.0 % |
| phase correlation | two overlapping single exposures | 95.0 % | 87.9 % |

"Correct" means a direction within 2 deg and a speed within 5 %. Only `fftvel` also marks its results as valid or not. [results.md](docs/results.md#baselines) has the results per speed and the 2025 peak-spacing rule in corrected form.

## How it works

Between two flashes the surface moves by the step $\mathbf v$ pixels. The image is the mean of $K$ shifted copies of the texture $b$, so its spectrum is the spectrum $B$ of the texture times a Dirichlet kernel $D$:

$$B_p(\mathbf f) = B(\mathbf f)\,D(\mathbf f\cdot\mathbf v),\qquad D(u) = \frac{\sin(K\pi u)}{K\sin(\pi u)}$$

The magnitude of $D$ is 1 wherever $\mathbf f\cdot\mathbf v$ is an integer, so the magnitude spectrum carries straight lines perpendicular to $\mathbf v$, spaced $1/|\mathbf v|$ cycles per pixel apart. The estimator finds them in four steps:

1. Spectrum. The image is multiplied by a 2-D Hann window and its FFT zero padded to twice the size. [thesis] The estimator works on the logarithm of the magnitude, which turns the product $B\cdot D$ into a sum. [rework] A non-square image is cropped to a central square, never resized, so both axes keep the camera's pixel scale.
2. Direction [thesis]. A Radon transform projects the log spectrum along parallel lines for every angle. Along the direction of motion the projection becomes a sharp comb, and the ripple profile WR (the sum of absolute differences between neighbouring bins of a projection) peaks there. The thesis's quality measure GEBA $= \max WR/\mathrm{mean}\,WR$ is reported. [rework] The projection averages along each line instead of summing, which removes a bias towards 45 and 135 deg.
3. Speed [rework]. The thesis takes the speed from the spacing of the repeated structures. The rework fits a closed-form model of the whole projection for $K$ flashes, including the autocorrelation of the Hann window, on a grid of speeds 0.5 % apart. The speed with the best correlation wins, and that correlation is the fit quality.
4. Validity [rework]. An estimate is valid if the fit quality reaches $\max\big(0.6,\ \tanh(9/\sqrt{m-5})\big)$ for a profile of $m$ independent points, the speed lies inside the measurable range (`in_range`) and the spectrum is not isotropic, as it is for defocus or a square repeating texture.

The options `exposures`, `duty`, `roi`, `clip_outliers` and `reject_blur` come from the rework as well. [docs/method.md](docs/method.md) gives the full model, the threshold calibration, all options and the 33 constants.

## Results

Ordered by how independent the data are of the estimator's own model. All numbers come from `experiments/results/*.csv`.

| Study | Data | Result |
|---|---|---|
| [Synthetic test sequences](docs/results.md#synthetic-test-sequences) | 4 sequences from the 2025 student project, made without `fftvel`; windows of 9 frames | Constant motion: 24 of 24 windows valid, 134.93 to 135.03 deg (true 135), speed within 0.005 px. Changing speed: 17 of 22 valid, errors up to 11.1 %. Irregular steps: 4 of 22 valid, errors up to 14.4 % |
| [2025 recordings](docs/results.md#the-2025-recordings) | 3 Basler frame sequences of a moving surface | Steps vary by a factor of 2.8 to 4.0 within a recording; all 6 windows rejected (quality 0.23 to 0.40, threshold 0.6) |
| [Single camera frames](docs/results.md#single-camera-frames-a-real-null-test) | 1360 crops of single frames of the 2025 recordings: one exposure each, so nothing repeats | 59 valid, all false alarms; 41 of the 47 on the moving surface look like motion blur |
| [Accuracy](docs/results.md#accuracy-and-noise-in-simulation) | in simulation: 256 px, 9 flashes, no noise, random and camera texture | Within 0.70 deg and 0.47 % from 0.4 to 64 px per flash; 0 of 24 valid at 0.3 px per flash, below the range |
| [Noise](docs/results.md#accuracy-and-noise-in-simulation) | in simulation: 60 images per noise level and texture | Noise at half the signal: 91.7 % correct (random texture) and 70.0 % (camera texture); noise as strong as the signal: 71.7 % and 50.0 %. 11 of 592 valid estimates wrong; 0 of 120 images without motion valid |
| [Baselines](docs/results.md#baselines) | in simulation: the 840 noise images | Correct (random / camera texture): `fftvel` 80.2 / 63.3 %, single-image baselines at most 72.1 / 59.8 %, phase correlation of two frames 95.0 / 87.9 % |

The simulated images come from the model the estimator fits. They check the implementation and say nothing about how well that model describes a real camera. Most numbers in the next section come from an independent forward model, which does not use the estimator's simulator and adds sensor noise, flash jitter and other violations of the model ([results.md](docs/results.md#model-violations-in-an-independent-forward-model)).

## Where it fails

Silent failures first: these give a wrong result that is still marked valid. All numbers are in simulation unless the row says otherwise.

| Cause | What happens | Remedy |
|---|---|---|
| Long flashes, `duty` not set | flashes of 0.25 and 0.5 of the period: 12.5 % and 18.8 % wrong but valid, speed near duty x speed / K | short flashes, or `duty` / `--duty`; with the true duty all images were correct and valid ([duty](docs/hardware.md#flash-duration-and-duty)) |
| Wrong flash count | 7, 8, 10 or 11 assumed for 9: 34.4 % wrong but valid, speed scaled by 9/K | the exact K with `L` or `exposures`; count the pulses once with a photodiode ([flash count](docs/hardware.md#the-exact-flash-count)) |
| Continuous light, or motion blur at slow speed | real single camera frames: 59 of 1360 crops valid. Continuous light for 3/4 of the signal: 34.4 % wrong but valid | dark enclosure or band-pass filter, control images with one flash; `reject_blur`, which costs correct results below about 2.2 px per flash ([ambient light](docs/hardware.md#ambient-light)) |
| Glints, saturated spots | two static glints of 3 times the mean signal: 96.9 % wrong but valid | `clip_outliers`, with which all were correct and valid; diffuse light, no saturated pixels ([glints](docs/hardware.md#saturation-and-glints)) |
| Strong static edges | an edge of twice the mean signal: 28.1 % wrong but valid | keep edges and masks out of the analysed square (`roi`), fill masked areas with the image mean ([edges](docs/hardware.md#saturation-and-glints)) |
| Periodic texture | a pattern repeating every 64 or 128 px: 100 % wrong but valid; every 160 px: 50.0 % | a non-repeating texture ([assumptions](docs/method.md#assumptions)) |
| Uneven flash intervals, changing speed | intervals alternating between 0.8 and 1.2 periods: 87 of 240 valid at twice the speed; one missing flash: 22 of 240 valid but wrong; accelerating test sequence: 17 of 22 windows valid, errors up to 11.1 % | flashes from a hardware timer, period checked with a photodiode ([flash count](docs/hardware.md#the-exact-flash-count)) |

These are flagged (marked not valid) in most cases:

| Cause | What happens | Remedy |
|---|---|---|
| Motion below the measurable range, about 3.35/K px per flash | 0 of 24 valid at 0.3 px per flash with 9 flashes | longer flash period or more flashes ([flash period](docs/hardware.md#choosing-the-flash-period)) |
| Heavy noise | camera texture with noise as strong as the signal: 50.0 % correct, 50.0 % valid; over all noise levels 11 of 592 valid estimates wrong | more light, more pixels per flash ([flash period](docs/hardware.md#choosing-the-flash-period)) |
| No texture or no motion | 0 of 120 images without motion valid, and 0 of 2268 simulated ones in the threshold calibration; real crops of a static room scene: 4 of 252 valid | keep saturated areas out of the image ([glints](docs/hardware.md#saturation-and-glints)) |
| Defocus | the isotropy check rejects the ring pattern of a static defocused surface; with sensor noise, correct and valid falls from 84.4 % in focus to 40.6 % at a blur radius of 2 px | focus at the working distance ([focus](docs/hardware.md#focus)) |

The direction is known only modulo 180 deg: one image cannot tell $\mathbf v$ from $-\mathbf v$.

## Use it

Python. The options after the image are keyword-only; the MATLAB name-value pairs correspond one to one ([options](docs/method.md#options-and-when-to-use-them)).

```python
from fftvel import estimate_velocity, simulate_pulsed_image, to_physical

bp = simulate_pulsed_image(256, (5, -3), 4)     # 9 flashes, 5 px right and 3 px up per flash
est = estimate_velocity(bp, L=4)                # or exposures=K for any flash count
print(est.angle_deg, est.speed, est.quality, est.valid)
p = to_physical(est, gsd=0.05, period=0.001)    # 0.05 mm per px, 1 ms between flashes
print(p.speed)                                  # mm/s
```

`fftvel.pulsed_from_stack` averages a stack of frames into one image. The command line reads PNG, TIFF, BMP, `.npy` and `.mat` files, writes the result as JSON to stdout and notes on what it read to stderr. A single image needs `--L` or `--exposures`. With `--gsd` (length per pixel) and `--period` (seconds between flashes) the JSON also holds `speed_per_s`, `vx_per_s` and `vy_per_s`: 200.11 mm/s for this synthetic image, made with 200 mm/s.

```bash
python -m fftvel estimate data/real/example/flash_a.png --L 4 --duty 0.02 --gsd 0.05 --period 0.001
```

Further options are `--roi` (side of the analysed square, default 512 px), `--clip-outliers`, `--stack` and `--plot`; `python -m fftvel estimate --help` lists them all, and [python/README.md](python/README.md) explains the install extras. The same estimate in MATLAB or GNU Octave:

```matlab
addpath('matlab')                                        % the folder that holds +fftvel
bp = double(imread('data/real/example/flash_a.png'));
est = fftvel.estimate_velocity(bp, 'L', 4, 'Duty', 0.02, 'Diagnostics', true);
p = fftvel.to_physical(est, 0.05, 1e-3);                 % 0.05 mm per px, 1 ms between flashes
fprintf('%.1f mm/s at %.1f deg, valid %d\n', p.speed, est.angle_deg, est.valid);
fftvel.plot_estimate(bp, est);
```

The live scripts are experimental. [`python/examples/live_camera.py`](python/examples/live_camera.py) streams from a Basler camera through pypylon or from a webcam or video file through OpenCV, [`matlab/examples/live_basler.m`](matlab/examples/live_basler.m) from a Basler camera through the Image Acquisition Toolbox. Neither has run on camera hardware. The Python script was checked with the pylon camera emulator and has 13 tests with fake camera modules; `live_basler.m` has not run in MATLAB. Their flash mode estimates every multi-flash frame on its own, while their average mode, which averages K consecutive frames under continuous light, is a test harness. At the default 512 px square one estimate takes about 4 s in Python ([live scripts](docs/hardware.md#live-scripts)).

To evaluate your own recordings, list the images in a `meta.csv` with flash count, flash period, duty, ground sample distance and the set speed and direction, and add images of the stopped surface as controls. [`data/real/README.md`](data/real/README.md) describes the format and the procedure, [docs/hardware.md](docs/hardware.md) what the set-up needs. On the synthetic example the tool marks both moving images valid with speed errors below 0.06 % and the static image not valid:

```bash
python experiments/real_measurements.py data/real/example/meta.csv --out results.csv
```

## Development

```bash
python -m pip install -e "./python[test]"
python -m pytest -q tests/python
cd tests/matlab && octave --no-gui --eval run_tests
```

The tests check the implementation: known motion in simulation and on the test sequences, rejection of images without motion, every option and its argument checks, and parity, which means that the Python and MATLAB/Octave packages give the same results to 1e-9 on 20 shared cases (so far checked in GNU Octave). They do not show that the model holds for a given camera; that takes real multi-flash images. Because the two languages are compared case by case, Python and MATLAB change together; [CONTRIBUTING.md](CONTRIBUTING.md) lists the steps.

```
python/src/fftvel/   Python package; python -m fftvel is the command line
python/examples/     demos and the live camera script
matlab/+fftvel/      MATLAB/Octave package, helpers in +internal/; matlab/examples/ the same demos
tests/               both test suites, the parity cases and constants.json
experiments/         one script per study, run.py runs them; results/ holds the CSV files
data/                test sequences, a camera texture, the 2025 recordings, real/ for new recordings
docs/                method, results, acquisition, history of the code; figures/
```

A full clone is about 120 MB, because the history also holds the files of the 2025 student project. A partial clone fetches only the files of the version you check out (about 70 MB, most of it the four recordings in `data/recordings/`) and older files when they are needed:

```bash
git clone --filter=blob:none https://github.com/onewmann/2D-FFT-Code
```

## How this was built

The scripts of the 2025 student project are kept unchanged at tag [`student-project-2025`](https://github.com/onewmann/2D-FFT-Code/tree/student-project-2025/2D-FFT), which points to commit `8f0aec4`. The tag also holds a Basler camera driver for MATLAB by another author (MIT licence), which this repository does not use. Everything after commit `8f0aec4` was done with the assistance of AI (Claude by Anthropic): first a code review of those scripts, then the rework in this repository, starting with the commit "Rewrite the 2D-FFT velocity estimator for MATLAB, Octave and Python". That covers the code in `matlab/` and `python/`, the tests, the experiments, the figures and this documentation. [docs/history.md](docs/history.md) indexes the 2025 scripts by function, lists what the review found and compares the 2025 estimator with the current one.

## Reference

Jan Horn: *Zweidimensionale Geschwindigkeitsmessung texturierter Oberflächen mit flächenhaften bildgebenden Sensoren* (two-dimensional velocity measurement of textured surfaces with area image sensors). Dissertation, Universität Karlsruhe (TH), 2006. Schriftenreihe des Instituts für Mess- und Regelungstechnik, Band 006, Universitätsverlag Karlsruhe, 2007. ISBN 978-3-86644-076-0. [Full text (CC BY-NC-ND 2.0 DE)](https://publikationen.bibliothek.kit.edu/1000006276)

## Citation and licence

[CITATION.cff](CITATION.cff) holds the citation of this code, and GitHub offers it under "Cite this repository". For the method, cite the thesis above. The code is under the MIT licence, see [LICENSE](LICENSE). The test sequences, the camera texture and the recordings in `data/` come from the 2025 student project; [data/README.md](data/README.md) describes them.
