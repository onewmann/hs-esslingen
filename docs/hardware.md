# Acquisition

What a camera set-up needs so that one image satisfies the model of
[method.md](method.md). No real multi-flash image has been recorded for this
repository yet. The numbers on this page come from the simulations and the
independent forward model in [results.md](results.md), from the 2025 camera
recordings, and from worked examples computed here.

## Timing of one image

A controller fires K flashes of duration $\tau$ at a period $T$ (start to
start) while the camera exposes one frame. The exposure has to begin before
the first flash and end after the last one, so it lasts at least
$(K-1)\,T + \tau$:

```
exposure   |<------------------ at least (K-1)*T + tau ------------------>|
flashes    |##|        |##|        |##|            ...              |##|
           tau   <----- T ----->                                (K flashes)
```

Between two flashes the surface moves by $`v_\text{px}`$ pixels, the speed the
estimator measures (px per flash). The duty $d = \tau/T$ is the fraction of
the period during which a flash is on; during each flash the surface moves by
$`d\,v_\text{px}`$ px, which blurs every copy along the motion.

The safest wiring lets the camera start a burst of exactly K pulses (an
exposure-active or frame-trigger output to the flash controller), or lets the
controller trigger the camera. With a free-running strobe the number of
flashes inside an exposure changes from frame to frame.

## Choosing the flash period

The motion per flash in pixels is

$$v_\text{px} = \frac{v\,T}{\mathrm{GSD}}$$

with $v$ the speed of the surface (mm/s), $T$ the flash period (s) and GSD the
ground sample distance, the length on the surface that one pixel covers
(mm/px). The estimator measures from about $3.35/K$ px per flash (0.37 px for
K = 9, 0.08 px for K = 41). Its speed grid reaches $N/2$ for an $N$ px square,
but the studies only test up to $N/4$, so treat $N/4$ as the upper end (128 px
for a 512 px square). Within that range, more pixels per flash help when there
is noise: in the noise study the share of correct estimates rose from 28.0 %
at 1 to 2 px per flash to 82.0 % at 4 to 8 px and 98.1 % at 16 to 32 px, in
simulation. Below about 2.2 px per flash a single blurred exposure and the
flash pattern cannot be told apart (blur check in
[method.md](method.md#blur-check-optional)). A few pixels per flash, about 4 to
10, is a good target, so $`T = v_\text{px}\,\mathrm{GSD}/v`$.

Worked example. A surface moves at 500 mm/s under a camera with
GSD = 0.05 mm/px (the 1936 px wide sensor then covers 96.8 mm). For 5 px per
flash, $T = 5 \cdot 0.05 / 500$ s = 0.5 ms, a flash rate of 2 kHz. With K = 9
the exposure must cover $(K-1+d)\,T$; with 20 µs flashes, $d$ = 0.04 and the
exposure needs at least 8.04 x 0.5 ms = 4.02 ms. The blur per flash is
$`d\,v_\text{px}`$ = 0.2 px, negligible. With the same T, speeds from 100 to
1000 mm/s give 1 to 10 px per flash; the lower end of 0.37 px corresponds to
37 mm/s. During the 9 flashes the surface moves 40 px, a small part of a
512 px square.

The same arithmetic shows why one image is used at all. At 30 m/s with a GSD
of 0.2 mm/px, a 256 px field is 51.2 mm wide. At 100 frames per second two
consecutive frames are 300 mm apart and do not overlap, while flashes 133 µs
apart give 20 px per flash inside one exposure.

## Flash duration and duty

Keep the flashes short, so that $`d\,v_\text{px}`$ stays well below about 3 px,
or give the duty to the estimator (`duty`, `--duty`). In simulation with the
forward model, flashes of a quarter and half of the period without the option
made 12.5 % and 18.8 % of the images wrong but valid, with the speed near
$`d\,v_\text{px}/K`$. With `duty` set to the true value, every image from duty
0 to 0.9 was correct and valid. Near duty 0.5, for example with a 50 % square
wave, every second main line of the pattern vanishes, so there the value must
be exact.

## The exact flash count

The estimator needs the true K (`L` for 2L+1 flashes, or `exposures` for any K,
also even). For slow motion a wrong K scales the speed by $`K_\text{true}/K`$
without lowering the fit quality. In simulation with the forward model, telling
the estimator 7, 8, 10 or 11 flashes for an image with 9 made 34.4 % of the
images wrong but valid, with the speed off by the factor 9/K.

The flashes must also be equally spaced. Intervals alternating between 0.8 and
1.2 periods gave 87 of 240 images valid at twice the true speed; one missing
flash gave 22 of 240 valid but wrong, in simulation. Random jitter of a clocked
trigger is less critical: in simulation with the forward model, a standard
deviation of 1 % of the period gave a median speed error of 0.15 %, and wrong
but valid results appeared only from 10 %. Check the number of pulses and the
period once with a photodiode and an oscilloscope. With the flash period given,
both live scripts and `experiments/real_measurements.py` warn when the exposure
is shorter than $(K-1+d)\,T$.

## Ambient light

The exposure lasts about $(K-1)\,T$, and any light other than the flashes adds
a motion-blurred image of the surface. A blur of $b$ px along one direction
looks like K flashes $b/K$ px apart and can be marked valid. This is the main
risk found so far on real images: single frames of the moving surface in the
2025 recordings (one exposure each, no flashes) came out valid in 19 and 28 of
340 crops of `vid3` and `vid4`. The valid crops of `vid4` point at about 90
deg, with K times the speed at a median of 4.67 px, a blur of about 4.7 px
([results.md](results.md#single-camera-frames-a-real-null-test)). In simulation
with the forward model, continuous light over an exposure of 10 flash periods
made 3.1 % of the images wrong but valid when it supplied half of the signal,
and 34.4 % from three quarters on. A constant offset such as flare only adds
shot noise: at three quarters of the signal, 43.8 % of the images were still
correct and valid.

Darken the room, or use an enclosure or a band-pass filter for the flash
colour, and check that an image taken with the flashes switched off is black.
Record controls that must come out not valid: images of the stopped surface,
and images of the moving surface with a single flash in the same exposure
time. A valid control shows blur from other light or from vibration.
`reject_blur` rejects blur-like results, at the cost of many correct ones
below about 2.2 px per flash.

## Focus

A defocused image has ring-shaped zeros in its spectrum, which look like slow
motion. The isotropy check is there to reject that pattern when nothing moves,
and both test suites check it on a defocused static texture. With motion,
defocus mostly costs valid results. In simulation with the forward model and
without noise, a disc blur of up to 8 px radius gave no wrong results marked
valid, and 31 or 32 of the 32 images per level were correct and valid; with the
noise of a 32000 e- sensor, the share of correct and valid images fell from
84.4 % in focus to 40.6 % at 2 px radius and 0 % at 5 px. Focus on the surface
at the working distance of the measurement, and check that an image of the
stopped surface comes out not valid.

## Saturation and glints

Static bright structures do not move with the surface and give fringes of their
own. In simulation with the forward model, two static glints with a peak of the
mean signal made 28.1 % of the images wrong but valid, and 96.9 % from three
times the mean on. `clip_outliers` (`--clip-outliers`) removed the effect
there: all images were correct and valid. A sharp static edge of 1.5 times the
mean signal gave 12.5 % wrong but valid results. Of the 12 valid single-frame
crops of the static room scene `vid1`, 7 lie in crops with at least 10 %
saturated pixels.

Set the exposure or the flash energy so that no pixel saturates; the estimate
reports the share of saturated pixels in `saturated`, and the command line and
the scripts warn above 1 %. Shiny surfaces need diffuse light or polarisers. Do
not expose too dark either: 8-bit quantisation caused no failures in simulation
with the forward model, 6 bits gave 6.2 % wrong but valid results, and 8 bits
at a quarter of full scale are about 6 bits. Store images in a lossless format
(PNG, TIFF, BMP, `.npy`, `.mat`); JPEG adds an 8 x 8 px block grid.

## Ground sample distance and physical units

The estimator measures pixels per flash. The speed in physical units is

$$v = \frac{v_\text{px}\,\mathrm{GSD}}{T}$$

`fftvel.to_physical(est, gsd, period)` does this in Python and MATLAB; the
command line takes `--gsd` and `--period` and then adds `speed_per_s`,
`vx_per_s` and `vy_per_s` to its JSON output, in the length unit of the GSD
per second:

```bash
python -m fftvel estimate data/real/example/flash_b.png --exposures 8 --duty 0.2 --gsd 0.05 --period 0.0005
```

```python
from fftvel import estimate_velocity, simulate_pulsed_image, to_physical
bp = simulate_pulsed_image(256, (5, -3), 4)        # or your own 2-D image
est = estimate_velocity(bp, L=4)
p = to_physical(est, gsd=0.05, period=0.0005)      # mm/s for the GSD in mm/px
```

`to_physical` converts every estimate, valid or not, so check `est.valid`
first.

Measure the GSD with a ruler or a checkerboard in the plane of the moving
surface, at the working distance, focus and zoom of the measurements, near the
centre of the image where the analysed square lies. Two marks far apart, each
read to within 1 px, give the GSD to 0.09 % when they are 1600 px apart and to
0.88 % at 160 px. Measure along x and y: a difference of more than about
0.5 % means that the camera looks at the surface at an angle. The axis that
moves the surface can measure it too: move it by a known distance between two
images of the stopped surface and take the displacement in pixels from
`fftvel.phase_correlation`. [`data/real/README.md`](../data/real/README.md)
gives both procedures step by step.

## Error budget

Each factor of $`v_\text{px}\,\mathrm{GSD}/T`$ adds its relative error to the
speed in mm/s; independent errors add in quadrature.

| Source | How it enters | Size |
|---|---|---|
| Estimator | relative error of $`v_\text{px}`$ | without noise, median 0.018 % and largest 0.47 % from 0.4 to 64 px per flash (K = 9, 256 px), in simulation; noise and model violations in [results.md](results.md) |
| GSD | one to one | 0.09 % with marks 1600 px apart read to 1 px, 0.88 % with marks 160 px apart |
| Flash period, systematic | one to one | measure it with a photodiode and an oscilloscope |
| Flash timing jitter, per flash | random error of $`v_\text{px}`$ | median speed error of the valid estimates 0.15, 0.35, 0.62, 1.38 and 3.16 % for a jitter of 1, 2, 5, 10 and 20 % of the period; wrong but valid results from 10 % on (forward model, in simulation) |
| Working distance with an entocentric lens | one to one through the GSD | 1.06 % for a height change of 3 mm at 300 mm distance with a 16 mm lens |
| Flash count | factor $`K_\text{true}/K`$ | for slow motion, 12.5 % if 9 flashes are taken for 8; not detected |

An entocentric lens is an ordinary lens: its magnification changes with the
object distance. For a thin lens of focal length $f$ at object distance $Z$,
$\mathrm{GSD} = p\,(Z-f)/f$ with $p$ the pixel pitch, so a change $\Delta Z$
changes the GSD by $\Delta Z/(Z-f)$. On a vehicle the body moves on its
suspension, and the height has to be measured or kept constant. A telecentric
lens keeps the magnification constant within its depth of field and removes
this term.

## Live scripts

Two scripts stream from a camera and estimate every image. Both need an
explicit mode:

* `flash`: every frame holds K strobe flashes and is estimated on its own; the
  speed is in px per flash, and in px/s with the flash period.
* `average`: continuous light; K consecutive frames are averaged into a pulsed
  image, and both scripts check that the frame intervals are equal. The 2025
  recordings show why this check matters
  ([history.md](history.md#what-the-review-found)). This mode is a test
  harness: frames under continuous light are blurred over their exposure, and
  the use case is the flash mode.

```bash
python python/examples/live_camera.py --mode flash --source pylon --L 4 --flash-period 0.0005 --duty 0.04
python python/examples/live_camera.py --mode average --source pylon --exposure-us 500
python python/examples/live_camera.py --mode flash --source recording.avi
```

```matlab
live_basler('Mode', 'flash', 'L', 4, 'FlashPeriod', 5e-4, 'Duty', 0.04)
live_basler('Mode', 'average', 'Exposure', 500)
```

[`python/examples/live_camera.py`](../python/examples/live_camera.py) reads
Basler cameras through pypylon, and webcams and video files through OpenCV.
[`matlab/examples/live_basler.m`](../matlab/examples/live_basler.m) uses the
Image Acquisition Toolbox with the GenTL adaptor. Both default to a 512 px
square, which takes about 4 s per estimate in Python (repeated call,
[results.md](results.md#run-time)), so they give a reading every few seconds.

Test status: neither script has run on camera hardware. The Python script is
covered by 13 tests with fake pypylon and OpenCV modules
(`tests/python/test_live_camera.py`), and its pylon path was checked with the
pylon camera emulator. `live_basler.m` has no automated test, because GNU
Octave has no Image Acquisition Toolbox, and it has not run in MATLAB.

## Evaluating recordings

[`experiments/real_measurements.py`](../experiments/real_measurements.py)
estimates every image listed in a `meta.csv`, converts the results to mm/s and
compares them with the set motion.
[`data/real/README.md`](../data/real/README.md) describes the columns (file,
flash count, flash period, duty, GSD, set speed and direction, static flag,
optional single-flash reference frames), the angle convention, the static
controls with their Clopper-Pearson bound, and how to report the result. Try it
on the synthetic example first:

```bash
python experiments/real_measurements.py data/real/example/meta.csv --out results.csv
```

On that example (simulated images with exactly known motion) it marks both
moving images valid with speed errors of 0.055 and 0.057 % and the static
image not valid, so any larger error on real images comes from the
acquisition or the reference. Real multi-flash recordings are being prepared
and are not in the repository yet.
