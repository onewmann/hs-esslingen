# Acquisition

What a camera set-up needs so that the mean of a burst of frames, or one
multi-flash image, satisfies the model of [method.md](method.md), and how to
record and evaluate your own measurements. The first part covers frame
bursts, the main mode of this repository; the multi-flash set-up of the thesis
follows as an extension. The first real recording that the estimator
measures is a hand-held phone video at 120 frames per second
([phone video](../data/measurements/2026-10-16_phone_120fps/README.md)); a
recording against a physical reference speed is still to come, and no
multi-flash image has been recorded.
The numbers on this page come from the simulations and the independent
forward model in [results.md](results.md), from the 2025 camera recordings,
and from worked examples computed here.

## Frame bursts

The camera takes K consecutive frames at a fixed frame rate under continuous
light while the surface moves, and the mean of these frames is the image the
estimator analyses. If the exposures are short and the steps between the
frames are equal, it is the image that K flashes in one exposure would give:
the frame interval takes the place of the flash period, and the exposure
divided by the frame interval the place of the duty. Both live scripts record
such bursts (`--mode average --save DIR` in Python,
`'Mode', 'average', 'SaveDir', DIR` in MATLAB), and
`experiments/real_measurements.py` evaluates them with mode `average` in
`meta.csv`: it averages the frames, takes the frame period from `meta.csv` or
from the files (frame rate or timestamps), and checks the steps between the
frames by phase correlation and the frame intervals by their timestamps.
[`data/measurements/README.md`](../data/measurements/README.md#averaged-frames-mode-average)
gives what the recording writes, the columns of `meta.csv` in average mode and
what the step check means; `data/measurements/synthetic_average` is a
synthetic example.

```bash
python python/examples/live_camera.py --mode average --source pylon --exposures 9 --frame-rate 80 --exposure-us 2000 --save data/measurements/2026-10-20_average
python experiments/real_measurements.py data/measurements/synthetic_average/meta.csv --out results.csv
```

### Timing of a burst

The camera exposes each frame for a time $e$, and the frames start a frame
interval $\Delta t$ apart (1 / frame rate):

```
frames     |##|          |##|          |##|          ...          |##|
            e  <--- dt --->                                   (K frames)
burst      |<----------------- (K-1)*dt + e ---------------------->|
```

Between two frames the surface moves by $`v_\text{px}`$ pixels, the speed the
estimator measures (px per frame). The duty $d = e/\Delta t$ is the fraction
of the frame interval during which the camera exposes; during each exposure
the surface moves by $`d\,v_\text{px}`$ px, which blurs every copy along the
motion. A burst of K frames spans $(K-1)\,\Delta t + e$, and the surface has to
move at constant speed during it.

A burst keeps its frames, so the steps can be checked, which a multi-flash
image does not allow: the frame intervals from the timestamps, and the steps
themselves by phase correlation of consecutive frames. A dropped frame shows
as one interval and one step twice as long.

### Camera settings

* A fixed frame rate: AcquisitionFrameRateEnable on and AcquisitionFrameRate
  set (`--frame-rate` in `live_camera.py`, `'FrameRate'` in `live_basler`).
  Both scripts print the frame rate the camera reaches (ResultingFrameRate);
  if it is lower than the one set, the exposure, the ROI or the bandwidth
  limits it.
* ExposureAuto and GainAuto off, so that all frames have the same exposure
  and gain. `--exposure-us` and `'ExposureUs'` set the exposure in µs and
  switch ExposureAuto off; set the gain and GainAuto off in the pylon Viewer
  and save them in the user set the camera starts with.
* An exposure of at most a quarter of the frame interval (duty at most
  0.25). During one exposure the surface moves duty times the step, which
  blurs every frame along the motion. Give the exposure to the estimator
  (`exposure_s` in `meta.csv`, or `duty`); in simulation with the forward
  model, exposures of a quarter and half of the interval without it made
  12.5 % and 18.8 % of the images wrong but valid, and with the true duty all
  images were correct and valid. `real_measurements.py` warns above 0.25, and
  in particular between 0.4 and 0.6, where the speed fit depends strongly on
  the exact exposure (the live scripts warn between 0.4 and 0.6 and above
  0.9).
* Continuous light without flicker, for example LEDs on a DC supply. Lamps
  on the 50 Hz mains flicker at 100 Hz, so short exposures at a fixed frame
  rate get different brightness, and the frames no longer count equally in
  the mean.
* A camera ROI, so that the camera reaches the frame rate: the maximum frame
  rate grows when fewer rows are read out and less data is sent per frame.
  Set Width and Height (for example 512 x 512, the analysed square) and
  centre it with OffsetX and OffsetY in the pylon Viewer. `live_basler` sets
  a centred square of `'ROI'` px itself.
* Enough buffers and no dropped frames. `live_camera.py` grabs every burst
  of K frames into K + 2 buffers (at least 10), reports grabs that failed
  and skips bursts that lost frames; on GigE enable jumbo frames. In the
  recordings, the timestamps and the step check show a dropped frame.

### Choosing the frame rate

The motion per frame in pixels is

$$v_\text{px} = \frac{v\,\Delta t}{\mathrm{GSD}}$$

with $v$ the speed of the surface (mm/s), $\Delta t$ = 1 / frame rate the
frame interval (s) and GSD the ground sample distance, the length on the
surface that one pixel covers (mm/px,
[measured as below](#ground-sample-distance-and-physical-units)). The
estimator measures from about $3.35/K$ px per frame (0.37 px for K = 9,
0.08 px for K = 41). Its speed grid reaches $N/2$ for an $N$ px square, but
the studies only test up to $N/4$, so treat $N/4$ as the upper end (128 px
for a 512 px square). Within that range, more pixels per frame help when
there is noise: in the noise study the share of correct estimates rose from
28.0 % at 1 to 2 px per step to 82.0 % at 4 to 8 px and 98.1 % at 16 to
32 px, in simulation. Below about 2.2 px per step a single blurred exposure
and the pattern of K exposures cannot be told apart (blur check in
[method.md](method.md#blur-check-optional)). A few pixels per frame, about 4
to 10, is a good target, so $`\Delta t = v_\text{px}\,\mathrm{GSD}/v`$.

Worked example: GSD = 0.05 mm/px, the axis at 20 mm/s, 5 px per frame:
$\Delta t$ = 5 x 0.05 / 20 s = 12.5 ms, a frame rate of 80 Hz. The exposure
must be at most $\Delta t/4$ = 3.1 ms; 2 ms gives duty 0.16 and a blur of
0.16 x 5 = 0.8 px per frame. With K = 9 a window spans 8 frame intervals,
100 ms, and the surface moves 40 px, a small part of the 512 px square. For
100 mm/s the same 5 px per frame need $\Delta t$ = 2.5 ms, 400 Hz, and an
exposure of at most 0.6 ms; that frame rate usually needs a small camera ROI.

### Test surface

The surface needs a texture that does not repeat. A texture that repeats every
160 px or less in the analysed square gave wrong results marked valid in 50 to
100 % of the simulated images, and a static edge, mask or glint in the square
can do the same (assumptions 1 and 9 in
[method.md](method.md#assumptions)). Suitable are a matte print of random noise
(for example Gaussian noise made with numpy, smoothed so that one grain spans
a few pixels in the image and stretched to the full grey range), matte paper
or cardboard with a fine irregular texture, or sandpaper. Not suitable are woven fabric, a
grid, a checkerboard or any other regular pattern, and glossy surfaces with
reflections. Keep sheet edges, clamps, tape, screw heads and dust spots out of
the analysed square during the whole burst, and make the sheet longer than the
move plus the field of view. Static bursts of the surface must come out not
valid ([static images](#static-images)). Graph paper shows what a regular pattern
does: in a phone video it pulled the direction towards the grid lines by up to
5.2 deg, and 3 of 89 windows of the phone held still came out valid, 13 when
decoded with another OpenCV version, because their quality lies at the
threshold ([graph paper](../data/measurements/2026-10-17_phone_grid/README.md)).
Use such paper only to measure the scale.

### Phone videos

A phone can record usable bursts. In the
[phone video](../data/measurements/2026-10-16_phone_120fps/README.md), slow
motion at 120 frames per second gave about 5 to 14 px per frame in the windows
with even steps, for a hand moving at an ordinary pace, and all those windows
were valid. Three earlier videos at 30 frames per second (not in the
repository) moved by 76 to 93 px per frame with strong motion blur, and few
windows were valid ([results](results.md#phone-video-at-120-frames-per-second)).

* Record in slow motion (120 or 240 frames per second), or in a manual mode
  with a short exposure (1/1000 s) and the lamp on.
* Move the phone flat, at a constant height and a steady pace, for example
  along a ruler or the edge of the table.
* The lamp makes a bright spot that stays at the same place in the image.
  Analyse a square away from it (`extract_frames.py` of the phone session
  shows how to cut one from a video).
* The phone stores a timestamp for every frame; dropped frames show up as
  longer gaps, and the step check flags the windows that hold them.
* Video files are compressed with loss. It did not prevent the measurement
  here, but keep the bit rate high.

### Reference speed

To check an estimate, the surface has to move at a known, constant speed and
direction during the burst. Possible set-ups:

* a linear axis driven at a set speed, with its encoder giving the actual
  speed;
* the bed or the head of a 3D printer, moved by a G-code command at a set
  feed rate (for example `G1 Y150 F1200`: 1200 mm/min, 20 mm/s). The printer
  reaches the feed rate only between its acceleration ramps, so record in the
  middle of a long move;
* a conveyor belt, with its speed measured by an encoder on a roller or by
  timing a mark over a known distance.

Fasten the textured surface flat to the moving part, at the working distance
of the camera. The set speed goes into `set_speed_mm_s` and the direction in
the image into `set_angle_deg`
([angle convention](../data/measurements/README.md#angle-convention)). The
step check of `real_measurements.py` gives a second reference from the frames
themselves, the mean step by phase correlation, which depends neither on the
axis nor on the GSD. It also flags windows whose steps spread by more than
10 %, for example while the axis still accelerates. The uncertainty of the
reference speed limits the speed error in mm/s that can be claimed
([reporting](../data/measurements/README.md#reporting-the-result)).

### Ground sample distance and physical units

The estimator measures pixels per step. The speed in physical units is

$$v = \frac{v_\text{px}\,\mathrm{GSD}}{T}$$

with $T$ the frame interval for averaged frames, or the flash period for a
multi-flash image. `fftvel.to_physical(est, gsd, period)` does this in Python
and MATLAB; the command line takes `--gsd` and `--period` and then adds
`speed_per_s`, `vx_per_s` and `vy_per_s` to its JSON output, in the length
unit of the GSD per second:

```bash
python -m fftvel estimate data/measurements/synthetic_average/moving_regular.npz --stack --exposures 9 --duty 0.1 --gsd 0.05 --period 0.005
```

```python
from fftvel import estimate_velocity, pulsed_from_stack, to_physical
bp, K = pulsed_from_stack(frames, exposures=9)     # frames: rows x columns x frames, e.g. from fftvel.io.load
est = estimate_velocity(bp, exposures=K, duty=0.1)
p = to_physical(est, gsd=0.05, period=0.005)       # mm/s for the GSD in mm/px and 5 ms between frames
```

`to_physical` converts every estimate, valid or not, so check `est.valid`
first.

The relative error of the GSD passes one to one into every speed in mm/s.
Measure it with a ruler or a checkerboard:

1. Put the checkerboard or the ruler in the plane of the moving surface, at
   the same height, with the camera, lens, focus, zoom and working distance of
   the measurements. Refocusing or moving the camera changes the GSD.
2. Take an image and find the pixel coordinates of two marks far apart, for
   example 80 mm on the ruler or 10 squares of the checkerboard, near the
   centre of the image where the analysed square lies (lens distortion
   changes the GSD towards the corners). In MATLAB:
   `imshow(I); [x, y] = ginput(2); d = hypot(diff(x), diff(y))`.
3. GSD = known length / pixel distance. Two marks, each read to within 1 px,
   give the GSD to 0.09 % when they are 1600 px apart and to 0.88 % at
   160 px.
4. Measure along x and along y. If the two differ by more than about 0.5 %,
   the camera looks at the surface at an angle, and one GSD does not fit both
   directions.

The axis that moves the surface can measure it too: take an image of the
stopped surface, move the axis by a known distance D (say 10 mm), take
another, and get the displacement d in px from `fftvel.phase_correlation`
([angle convention](../data/measurements/README.md#angle-convention) shows the
call). Then GSD = D / |d|. Keep |d| below a quarter of the analysed square
(with a 1216 px square, below about 300 px).

### Static images

Record a number of bursts of the stopped surface with the same settings (K,
frame rate, exposure, ROI) and give them `static` 1 in `meta.csv`. All of them
must come out not valid. The step check warns about a static burst whose
frames move by more than 0.25 px. Move the surface to a new patch (or change
the set-up) between static bursts: repeated recordings of one patch differ
only by sensor noise, so they count as one in the bound below. `static` 1
means "must come out not valid", so it also fits a second control for
multi-flash images: images of the moving surface with a single flash in the
same exposure time (listed with the `exposures` of the real measurements).
They hold no flash pattern, and a valid one shows blur from other light or
vibration. With none of N independent static images valid, the one-sided
95 % Clopper-Pearson upper bound of the rate of static images marked valid is
$1 - 0.05^{1/N}$: 26 % for N = 10, 9.5 % for N = 30, 4.9 % for N = 60, 3.0 %
for N = 100. N is the number of distinct patches or set-ups, not of frames.
The summary of `real_measurements.py` prints this bound and counts every
static image it is given as independent.

### First recording, step by step

A plan for a first session with frame bursts. Fix the exclusion rule and the
success criteria below before looking at any estimate.

1. Camera: Mono8, ExposureAuto and GainAuto off, a fixed gain, DC LED light,
   focus at the working distance ([camera settings](#camera-settings)).
2. GSD: on a full-sensor image with the same lens, focus and working distance,
   measure it along x and y with marks far apart
   ([GSD](#ground-sample-distance-and-physical-units)). Then set a centred
   camera ROI of 512 x 512 px. Inside that ROI, marks are at most about 500 px
   apart, which gives the GSD only to about 0.3 %.
3. Direction: measure `set_angle_deg` once by phase correlation of two images
   of the stopped surface before and after a move
   ([angle convention](../data/measurements/README.md#angle-convention)).
4. Speeds: two or three set speeds that give about 4, 6 and 10 px per frame,
   with $`\Delta t = v_\text{px}\,\mathrm{GSD}/v`$ and an exposure of at most
   $\Delta t/4$ ([frame rate](#choosing-the-frame-rate)). Check that the frame
   rate the camera reports equals the set one. On a 3D printer, check that its
   maximum feed rate allows the speed.
5. Dry run: 3 moving and 3 static bursts, then `real_measurements.py`. Expect
   the frame period from the timestamps to equal 1 / frame rate, no uneven
   steps and no saturated pixels.
6. Moving bursts: at least 20 per speed, in the middle of long moves, with one
   note per speed (`--save DIR --note "20 mm/s"`).
7. Static bursts: at least 30, each at a different patch, with the surface
   moved by at least one field of view in between (`--save DIR --count 1
   --static --note "patch 7"`; in MATLAB, Stop after one burst). With none of
   30 valid, the bound above is 9.5 %.
8. Fill in `gsd_mm_per_px`, `set_speed_mm_s` and `set_angle_deg` in `meta.csv`
   and run `real_measurements.py` with `--figure`. Then comment out the rows
   with `uneven_steps` 1 in `meta.csv` with `#`, give the reason in `note`, and
   run it again.

Exclusion rule: the main numbers, from the second run, leave out windows with
`uneven_steps` 1 (a dropped frame or a change of speed, which does not depend
on the estimate). The first run, with all windows, is reported as well.
Windows that the summary names as not moving at the set speed stay in.

Success criteria, against two references:

* The estimator against the mean step between frames, which depends neither
  on the axis nor on the GSD: median, 95th percentile and maximum of
  $|v_\text{est}/v_\text{step} - 1|$ and of the angle difference, and the
  number of valid windows off by more than 5 % or 2 deg. The summary prints
  these.
* The whole chain against the set speed: at least 90 % of the moving windows
  within the measurable range valid, and no valid window off by more than 5 %
  or 2 deg. Windows that the summary names as not moving at the set speed
  (column `step_vs_set`) point to the axis, not to the estimator.
* None of the static bursts valid.

If the GSD comes from moving the same axis by a known distance, a scale error
of the axis cancels in the speed in mm/s: the comparison with the set speed
then tests the estimator, the timing and the GSD reading, but not the axis.

## Focus, saturation and error budget

### Focus

A defocused image has ring-shaped zeros in its spectrum, which look like slow
motion. The isotropy check is there to reject that pattern when nothing moves,
and both test suites check it on a defocused static texture. With motion,
defocus mostly costs valid results. In simulation with the forward model and
without noise, a disc blur of up to 8 px radius gave no wrong results marked
valid, and 31 or 32 of the 32 images per level were correct and valid; with the
noise of a 32000 e- sensor, the share of correct and valid images fell from
84.4 % in focus to 40.6 % at 2 px radius and 0 % at 5 px. Focus on the surface
at the working distance of the measurement, and check that a recording of the
stopped surface comes out not valid.

### Saturation and glints

Static bright structures do not move with the surface and give fringes of their
own. In simulation with the forward model, two static glints with a peak of the
mean signal made 28.1 % of the images wrong but valid, and 96.9 % from three
times the mean on. `clip_outliers` (`--clip-outliers`) removed the effect
there: all images were correct and valid. A sharp static edge of 1.5 times the
mean signal gave 12.5 % wrong but valid results. Of the 12 valid single-frame
crops of the room scene `vid1`, 7 lie in crops with at least 10 %
saturated pixels.

Set the exposure, the light or the flash energy so that no pixel saturates;
the estimate reports the share of saturated pixels in `saturated`, and the
command line and the scripts warn above 1 %. Shiny surfaces need diffuse
light or polarisers. Do not expose too dark either: 8-bit quantisation caused
no failures in simulation with the forward model, 6 bits gave 6.2 % wrong but
valid results, and 8 bits at a quarter of full scale are about 6 bits. Store
images in a lossless format (PNG, TIFF, BMP, `.npy`, `.npz`, `.mat`); JPEG
adds an 8 x 8 px block grid.

### Error budget

Each factor of $`v_\text{px}\,\mathrm{GSD}/T`$ adds its relative error to the
speed in mm/s; independent errors add in quadrature.

| Source | How it enters | Size |
|---|---|---|
| Estimator | relative error of $`v_\text{px}`$ | without noise, median 0.018 % and largest 0.47 % from 0.4 to 64 px per step (K = 9, 256 px), in simulation; noise and model violations in [results.md](results.md) |
| GSD | one to one | 0.09 % with marks 1600 px apart read to 1 px, 0.88 % with marks 160 px apart |
| Frame interval or flash period, systematic | one to one | frames: the camera clock, against which the timestamps check the intervals; flashes: measure the period with a photodiode and an oscilloscope |
| Timing jitter, per frame or flash | random error of $`v_\text{px}`$ | median speed error of the valid estimates 0.15, 0.35, 0.62, 1.38 and 3.16 % for a jitter of 1, 2, 5, 10 and 20 % of the period; wrong but valid results from 10 % on (forward model, in simulation) |
| Working distance with an entocentric lens | one to one through the GSD | 1.06 % for a height change of 3 mm at 300 mm distance with a 16 mm lens |
| K | factor $`K_\text{true}/K`$ | frames: the number averaged, known from the stack; flashes: for slow motion, 12.5 % if 9 flashes are taken for 8, not detected |

An entocentric lens is an ordinary lens: its magnification changes with the
object distance. For a thin lens of focal length $f$ at object distance $Z$,
$\mathrm{GSD} = p\,(Z-f)/f$ with $p$ the pixel pitch, so a change $\Delta Z$
changes the GSD by $\Delta Z/(Z-f)$. On a vehicle the body moves on its
suspension, and the height has to be measured or kept constant. A telecentric
lens keeps the magnification constant within its depth of field and removes
this term.

## Multi-flash images

One image with K flashes is the setting of the thesis. It measures motion
against the flash period, which a controller can set in microseconds, so it
works where separate frames share no content ([flash period](#choosing-the-flash-period)).
It needs a flash controller and a dark set-up, and nothing in the image shows
a wrong flash count or uneven flash intervals. Focus, saturation, the GSD and
the static controls are handled as for frame bursts.

### Timing of one image

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

### Choosing the flash period

The motion per flash in pixels is

$$v_\text{px} = \frac{v\,T}{\mathrm{GSD}}$$

with $v$ the speed of the surface (mm/s), $T$ the flash period (s) and GSD the
ground sample distance (mm/px). The measurable range is the same as for
frames, from about $3.35/K$ px per flash up to about $N/4$, and a few pixels
per flash, about 4 to 10, is a good target ([frame rate](#choosing-the-frame-rate)),
so $`T = v_\text{px}\,\mathrm{GSD}/v`$.

Worked example. A surface moves at 500 mm/s under a camera with
GSD = 0.05 mm/px (the 1936 px wide sensor then covers 96.8 mm). For 5 px per
flash, $T = 5 \cdot 0.05 / 500$ s = 0.5 ms, a flash rate of 2 kHz. With K = 9
the exposure must cover $(K-1+d)\,T$; with 20 µs flashes, $d$ = 0.04 and the
exposure needs at least 8.04 x 0.5 ms = 4.02 ms. Set, for example, 4.5 ms and
make sure that exactly 9 flashes fall into it. The blur per flash is
$`d\,v_\text{px}`$ = 0.2 px, negligible. With the same T, speeds from 100 to
1000 mm/s give 1 to 10 px per flash; the lower end of 0.37 px corresponds to
37 mm/s. During the 9 flashes the surface moves 40 px, a small part of a
512 px square.

The same arithmetic shows where a multi-flash image is needed. At 30 m/s with
a GSD of 0.2 mm/px, a 256 px field is 51.2 mm wide. At 100 frames per second
two consecutive frames are 300 mm apart and do not overlap, while flashes
133 µs apart give 20 px per flash inside one exposure.

### Flash duration and duty

Keep the flashes short, so that $`d\,v_\text{px}`$ stays well below about 3 px,
or give the duty to the estimator (`duty`, `--duty`). In simulation with the
forward model, flashes of a quarter and half of the period without the option
made 12.5 % and 18.8 % of the images wrong but valid, with the speed near
$`d\,v_\text{px}/K`$. With `duty` set to the true value, every image from duty
0 to 0.9 was correct and valid. Near duty 0.5, for example with a 50 % square
wave, every second main line of the pattern vanishes, so there the value must
be exact.

### The exact flash count

The estimator needs the true number of flashes K: `exposures` (`--exposures`,
`'Exposures'`) for any K of at least 2, also even, or `L` for K = 2L+1, the
half-count used in the thesis. For slow motion a wrong K scales the speed by
$`K_\text{true}/K`$ without lowering the fit quality. In simulation with the
forward model, telling the estimator 7, 8, 10 or 11 flashes for an image with
9 made 34.4 % of the images wrong but valid, with the speed off by the factor
9/K.

The flashes must also be equally spaced. Intervals alternating between 0.8 and
1.2 periods gave 87 of 240 images valid at twice the true speed; one missing
flash gave 22 of 240 valid but wrong, in simulation. Random jitter of a clocked
trigger is less critical: in simulation with the forward model, a standard
deviation of 1 % of the period gave a median speed error of 0.15 %, and wrong
but valid results appeared only from 10 %. Check the number of pulses and the
period once with a photodiode and an oscilloscope. With the flash period given,
both live scripts and `experiments/real_measurements.py` warn when the exposure
is shorter than $(K-1+d)\,T$.

### Ambient light

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
Record controls that must come out not valid
([static images](#static-images)): images of the stopped surface, and images
of the moving surface with a single flash in the same exposure time. A valid
control shows blur from other light or from vibration. `reject_blur` rejects
blur-like results, at the cost of many correct ones below about 2.2 px per
flash.

## Live scripts

Two scripts stream from a camera and estimate every burst or image. Both need
an explicit mode:

* `average`: continuous light; K consecutive frames are averaged into one
  image, and both scripts check that the frame intervals are equal. The 2025
  recordings show why this check matters
  ([history.md](history.md#what-the-review-found)). With a Basler camera both
  scripts read the exposure time and pass it on as the duty. The speed is in
  px per frame, and in px/s from the frame times.
* `flash`: every frame holds K strobe flashes and is estimated on its own; the
  speed is in px per flash, and in px/s with the flash period.

```bash
python python/examples/live_camera.py --mode average --source pylon --exposures 9 --exposure-us 500
python python/examples/live_camera.py --mode flash --source pylon --exposures 9 --flash-period 0.0005 --duty 0.04
python python/examples/live_camera.py --mode flash --source recording.avi --exposures 9
```

```matlab
addpath('matlab', 'matlab/examples')       % from the repository root
live_basler('Mode', 'average', 'Exposures', 9, 'ExposureUs', 500)
live_basler('Mode', 'flash', 'Exposures', 9, 'FlashPeriod', 5e-4, 'Duty', 0.04)
```

[`python/examples/live_camera.py`](../python/examples/live_camera.py) reads
Basler cameras through pypylon, and webcams and video files through OpenCV.
[`matlab/examples/live_basler.m`](../matlab/examples/live_basler.m) uses the
Image Acquisition Toolbox with the GenTL adaptor; `'ExposureUs'` is the
exposure time in µs and `'UpdateInterval'` the time between estimates in s.
Both default to a 512 px square, which takes about 4 s per estimate in Python
(repeated call, [results.md](results.md#run-time)), so they give a reading
every few seconds.

Test status: neither script has run on camera hardware. The Python script is
covered by tests with fake pypylon and OpenCV modules
(`tests/python/test_live_camera.py`), and its pylon path was checked with the
pylon camera emulator. `live_basler.m` has no automated test, because GNU
Octave has no Image Acquisition Toolbox, and it has not run in MATLAB.

## Evaluating your own measurements

[`experiments/real_measurements.py`](../experiments/real_measurements.py)
estimates every burst or image listed in a `meta.csv`, converts the results
to mm/s and compares them with the set motion.
[`data/measurements/README.md`](../data/measurements/README.md) describes the
folder layout, the columns of `meta.csv` (file, mode, number of frames or
flashes, frame interval or flash period, exposure or duty, GSD, set speed and
direction, static flag, optional single-flash reference frames), the angle
convention, the columns of `results.csv` and how to report the result. Try it
on the synthetic examples first:

```bash
python experiments/real_measurements.py data/measurements/synthetic_average/meta.csv --out results.csv
python experiments/real_measurements.py data/measurements/synthetic_flash/meta.csv --out results.csv
```

On the burst example (simulated frames with exactly known motion) it marks
both moving windows valid with speed errors of 0.08 and 0.007 %, flags the
window with a dropped frame through the step check and marks the static
window not valid. On the multi-flash example it marks both moving images
valid with speed errors of 0.055 and 0.057 % and the static image not valid.
These errors come from the estimator alone, without the errors of a real
axis, GSD or timing. The first real recording, a hand-held phone video at 120
frames per second, is in `data/measurements/2026-10-16_phone_120fps`; a
recording against a physical reference speed is still to come.

### Checklist

For frame bursts:

* A fixed frame rate, ExposureAuto and GainAuto off, enough buffers and no
  dropped frames ([camera settings](#camera-settings)). The frame-interval
  check and the step check flag windows with uneven steps.
* An exposure of at most a quarter of the frame interval, given in
  `exposure_s` (or as `duty`).
* Continuous light without flicker, for example LEDs on a DC supply.
* A constant speed during each burst and a known reference speed
  ([reference speed](#reference-speed)).
* A random, non-repeating texture with no static edges in the analysed square
  ([test surface](#test-surface)).
* The number of averaged frames in the column `frames`.

For multi-flash images:

* Exactly K flashes in every image, with the camera starting the burst or the
  controller triggering the camera ([timing](#timing-of-one-image)). K goes
  into the `exposures` column; a wrong K scales the speed by $`K_\text{true}/K`$
  for slow motion, and the result can still be marked valid
  ([flash count](#the-exact-flash-count)).
* Equal flash intervals, with the number of pulses and the flash period
  checked once with a photodiode and an oscilloscope. An error in the period
  passes one to one into the speed in mm/s.
* An exposure of at least $(K-1+d)\,T$ that holds exactly K flashes, and no
  other light during it: an image with the flashes switched off is black
  ([ambient light](#ambient-light)).

For both:

* No saturated pixels. The script warns when more than 1 % of the analysed
  pixels are saturated (column `saturated`); `--clip-outliers` helps against
  small glints ([saturation](#saturation-and-glints)).
* The surface in focus at the working distance ([focus](#focus)).
* The GSD measured at the working distance
  ([GSD](#ground-sample-distance-and-physical-units)) and the set direction
  measured once
  ([angle convention](../data/measurements/README.md#angle-convention)).
* Static controls recorded with the same settings
  ([static images](#static-images)).
* Lossless files: PNG, TIFF, BMP, `.npy`, `.npz` or `.mat`.
