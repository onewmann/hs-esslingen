# Method

This page explains how `fftvel` gets direction and speed from one image that
holds $K$ equally spaced copies of a moving surface, the mean of $K$ camera
frames or one exposure with $K$ flashes, step by step, and where each step
comes from. Two tags mark the origin:

* [thesis]: taken from the doctoral thesis of Jan Horn (full reference at the
  end of this page).
* [rework]: added in the rework of the 2025 student project, that is after the
  tag `student-project-2025` (see [history.md](history.md)).

The code is [`python/src/fftvel/core.py`](../python/src/fftvel/core.py) and
[`matlab/+fftvel/estimate_velocity.m`](../matlab/+fftvel/estimate_velocity.m),
which compute the same numbers (checked in GNU Octave and in MATLAB). What the method achieves is in
[results.md](results.md), what an acquisition needs in
[acquisition.md](acquisition.md).

## The image model

The thesis and the function names of `fftvel` call the analysed image a
pulsed-exposure image (gepulste Belichtung, as the 2025 scripts write it), and
the option `exposures` counts its $K$ exposures: the averaged frames, or the
flashes. $L$ is the half-count of the thesis, with $K = 2L+1$.

Two acquisitions give the same image. With averaged frames, a camera takes
$K$ frames at a fixed frame rate while the textured surface moves, and their
mean is computed afterwards; the step $\mathbf v$ is the motion in pixels
between two frames. With flashes, a strobe fires $K$ times during one camera
exposure and the sensor adds up the copies; the step is the motion between
two flashes. The second is the setting of the thesis. In both cases the $K$
exposures follow each other at equal steps.

[thesis] The image is the mean of the $K = 2L+1$ shifted copies of the
texture $b$:

$$b_p(\mathbf x) = \frac{1}{2L+1}\sum_{k=-L}^{L} b(\mathbf x - k\mathbf v)$$

By the shift theorem of the Fourier transform each shift becomes a phase
factor, and the sum of these factors is a Dirichlet kernel $D$:

$$B_p(\mathbf f) = B(\mathbf f)\,D(\mathbf f\cdot\mathbf v),\qquad D(u) = \frac{\sin\big((2L+1)\pi u\big)}{(2L+1)\sin(\pi u)}$$

Here $\mathbf f$ is the spatial frequency in cycles per pixel and $B$ the
spectrum of the texture. $D$ equals 1 wherever $\mathbf f\cdot\mathbf v$ is an
integer, has $2L$ zeros and $2L-1$ weak side lobes between two such places,
and is small elsewhere. The magnitude spectrum of the texture is therefore
covered by straight lines perpendicular to $\mathbf v$, spaced $1/|\mathbf v|$
cycles per pixel apart. Their orientation gives the direction of motion, their
spacing the length of the step.

For an even number of exposures, or a sum from $k = 0$ to $K-1$, only a phase
factor changes. The magnitude is $|\sin(K\pi u)/(K\sin\pi u)|$ for any $K$,
which is why the estimator accepts any $K \ge 2$ (option `exposures`). The
magnitude is the same for $\mathbf v$ and $-\mathbf v$, so one image gives the
direction only modulo 180 deg. With averaged frames, phase correlation of two
consecutive frames gives the sign.

The figure shows every stage for a simulated image: the camera texture of
[`data/surface_texture.png`](../data/surface_texture.png), moving by (5, -3) px
per frame, as the mean of 9 frames with short exposures.

![Every stage of one estimate](figures/pipeline.png)

## Step 1: window and log spectrum

[rework] The analysed region is the central square of the image, of side $N$:
the largest square, or a smaller one with `roi` (the command line, the live
scripts and `experiments/real_measurements.py` default to 512 px). Non-square
images are cropped, never resized, so both axes keep the camera's pixel scale.
The live and video scripts of 2025 resized every frame to 256 x 256 px, which
scaled the two axes differently
([history.md](history.md#what-the-review-found)).

[rework] With `clip_outliers`, every pixel more than 6 robust standard
deviations from the median (1.4826 times the median absolute deviation), and
every pixel within 3 px of such a pixel, is set to the median first.

The mean is subtracted and the square is multiplied by a separable 2-D Hann
window, a raised-cosine taper that falls to 0 at the borders so that the edges
of the square add no lines of their own (the 2025 scripts used the same
window). The FFT is zero padded to $\tilde N = 2N$ points per axis (`npad` in
the code, `pad_factor` 2), which samples the spectrum twice as densely.

[thesis] The estimator works on the logarithm of the magnitude spectrum. The
logarithm compresses the large range of $`|B_p|`$ and turns the product
$B\cdot D$ into a sum, so the line pattern adds to the texture spectrum
instead of scaling it.

[rework] The floor inside the logarithm is relative to the image,
$`S = \log\big(|B_p| + 0.05\,\mathrm{median}\,|B_p|\big)`$, so scaling the image
shifts $S$ by a constant and changes nothing else.

## Step 2: direction

[thesis] A Radon transform integrates the log spectrum along parallel
projection paths, one projection for every angle $\theta$. Along the direction
of motion the projection crosses the spectral lines at right angles and
becomes a sharp comb; at other angles the lines are smeared over many bins.
The ripple profile WR (Welligkeitsprofil in the 2025 scripts) measures how
strongly each projection oscillates,

$$WR(\theta) = \sum_{j} \big|M(j+1,\theta) - M(j,\theta)\big|,$$

taken over the bins $|j| \le 0.9\,(\tilde N/2 - 1)$, and the direction of
motion is the angle where WR is largest.

[rework] $M(j,\theta)$ is the mean of $S$ along each projection path instead
of the sum. Every sample inside the disc of radius $\tilde N/2 - 1$ is split
linearly between its two nearest bins, and every bin is divided by its weight.
With sums over the square spectrum, as in the 2025 scripts, the projection
paths through the corners at 45 and 135 deg are longer, which biases the
maximum towards the diagonals. The angle grid runs from 0 to 179 deg in steps
of 1 deg, and a parabola through the maximum and its two neighbours (periodic
over 180 deg) refines it.

[thesis] GEBA $= \max WR / \mathrm{mean}\,WR$ is the quality measure of the
thesis: it says how clearly one direction stands out. The 2025 scripts expand
the name as Güte der Einzelbild-Analyse (quality of the single-image
analysis) in the header of `main.m`. `fftvel` reports it as `geba` and does
not use it for `valid`. In the threshold calibration the fit
quality of step 3, scaled by the profile length, separated correct estimates
from images without motion better than GEBA, in simulation: the area under the
ROC curve (the probability that a randomly chosen correct estimate scores
higher than a randomly chosen image without motion) was 0.988 against 0.947
([results.md](results.md#threshold-calibration)).

## Step 3: speed

[thesis] The thesis takes the speed from the spacing of the repeated
structures in the projection along the direction of motion: the main lines lie
$\tilde N/|\mathbf v|$ bins apart.

[rework] `fftvel` fits a model of the whole projection instead of locating
peaks. The profile $g(j)$ is the mean of the projection at $+j$ and $-j$, for
$j$ from 2 to $0.9\,(\tilde N/2 - 1)$; the two bins next to zero frequency
hold the central peak and are left out. For a random texture the expected
power of the windowed pulsed image at bin $j$ is

$$P(j) = \frac{1}{K^2}\Big[K + 2\sum_{s=1}^{K-1}(K-s)\,\rho(s\mathbf v)\,\cos\Big(\frac{2\pi j s |\mathbf v|}{\tilde N}\Big)\Big]$$

Each term pairs two exposures $s$ steps apart, and there are $K-s$ such pairs.
$\rho(s\mathbf v)$ is the autocorrelation of the window at the offset between
the two copies: copies shifted by $s\mathbf v$ overlap only where the window
covers both, so their correlation falls as the offset grows. For the
separable Hann window of $N$ points it has a closed form:

$$\rho(s\mathbf v) = \rho_w(s v_x)\,\rho_w(s v_y),\qquad \rho_w(\tau) = (1-x)\,\frac{2+\cos 2\pi x}{3} + \frac{\sin 2\pi x}{2\pi},\qquad x = \frac{|\tau|}{N-1}$$

for $x < 1$, and $`\rho_w = 0`$ beyond. For $N = 256$ it agrees with the
discrete autocorrelation of the window to $3.3\cdot 10^{-10}$. Since $P$ is a
finite sum of cosines, it contains the blur that the window causes in the
spectrum and cannot alias.

If each exposure lasts the fraction $d$ of the time between two exposures
(option `duty`), the texture moves by $d\,|\mathbf v|$ during it, and $P(j)$ is
multiplied by
$\mathrm{sinc}^2(j\,d\,|\mathbf v|/\tilde N)$ with
$\mathrm{sinc}(x) = \sin(\pi x)/(\pi x)$. For averaged frames the duty is the
exposure time over the frame interval, for flashes the flash duration over
the flash period.

The template for one speed is $\tfrac12\log\big(P(j) + 10^{-3}\big)$, on a
logarithmic grid of speeds with steps of 0.5 % (factor $e^{0.005}$). A cubic
polynomial in $j$ is removed from the profile and from every template; it
takes up the smooth spectrum of the texture itself. The fit quality is the
correlation coefficient of the detrended profile and the detrended template.
The speed with the highest quality wins, and a parabola on the log-speed grid
refines it. The grid starts at

$$v_{\min} = \frac{1.5\,\tilde N}{K\,j_{\max}} \approx \frac{3.35}{K}\ \text{px per step},\qquad j_{\max} = \big\lfloor 0.9\,(\tilde N/2 - 1)\big\rfloor,$$

where the first zero of the Dirichlet pattern still lies well inside the
profile. That is 0.37 px for $K = 9$ and 0.08 px for $K = 41$, for any square
from 128 to 1024 px. Below it the grid continues down to a quarter of that
value as a guard band. The top of the grid is $\tilde N/4$, where the main
lines are 4 bins apart. The studies in [results.md](results.md) test speeds up
to $N/4$ (64 px per step for a 256 px square).

## Step 4: validity

[rework] An estimate is valid when the fit quality reaches a threshold, the
speed lies inside the measurable range and the spectrum is not isotropic.
With the option `reject_blur` it must also not look like a single blurred
exposure. The result carries every part of this decision (`quality`,
`threshold`, `in_range`, `isotropic`, `blur_like`), so a front end can report
why an estimate was rejected.

### Quality threshold

Even an image without motion matches one of the many templates to some
degree, and a short profile matches better by chance than a long one. The
Fisher transform $\operatorname{atanh}(q)$ maps a correlation coefficient $q$
to a variable whose spread no longer depends on $q$ itself. If the profile has
$m$ independent points, the transform of the best chance correlation scales
with $1/\sqrt{m-5}$: the cubic takes four points and the fit one more. The
threshold is therefore

$$q_{\text{thr}} = \max\Big(0.6,\ \tanh\frac{9}{\sqrt{m-5}}\Big),\qquad m = \min\big(n_j,\ \lfloor 2N n_j/\tilde N\rfloor\big)$$

with $`n_j`$ the number of profile points. Zero padding beyond the factor 2 only
interpolates between points, so it does not raise $m$ and does not lower the
threshold.

The constant 9 was checked with
[`experiments/calibrate_threshold.py`](../experiments/calibrate_threshold.py).
Each image gives $k = \operatorname{atanh}(q)\,\sqrt{m-5}$, and the quality
alone passes exactly when $k \ge 9$. On 1008 simulated images without motion
(white noise, a static random texture and static crops of the camera texture
with a little noise; 32 to 256 px; K = 3, 5, 9 and 21, that is L = 1, 2, 4
and 10; 12 seeds each) the
largest $k$ was 8.41 and the 99th percentile 7.47. On all 2268 simulated
images without motion, held-out textures included, the largest $k$ was 8.64,
and none was valid. None of these textures is strongly anisotropic; a static
texture of streaks passes in about half of the images (assumption 1). These
sets also served to check the constant, so the
held-out real images are the better guide to real false alarms: of 252 crops
of single frames of the room scene `vid1`, 4 were valid. These crops come from
the 10 frames of one scene, so they are not independent trials and say little
about other scenes. The cost of the threshold is that some correct estimates on small images are rejected: with 9, 52.7 % of
the correct estimates at 32 px are valid and 98.5 % at 256 px, in simulation.

| Square (px) | m | Threshold | Largest k without motion | Largest quality without motion |
|---:|---:|---:|---:|---:|
| 32 | 26 | 0.961 | 8.30 | 0.948 |
| 48 | 41 | 0.905 | 7.97 | 0.869 |
| 64 | 55 | 0.855 | 7.93 | 0.808 |
| 96 | 84 | 0.767 | 8.00 | 0.717 |
| 128 | 113 | 0.699 | 8.41 | 0.669 |
| 192 | 170 | 0.605 | 7.97 | 0.551 |
| 256 | 228 | 0.600 | 7.74 | 0.476 |

*Threshold per square size and the largest values among the 1008 calibration
images, from `experiments/results/threshold_calibration.csv`. At 256 px the
floor of 0.6 decides; at 192 px the two are almost equal.*

### Measurable range

`in_range` is false when the best speed lies in the guard band below
$`v_{\min}`$ or at the top end of the grid. Motion slower than the range is
then reported as not valid instead of as a wrong speed: in simulation, all 24
images at 0.3 px per step with $K = 9$ came out not valid.

### Isotropy check

The same speed fit is repeated on the projection perpendicular to the
direction found. Motion makes the spectrum anisotropic, so that fit is
clearly worse or gives another speed. If its quality is less than 0.15 lower
and its speed lies within a factor $e^{0.12} \approx 1.13$, the structure is
isotropic and the estimate is not valid. It is meant to catch defocus, whose
disc blur puts ring-shaped zeros into the spectrum that look like slow motion
in every direction, and square repeating textures, which have the same period
along both axes; the test suite checks it on a defocused camera texture and on
a repeating texture. On real graph paper held still it flagged three of the
six windows that reached the quality threshold and let the other three
through as valid ([results.md](results.md#graph-paper-a-periodic-texture)). A
static texture that is itself strongly anisotropic, such as streaks, passes it
(assumption 1). In the noise study it flagged 10 of 840 images, all of which
also failed the quality threshold, in simulation.

### Blur check (optional)

A single exposure during which the surface moves by $d$ px (motion blur, one
exposure instead of $K$) multiplies the spectrum by
$\mathrm{sinc}^2(j d/\tilde N)$. Its main lobe looks like that of $K$
exposures $d/K$ apart. Only the repeated main lines of the $K$ exposures,
$\tilde N/|\mathbf v|$ bins apart, tell the two apart, and
below $`|\mathbf v| = \tilde N/j_{\max}`$ (2.24 px for a 256 px square) they lie
beyond the fitted profile. `fit_blur` therefore fits a template for one
continuous exposure of free length $d$, sampled every $h = 0.5$ px along the
direction $\hat{\mathbf e}$ found in step 2:

$$P_{\text{blur}}(j) = \frac{h}{d}\Big[1 + 2\sum_{q\ge1}\max\Big(1-\frac{qh}{d},\,0\Big)\,\rho(qh\,\hat{\mathbf e})\cos\frac{2\pi j q h}{\tilde N}\Big]$$

`blur_like` is true when this template fits better than the $K$ exposures and
the speed is below $`\tilde N/j_{\max}`$. The result reports `quality_blur`,
`blur_length` and `blur_like` in any case; only `reject_blur` makes
`blur_like` count against `valid`. On single camera frames of the recording
`vid4`, all 28 crops that came out valid were blur-like. In the noise study, 21
of the 47 correct valid estimates below 2.24 px were blur-like as well, in
simulation, so `reject_blur` would have rejected them.

## Saturation warning

[rework] `saturated` is the share of analysed pixels at the maximum of the
square, measured before any processing. For an integer image whose maximum is
the largest value of its type (255 for 8 bit), any share counts; otherwise the
share counts from 0.5 % on. It does not change `valid`. The command line, the
live scripts and `experiments/real_measurements.py` warn above 1 %, because
glints and saturated areas are static bright structures: in the simulated
forward model of
[results.md](results.md#model-violations-in-an-independent-forward-model), two
glints of 3 times the mean signal or more made 96.9 % of the images wrong but
valid.

## Assumptions

The estimator gives a correct result only if the image satisfies the model.
These are all the assumptions it makes; [results.md](results.md) shows what
happens when some of them fail, and [acquisition.md](acquisition.md) how an
acquisition meets them.

1. The surface carries a random, non-repeating texture over the whole analysed
   square. A texture that repeats every 160 px or less gave wrong results
   marked valid in 50 to 100 % of the images, in simulation with the forward
   model. A static texture of streaks, with a correlation length along the
   streaks 16 times or more that across them (1 to 4 px across, 16 to 64 px
   along), came out valid in 39 of 72 simulated 256 px images, as motion of
   4.4 to 4.7 px per step along the streaks, and `reject_blur` and
   `clip_outliers` do not change that; the isotropy check does not reject it,
   because the texture itself is anisotropic. The same textures in motion were
   measured correctly in all 32 cases, so the risk is a static brushed or
   grained surface reported as moving
   ([static textures of streaks](results.md#static-textures-of-streaks)).
2. The motion is a pure translation in the image plane during the burst: no
   rotation, no change of scale (constant working distance) and no change of
   perspective.
3. The velocity is constant during the $K$ frames or flashes.
4. The frames or flashes are equally spaced in time: a fixed frame rate
   without dropped frames, or flashes from a clocked trigger. Alternating
   intervals and a missing frame or flash both gave wrong results marked
   valid in simulation.
5. The number $K$ of averaged frames or of flashes is known exactly.
6. Each exposure is short compared with the time between two exposures, or
   the duty is known (`duty`): the exposure time over the frame interval, or
   the flash duration over the flash period.
7. All frames or flashes have about the same energy, so the light must not
   flicker from frame to frame.
8. For a multi-flash image, no light other than the flashes reaches the
   sensor during the exposure. Continuous light adds a motion-blurred image of
   the surface. Averaged frames are taken under continuous light; there the
   exposure of each frame is the duty of assumption 6.
9. Nothing in the analysed square stands still: no sharp static edges, masks,
   glints or hot pixels.
10. The camera is in focus.
11. The sensor responds linearly and does not saturate, and the image is
    stored without lossy compression (JPEG adds an 8 x 8 px block grid).
12. The pixels are square and the image is not rescaled.
13. The speed lies between about $3.35/K$ and $N/4$ px per step.
14. The analysed square has at least about 64 px; below 128 px the stricter
    threshold rejects many correct estimates.
15. The texture spectrum is smooth enough for a cubic polynomial in $j$ to
    absorb it along the direction of motion.

The sign of the motion is never known from the image: it gives the direction
only modulo 180 deg. With averaged frames, phase correlation of two
consecutive frames gives the sign.

## Options and when to use them

The options are keyword arguments in Python
(`estimate_velocity(bp, exposures=9)`) and name-value pairs in MATLAB
(`'Exposures', 9`).

| Python | MATLAB | Default | Use it when |
|---|---|---|---|
| `exposures` | `'Exposures'` | 2L+1 | always: the image holds K averaged frames or K flashes, also an even number. Give the true count: a wrong K scales slow speeds by K_true/K and can stay valid. In simulation with the forward model, telling the estimator 7, 8, 10 or 11 instead of 9 made 34.4 % of the images wrong but valid |
| `L` | `'L'` | 4 | the half-count of the thesis: the image holds 2L+1 frames or flashes. `exposures` overrides it |
| `duty` | `'Duty'` | 0 | the motion during one exposure (duty times speed) is not well below about 3 px. The duty is the exposure time over the frame interval for averaged frames, the flash duration over the flash period for flashes. In simulation with the forward model without the option, duty 0.25 and 0.5 gave 12.5 % and 18.8 % wrong but valid results, with the speed near duty times speed divided by K; with `duty` set to the true value all images from 0 to 0.9 were correct and valid. Near 0.5 every second main line vanishes, so the value must be known |
| `roi` | `'ROI'` | largest square; 512 px in the command line and the live scripts | side of the analysed central square. Set it when the image is larger than needed. Each doubling of the side multiplies the run time by 4 to 6.5: 0.61 s for 256 px and 3.96 s for 512 px per repeated call in Python. 256 to 512 px is a good range |
| `clip_outliers` | `'ClipOutliers'` | false | static bright spots (glints of the light, hot pixels) may be in the image. In simulation with the forward model two glints of 3 times the mean signal made 96.9 % of the images wrong but valid; with the option all were correct and valid. It is off by default because it also flattens bright details of the texture |
| `reject_blur` | `'RejectBlur'` | false | a single blurred exposure may be estimated instead of K exposures, for example a multi-flash image with continuous light or vibration. Below about 2.2 px per step it also rejects many correct results (see the blur check) |
| `pad_factor` | `'PadFactor'` | 2 | rarely; larger values interpolate the spectrum and cost time |
| `min_quality` | `'MinQuality'` | 0.6 | a stricter floor than 0.6 is wanted; small images need more anyway |
| `diagnostics` | `'Diagnostics'` | false | the intermediate results are needed, for example for `plot_estimate` |

## Constants

Every number of the algorithm is a named constant. The table is a copy of
[`tests/constants.json`](../tests/constants.json), which also lists where each
constant is used; tests in both languages check that `core.py` and
`matlab/+fftvel/+internal/constants.m` hold exactly these names and values,
and `tests/python/test_docs_constants.py` checks this table against the file
([CONTRIBUTING.md](../CONTRIBUTING.md#constants)).

<!-- constants table: copy of tests/constants.json, checked by tests/python/test_docs_constants.py -->
| Name | Value | Meaning |
|---|---|---|
| `DEFAULT_L` | 4 | Default L: half the number of exposures, K = 2L+1 = 9. |
| `DEFAULT_PAD_FACTOR` | 2 | Default zero padding factor of the FFT: npad = pad\_factor \* n. |
| `DEFAULT_MIN_QUALITY` | 0.6 | Default lowest fit quality accepted for valid (min\_quality / MinQuality). |
| `MIN_SIZE` | 16 | Smallest side of the analysed square, and smallest roi, in px. |
| `MIN_EXPOSURES` | 2 | Smallest number of exposures K accepted by the exposures option. |
| `LOG_FLOOR` | 0.05 | Floor of the log spectrum as a fraction of the median of \|B\|: S = log(\|B\| + LOG\_FLOOR\*median\|B\|); the mean of \|B\| replaces the median if the median is 0. Keeps S independent of the image scale. |
| `BLANK_FLOOR` | 1.0 | Floor used instead when \|B\| is 0 everywhere (blank image); the estimate then comes out not valid. |
| `CLIP_SIGMA` | 6 | clip\_outliers: pixels more than CLIP\_SIGMA robust standard deviations from the median are set to the median. |
| `MAD_SCALE` | 1.4826 | Robust standard deviation = MAD\_SCALE \* median absolute deviation (1.4826 for normal noise). |
| `CLIP_RADIUS` | 3 | clip\_outliers: the pixels within CLIP\_RADIUS px (dy^2 + dx^2 &lt;= CLIP\_RADIUS^2) of a clipped pixel are set to the median as well, so that the faint rim of a glint goes with its core. |
| `ANGLE_STEP` | 1.0 | Step of the angle grid in degrees: 0, ANGLE\_STEP, ..., 180 - ANGLE\_STEP (0..179). The maximum of WR is refined by a parabola through three neighbouring angles (periodic); the speed maximum likewise on the speed grid. |
| `RFRAC` | 0.9 | Fraction of the spectrum radius rmax = npad/2 - 1 used for the ripple profile WR and for the speed fit. |
| `WEIGHT_FLOOR` | 1e-12 | Lower bound of the bin weights in the mean projection, so empty bins do not divide by 0. |
| `PROFILE_START` | 2 | First bin of the half profile used for the speed fit (bins 0 and 1 hold the peak at zero frequency). |
| `SPEED_STEP` | 0.005 | Step of the logarithmic speed grid: v(k+1)/v(k) = exp(SPEED\_STEP). |
| `VMIN_FACTOR` | 1.5 | Lower end of the measurable range: vmin = VMIN\_FACTOR\*npad/(K\*max(jpos)), where the first spectral zero still lies well inside the profile. |
| `MIN_LINE_SPACING` | 4 | Upper end of the measurable range: vmax = npad/MIN\_LINE\_SPACING, where the motion lines are this many bins apart. |
| `GUARD_FACTOR` | 4 | The speed grid extends down to vmin/GUARD\_FACTOR as a guard band; a best speed there gives in\_range = false. |
| `TREND_DEGREE` | 3 | Degree of the polynomial in jpos that absorbs the smooth texture spectrum before the correlation (cubic). |
| `TEMPLATE_FLOOR` | 0.001 | The Dirichlet and box templates are 0.5\*log(P + TEMPLATE\_FLOOR). |
| `FLAT_MODEL` | 1e-06 | Speeds (or blur lengths) whose detrended template has a norm below FLAT\_MODEL times the norm of the template score 0 (flat model, no information). |
| `NULL_K` | 9 | Null quality tanh(NULL\_K/sqrt(m - NULL\_OFFSET)): the largest quality was tanh(8.41/sqrt(m-5)) on the 1008 calibration images without motion and tanh(8.64/sqrt(m-5)) on all 2268 simulated images without motion (experiments/calibrate\_threshold.py); none of them is valid with 9. |
| `NULL_OFFSET` | 5 | Points taken from the m independent profile points: TREND\_DEGREE + 1 by the polynomial and 1 by the fit; the threshold is infinite for m &lt;= NULL\_OFFSET. |
| `INDEPENDENT_PAD` | 2 | Independent profile points m = min(points, floor(INDEPENDENT\_PAD\*n\*points/npad)): padding beyond this factor only interpolates and does not lower the threshold. |
| `ISO_GAP` | 0.15 | Isotropy check: isotropic if quality - quality\_perp &lt; ISO\_GAP and ... |
| `ISO_LOG_SPEED` | 0.12 | ... \|log(speed\_perp/speed)\| &lt; ISO\_LOG\_SPEED (defocus blur or a square repeating texture, not motion). |
| `BLUR_STEP` | 0.5 | Blur check: the continuous blur of box\_template is sampled every BLUR\_STEP px (for a blur of K\*BLUR\_STEP px it equals the Dirichlet template of K exposures BLUR\_STEP apart). |
| `BLUR_GAP` | 0.0 | Blur check: blur\_like if quality - quality\_blur &lt; BLUR\_GAP (0: one continuous exposure, motion blur without flashes, fits better than the K flashes) and the speed is below npad/max(profile\_j), where the repeated lines of the flashes lie beyond the fitted profile. |
| `SAT_SHARE` | 0.005 | saturated: the share of pixels at the maximum counts only when it is at least SAT\_SHARE, unless the maximum is the largest value of an integer type (then any share counts). |
| `SAT_WARN` | 0.01 | Front ends (command line, stream, matlab/examples/live\_basler.m, experiments/real\_measurements.py) warn when saturated exceeds SAT\_WARN. |
| `SIM_MARGIN` | 32 | Simulator: the texture is at least n + (K-1+duty)\*max\|v\| + SIM\_MARGIN px in both directions. |
| `DEFAULT_CORR_LEN` | 1.5 | Simulator: default correlation length of the random texture in px. |
| `GREY_WEIGHTS` | 0.2989, 0.587, 0.114 | R, G, B weights of the conversion of colour images to grey (those of MATLAB rgb2gray). |

[rework] With these constants the Python and the MATLAB/Octave package give
the same results to $10^{-9}$ on shared parity cases (checked in GNU Octave
and in MATLAB), which
[CONTRIBUTING.md](../CONTRIBUTING.md#parity-cases) describes together with how
the cases are regenerated. A change to the algorithm therefore changes both
packages together.

## Reference

Jan Horn: *Zweidimensionale Geschwindigkeitsmessung texturierter Oberflächen
mit flächenhaften bildgebenden Sensoren* (two-dimensional velocity measurement
of textured surfaces with area image sensors). Doctoral thesis, Universität
Karlsruhe (TH), 2006. Schriftenreihe des Instituts für Mess- und
Regelungstechnik, Band 006, Universitätsverlag Karlsruhe, 2007. ISBN
978-3-86644-076-0.
[Full text (CC BY-NC-ND 2.0 DE)](https://publikationen.bibliothek.kit.edu/1000006276)
