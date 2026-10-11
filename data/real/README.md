# Real multi-flash measurements

This folder holds real multi-flash images, one subfolder per measurement
session, each with a `meta.csv` that describes every image.
`experiments/real_measurements.py` estimates all of them, converts the
results to mm/s and compares them with the set motion:

```
data/real/
  README.md
  example/                 synthetic example, written by example/make_example.py
  2026-10-14_axis/         one folder per session (any name)
    meta.csv
    img_0001.png
    img_0002.png
    static_0001.png
    ...
```

```bash
python experiments/real_measurements.py data/real/2026-10-14_axis/meta.csv --figure data/real/2026-10-14_axis/results.png
python experiments/real_measurements.py data/real/example/meta.csv --out results.csv
```

On Windows write the paths with backslashes if you like
(`python experiments\real_measurements.py data\real\2026-10-14_axis\meta.csv`),
and put a path that contains spaces in double quotes
(`python experiments\real_measurements.py "D:\Messungen\14 Oct\meta.csv"`).
The script needs the package and its extras:
`python -m pip install -e "./python[plot,data]"` (scipy for the confidence
bounds and .mat files, Pillow for PNG and TIFF, h5py for MATLAB v7.3 files,
matplotlib for `--figure`).

It writes `results.csv` next to `meta.csv` (or to `--out FILE`) and prints a
summary: the share of moving images marked valid, the median and 95th
percentile of the speed and angle errors of the valid ones, the number of
valid images that are off by more than 5 % or 2 deg, the number of static
images marked valid (it should be 0) with a Clopper-Pearson upper bound, and
the worst images. `--figure PNG` plots estimated against set speed, the speed
errors and the angle errors. Other options: `--roi N` (default side of the
analysed square, 512), `--clip-outliers` (against glints), `--tol-speed` and
`--tol-angle` (what counts as off, default 0.05 and 2 deg), `--worst N`,
`--jobs N` (worker processes; a 512 px square takes several seconds).

Try it on `example/` first. Without `--out` the results land in
`example/results.csv`.

If a row of `meta.csv` has a problem (see below), the script prints it with
its line, skips that row and analyses the others; the skipped rows appear in
`results.csv` with the reason in `error` and are listed again at the end.
The exit status is 0 when every row was analysed, 1 when some were not
(also when none were; `results.csv` is then not written), and 2 when
`meta.csv` or the output cannot be used at all (no `file` or `exposures`
column, or `results.csv` cannot be written). On Windows close
`results.csv` in Excel before running the script again: Excel locks the
file, and the script then stops before computing anything.

## meta.csv

One row per image, with a header row. Only `file` and `exposures` are
required; leave a cell empty when the value is not known.

| Column | Unit | Meaning |
|---|---|---|
| `file` | | image file, relative to the folder of `meta.csv` (or an absolute path): `.png`, `.tif`, `.bmp`, `.npy` or `.mat` (see [Files](#files)) |
| `exposures` | | number of flashes K in the image, at least 2 (also even) |
| `flash_period_s` | s | time from the start of one flash to the start of the next. Empty: speeds in px per flash only |
| `duty` | | flash pulse length divided by the flash period, in [0, 1). 0 or empty if unknown |
| `gsd_mm_per_px` | mm/px | ground sample distance: length on the surface covered by one pixel (see [gsd](#ground-sample-distance)). Empty: speeds in px per flash only |
| `set_speed_mm_s` | mm/s | speed set on the axis or turntable. Empty if unknown |
| `set_angle_deg` | deg | set direction of motion in the image (see [Angle convention](#angle-convention)). Empty if unknown |
| `static` | | 1 for an image of the stopped surface (or another image that must come out not valid, see [Static images](#static-images)), else 0 (empty: 0) |
| `roi` | px | side of the analysed central square; empty: `--roi` (default 512); 0: the largest central square |
| `note` | | free text |
| `exposure_s` | s | optional: camera exposure time; the script warns if it is shorter than (K - 1 + duty) flash periods |
| `ref_frame_a`, `ref_frame_b` | | optional: two single-flash frames for the [phase-correlation cross-check](#phase-correlation-cross-check) |
| `ref_dt_s` | s | optional: time from `ref_frame_a` to `ref_frame_b`; empty: only the directions are compared |

```
file,exposures,flash_period_s,duty,gsd_mm_per_px,set_speed_mm_s,set_angle_deg,static,roi,exposure_s,note
img_0001.png,9,0.0005,0.04,0.0512,500,0,0,,0.0045,axis at 500 mm/s
img_0002.png,9,0.0005,0.04,0.0512,250,0,0,,0.0045,axis at 250 mm/s
static_0001.png,9,0.0005,0.04,0.0512,,,1,,0.0045,axis stopped
```

* Comma-separated, as above. Semicolon-separated files, as German Excel
  saves them, work as well, and then decimal commas (`0,0512`) are read too.
  UTF-8 and the Windows encoding of Excel are both read.
* Empty rows are skipped. A row whose `file` starts with `#` is skipped and
  counted as commented out in the summary; give the reason in `note`.
* Columns with other names are allowed; the script lists them as not used,
  which also catches typos in the names above.
* `set_speed_mm_s` 0 is only accepted for `static` 1.
* Before anything is computed, all rows are checked (numbers, ranges, files
  that do not exist), and every problem is reported with its line number.
  A row with more cells than the header is a problem too: in a
  comma-separated file a decimal comma (`0,2`) or a comma in an unquoted
  note shifts the following values into the wrong columns. Rows with
  problems are skipped; the others are analysed.
* The same image listed twice (for example `rec1.mat` and `rec1.mat#3`)
  gets a warning, because it then counts twice in the summary.

### Files

* Lossless formats only: PNG or TIFF with 8 or 16 bit, BMP, `.npy`, or `.mat`
  with one numeric array (MATLAB v7.3 files need h5py). Colour images are
  converted to grey. JPEG compression adds an 8 x 8 px block grid that can
  give wrong results marked valid; the script warns about JPEG files.
* A file that holds several frames gives one result per frame, each
  estimated on its own and named `file#1`, `file#2`, ...: a MATLAB `getdata`
  recording (rows x columns x 1 x frames, saved with
  `save('rec1.mat', 'recording')`), a 3-D `.npy` or a multi-page TIFF.
  Every frame must then be a multi-flash image. `rec1.mat#3` in the `file`
  column selects frame 3 (counted from 1, as in MATLAB). Frames are never
  averaged here; to average single exposures into a pulsed image use
  `python -m fftvel estimate FILE --stack` instead.
* The analysed square is the central `roi` x `roi` px of the image. A
  camera ROI that is already square and centred on the surface does the
  same and saves transfer time.

### Angle convention

Directions are in degrees, counter-clockwise from the +x axis (the column
index, to the right) with y pointing up on the screen, the image shown with
its first row at the top (as `imshow` in MATLAB shows it). One
image cannot tell a direction from its opposite, so all angles are taken
modulo 180 deg:

```
                    90
                    |
          135  \    |    /  45
                \   |   /
                 \  |  /
     0 = 180 ------ + ------> x  0       first image row at the top
                 /  |  \
                /   |   \
           45  /    |    \  135
                    |
                    90
```

Motion to the right or to the left is 0, up or down on the screen 90, up
and to the right 45, down and to the right 135. The set direction depends
on how the camera is mounted relative to the axis. Measure it once with two
images of the surface between which the axis moved (two single-flash
frames of the moving surface, or two images of the stopped surface before
and after a move):

```python
import numpy as np
from fftvel import phase_correlation
dx, dy, peak = phase_correlation(a, b)          # a, b: 2-D arrays of the same size
angle = np.degrees(np.arctan2(-dy, dx)) % 180   # set_angle_deg
```

In MATLAB: `[dx, dy] = fftvel.phase_correlation(a, b); angle = mod(atan2d(-dy, dx), 180)`.

### Ground sample distance

`gsd_mm_per_px` is the length on the moving surface that one pixel covers.
Its relative error passes one to one into every speed in mm/s, so measure it
carefully:

1. Put a checkerboard or a ruler in the plane of the moving surface, at the
   same height, with the camera, lens, focus and working distance of the
   measurements. Refocusing or moving the camera changes the gsd.
2. Take an image and find the pixel coordinates of two marks far apart, for
   example 80 mm on the ruler or 10 squares of the checkerboard, near the
   image centre where the analysed square lies (lens distortion changes the
   gsd towards the corners). In MATLAB: `imshow(I); [x, y] = ginput(2);
   d = hypot(diff(x), diff(y))`.
3. gsd = known length / pixel distance. Measure along x and along y; if the
   two differ by more than about 0.5 %, the camera looks at the surface at
   an angle, and one gsd does not fit both directions.
4. Reading each mark to within 1 px over 1600 px gives about 0.1 % error;
   over 160 px it is about 1 %.

The axis can also measure it: take an image of the stopped surface, move the
axis by a known distance D (say 10 mm), take another, and get the
displacement d in px by phase correlation as above. Then gsd = D / |d|. Keep
|d| below a quarter of the analysed square (with a 1216 px square, below
about 300 px).

### Flash period

The motion per flash in px is

    v_px = v * T / gsd

with v the speed in mm/s, T the flash period in s and gsd in mm/px. The
estimator measures from about 3.3/K px per flash (0.37 px for K = 9) up to
about n/4 px for an n px square (128 px for the default 512 px); the
simulations test up to n/4. Aim for a few px per flash, about 2 to 10, so
that T = v_px * gsd / v.

Worked example: v = 500 mm/s, gsd = 0.05 mm/px (the 1936 px wide sensor
covers 97 mm), 5 px per flash: T = 5 * 0.05 / 500 s = 0.5 ms, a flash rate
of 2 kHz. With K = 9 the exposure has to cover (K - 1 + duty) * T, so at
least 8.04 * 0.5 ms = 4.0 ms; set 4.5 ms and make sure that exactly 9
flashes fall into it. A 20 us pulse gives duty = 20 us / 0.5 ms = 0.04 and a
blur of duty * v_px = 0.2 px per flash, which is negligible. With the same
T an axis at 100 to 1000 mm/s gives 1 to 10 px per flash. During the 9
flashes the surface moves (K - 1) * v_px = 40 px, a small part of the 512 px
square.

### Acquisition

* Exactly K flashes per image. The safest set-up lets the camera start a
  burst of exactly K pulses (exposure active or frame trigger output to the
  flash controller), or lets the controller trigger the camera. With a
  free-running strobe the number of flashes in an exposure changes from
  frame to frame. A wrong K scales the speed by K_true/K for slow motion,
  and the result can still be marked valid.
* Equal flash intervals: uneven ones can give twice the speed, marked valid.
  Check the number of pulses and the flash period once with a photodiode
  and an oscilloscope; an error in the period passes one to one into the
  speed in mm/s.
* No other light during the exposure. The exposure lasts about (K - 1) * T,
  and continuous light during it adds a motion-blurred image of the
  surface. Any blur of b px along one direction, from such light or from
  vibration, looks like K flashes b/K px apart and can be marked valid.
  The smoke test on the 2025 recordings ran into this: single exposures of
  `vid4.mat`, read as if they held 9 flashes, gave 3 of 10 frames valid with
  0.48 to 0.54 px per flash, all at 90 deg, that is a blur of about 4.3 to
  4.8 px, vertical on the screen. The cause of that blur is not known; from
  frame to frame the surface moved 1.6 to 6.2 px in directions of 89 to
  128 deg (phase correlation of the central 512 px square). Darken
  the room or use an enclosure or a band-pass filter for the flash colour,
  check that an image with the flash switched off is black, and record the
  single-flash controls described under [Static images](#static-images).
* No saturation: the script warns when more than 0.5 % of the analysed
  pixels sit at the image maximum. Glints of the flash on shiny surfaces
  can be reported as motion; `--clip-outliers` helps against small ones.
* Focus: a defocused image can look like slow motion. The static images
  below check this.

### Static images

Record a number of images of the stopped surface with the same settings (K,
flash period, exposure, ROI) and give them `static` 1. All of them must come
out not valid. `static` 1 means "must come out not valid", so it also fits
a second control: images of the moving surface with a single flash in the
same exposure time (listed with the `exposures` of the real measurements).
They hold no flash pattern, and a valid one shows blur from other light or
vibration. With none of N static images valid, the one-sided 95 %
Clopper-Pearson upper bound of the rate of static images marked valid is
1 - 0.05^(1/N): 26 % for N = 10, 9.5 % for N = 30, 4.9 % for N = 60, 3.0 %
for N = 100. The summary prints this bound for the images given.

### Phase-correlation cross-check

`ref_frame_a` and `ref_frame_b` name two single-flash frames (one flash per
frame) of the moving surface at the same speed, `ref_dt_s` apart. Phase
correlation of their central squares (the same size as the analysed square)
gives the displacement (dx, dy) in px. The script reports it, its direction,
the reference speed in px per flash, |d| * T / `ref_dt_s`, and in mm/s,
|d| * gsd / `ref_dt_s`, and the differences to the estimate. The comparison
in px per flash does not depend on the gsd or on the accuracy of the axis,
only on the timing. Keep |d| well below half of the analysed square; the
script warns above a quarter. The frames can be in a stack: `ref.mat#1` and
`ref.mat#2`.

## results.csv

One row per image (per frame for stacks), in the order of `meta.csv`, with
the columns of `meta.csv` and the ones below. The first columns are `file`,
`line`, `valid`, `reason`, `speed_mm_s`, `set_speed_mm_s`, `speed_err`,
`angle_deg`, `set_angle_deg` and `angle_err_deg`; the settings, the
cross-check, `warnings`, `error` and `note` follow. The file is
comma-separated with a decimal point, in UTF-8. German Excel splits it
into columns when it is opened through Data > From Text/CSV (delimiter
comma); a double click puts each row into one cell.

| Column | Meaning |
|---|---|
| `line` | line of the image in `meta.csv` |
| `roi` | side of the analysed square in px |
| `angle_deg`, `speed_px`, `vx_px`, `vy_px` | estimate in px per flash (vx along the column index, vy along the row index, pointing down). One image cannot tell a direction from its opposite, so the motion may be (vx, vy) or (-vx, -vy); `vx_px` can be negative for motion to the right |
| `speed_mm_s`, `vx_mm_s`, `vy_mm_s` | the same in mm/s, from `fftvel.to_physical` |
| `quality`, `quality_perp`, `threshold`, `geba` | fit quality, quality of the perpendicular fit, quality needed, ripple ratio |
| `in_range`, `isotropic`, `valid` | 1 or 0, as in the estimate |
| `reason` | why the estimate is not valid |
| `expected_px` | set speed in px per flash, `set_speed_mm_s * flash_period_s / gsd_mm_per_px` |
| `speed_err` | `speed_mm_s / set_speed_mm_s - 1` (signed) |
| `angle_err_deg` | estimated minus set direction, modulo 180, in [-90, 90) |
| `ref_dx_px`, `ref_dy_px`, `ref_peak`, `ref_angle_deg` | phase correlation of the reference frames |
| `ref_speed_px`, `ref_speed_mm_s` | reference speed in px per flash and mm/s |
| `ref_speed_diff`, `ref_angle_diff_deg` | `speed_px / ref_speed_px - 1` and estimated minus reference direction |
| `warnings` | exposure too short, set speed outside the measurable range, saturation, JPEG, large reference displacement |
| `error` | why the row was skipped (`row skipped: ...`, a problem in `meta.csv`) or the file could not be read or estimated |

## Reporting the result

* Report every image that was recorded. Decide about exclusions before
  looking at the estimates and only for reasons that do not depend on them
  (for example, the axis was still accelerating); comment those rows out
  with `#` and give the reason in `note`. The summary counts them.
* Give the share of moving images marked valid and the reasons for the
  others (column `reason`), the number of images behind each statistic, the
  median, 95th percentile and maximum of the errors of the valid images
  (with fewer than 20 images the 95th percentile lies between the largest
  values), the valid images that are off by more than 5 % or 2 deg, and the
  static images marked valid, each with its upper bound.
* State the uncertainty of the reference: the set speed against the actual
  speed of the axis (encoder), the gsd and the flash period. Their relative
  errors add to the speed error in mm/s; a smaller error than theirs cannot
  be claimed.
* Give the settings: K, flash period, pulse length, exposure, ROI, camera,
  lens, working distance, illumination and surface.

## Example

`example/` holds a synthetic example to try the script on before real images
exist: three multi-flash images (9 flashes at 4 px per flash; 8 flashes at
5.6 px per flash with duty 0.2 and noise; a stopped surface) and two
single-flash reference frames, 256 x 256 px, 16-bit PNG.
`example/make_example.py` writes them with `fftvel.simulate_pulsed_image`
from `data/surface_texture.png` (a real Basler texture) with exactly known
motion, so the errors it shows (below 0.1 % and 0.2 deg) come from the
estimator alone, without the errors of a real axis, gsd or flash timing.
