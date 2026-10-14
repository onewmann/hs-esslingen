# Measurements

This folder is for measurement sessions: one subfolder per session, each
with a `meta.csv` that describes every burst of frames or image. Rows with
mode `average` average a burst of K frames into one image; rows with mode
`flash` (the default) hold one multi-flash image each. The two sessions
`synthetic_average/` and `synthetic_flash/` are simulated examples of the
format, not measurements. `experiments/real_measurements.py` estimates every
row of a session, converts the results to mm/s and compares them with the set
motion. What the camera set-up needs, how to measure the ground sample
distance and which controls to record is in
[docs/acquisition.md](../../docs/acquisition.md).

```
data/measurements/
  README.md
  synthetic_average/       simulated frame bursts for average mode, written by synthetic_average/make_example.py
  synthetic_flash/         simulated multi-flash images, written by synthetic_flash/make_example.py
  2026-10-20_average/      one folder per session, named <date>_<mode> (any name works)
    meta.csv
    burst_0001.npz
    burst_0002.npz
    ...
  2026-10-14_flash/
    meta.csv
    img_0001.png
    static_0001.png
    ...
```

```bash
python experiments/real_measurements.py data/measurements/2026-10-20_average/meta.csv --figure data/measurements/2026-10-20_average/results.png
python experiments/real_measurements.py data/measurements/synthetic_average/meta.csv --out results.csv
```

On Windows write the paths with backslashes if you like
(`python experiments\real_measurements.py data\measurements\2026-10-20_average\meta.csv`),
and put a path that contains spaces in double quotes
(`python experiments\real_measurements.py "D:\Messungen\20 Oct\meta.csv"`).
The script needs the package and its extras:
`python -m pip install -e "./python[plot,data]"` (scipy for the confidence
bounds and .mat files, Pillow for PNG and TIFF, h5py for MATLAB v7.3 files,
matplotlib for `--figure`).

It writes `results.csv` next to `meta.csv` (or to `--out FILE`) and prints a
summary: the share of moving images marked valid, the median and 95th
percentile of the speed and angle errors of the valid ones, the number of
valid images that are off by more than 5 % or 2 deg, the number of static
images marked valid (it should be 0) with a Clopper-Pearson upper bound, and
the worst images. In average mode it adds the windows with uneven steps, the
comparison of the valid windows with the mean step between their frames
(median, 95th percentile and maximum, and the number off by more than the
tolerance) and the windows that did not move at the set speed. `--figure
PNG` plots estimated against set speed, the speed errors and the angle
errors; without a set speed or angle, it plots the estimate against the mean
step. Other options: `--roi N` (default side of the analysed square, 512),
`--clip-outliers` (against glints), `--tol-speed` and `--tol-angle` (what
counts as off, default 0.05 and 2 deg), `--worst N`, `--jobs N` (worker
processes; a 512 px square takes several seconds).

Try it on `synthetic_average/` first. Without `--out` the results land in
`synthetic_average/results.csv`. For multi-flash images, see
[Multi-flash images](#multi-flash-images-mode-flash) and `synthetic_flash/`.

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

One row per burst or image, with a header row. Only `file` and `exposures`
are required (`exposures` may be left out when there is a `mode` column);
leave a cell empty when the value is not known. These columns mean the same
in both modes:

| Column | Unit | Meaning |
|---|---|---|
| `file` | | burst or image file, relative to the folder of `meta.csv` (or an absolute path): `.npz`, `.npy`, `.mat`, `.png`, `.tif` or `.bmp` (see [Files](#files)) |
| `mode` | | optional: `average` for a burst of [averaged frames](#averaged-frames-mode-average), or `flash` (empty: flash) for a [multi-flash image](#multi-flash-images-mode-flash) |
| `gsd_mm_per_px` | mm/px | ground sample distance: length on the surface covered by one pixel (how to measure it: [docs/acquisition.md](../../docs/acquisition.md#ground-sample-distance-and-physical-units)). Empty: speeds in px per step only |
| `set_speed_mm_s` | mm/s | speed set on the axis, printer or conveyor. Empty if unknown |
| `set_angle_deg` | deg | set direction of motion in the image (see [Angle convention](#angle-convention)). Empty if unknown |
| `static` | | 1 for a recording of the stopped surface (or another one that must come out not valid, see [static images](../../docs/acquisition.md#static-images)), else 0 (empty: 0) |
| `roi` | px | side of the analysed central square; empty: `--roi` (default 512); 0: the largest central square |
| `note` | | free text |

The columns that describe the timing depend on the mode: `frames`,
`first_frame`, `frame_period_s`, `exposure_s` and `duty` for
[averaged frames](#metacsv-in-average-mode), and `exposures`,
`flash_period_s`, `duty`, `exposure_s` and the reference frames for
[multi-flash images](#metacsv-in-flash-mode).

* Comma-separated. Semicolon-separated files, as German Excel saves them,
  work as well, and then decimal commas (`0,0512`) are read too. UTF-8 and
  the Windows encoding of Excel are both read. German Excel opens a
  comma-separated file in columns through Data > From Text/CSV (a double
  click puts each row into one cell); saved as CSV again, it has semicolons
  and decimal commas, and the live scripts append to it in that form.
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
  gets a warning, because it then counts twice in the summary. In average
  mode, so do two rows whose windows share frames of the same stack.

### Files

* Lossless formats only: `.npz` (the image or stack under the name
  `frames`, or its only array), `.npy`, `.mat` with one numeric array
  (MATLAB v7.3 files need h5py; a `.mat` may also hold timestamps, see
  below), PNG or TIFF with 8 or 16 bit, or BMP. Colour image files are
  converted to grey (16-bit colour only from TIFF); convert a colour array in
  a `.npy`, `.npz` or `.mat` file to grey first, otherwise its channels count
  as frames (the script warns). JPEG compression adds an 8 x 8 px block grid
  that can give wrong results marked valid; the script warns about JPEG
  files.
* Rows with mode `average` average consecutive frames of a stack into one
  pulsed image ([averaged frames](#averaged-frames-mode-average)): a MATLAB
  `getdata` recording (rows x columns x 1 x frames, saved with
  `save('rec1.mat', 'recording')`), a 3-D `.npy` or `.npz` or a multi-page
  TIFF.
* In flash mode, a file that holds several frames gives one result per
  frame, each estimated on its own and named `file#1`, `file#2`, ...; every
  frame must then be a multi-flash image. `rec1.mat#3` in the `file` column
  selects frame 3 (counted from 1, as in MATLAB).
* The analysed square is the central `roi` x `roi` px of the image. A
  camera ROI that is already square and centred on the surface does the
  same and saves transfer time.

### Angle convention

Directions are in degrees, counter-clockwise from the +x axis (the column
index, to the right) with y pointing up on the screen, the image shown with
its first row at the top (as `imshow` in MATLAB shows it). One averaged or
multi-flash image cannot tell a direction from its opposite, so all angles
are taken modulo 180 deg:

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
images of the surface between which the axis moved (two frames of a burst
of the moving surface, or two images of the stopped surface before and
after a move):

```python
import numpy as np
from fftvel import phase_correlation
dx, dy, peak = phase_correlation(a, b)          # a, b: 2-D arrays of the same size
angle = np.degrees(np.arctan2(-dy, dx)) % 180   # set_angle_deg
```

In MATLAB: `[dx, dy] = fftvel.phase_correlation(a, b); angle = mod(atan2d(-dy, dx), 180)`.

## Averaged frames (mode average)

The camera takes K consecutive short exposures at a fixed frame rate under
continuous light, and the mean of these K frames is the pulsed image. If the
exposures are short and the steps between the frames are equal, it is the
image that K flashes, one per frame, would give: the frame period takes the
place of the flash period, and the exposure divided by the frame period the
place of the duty. Rows of `meta.csv` with mode `average` average the frames;
both live scripts record them. The camera settings and how to choose the
frame rate from the speed are in
[docs/acquisition.md](../../docs/acquisition.md#frame-bursts).

Results of this mode test the estimator on a real camera, lens, texture and
noise, with the real GSD and set motion. They do not test what only a
multi-flash image involves: the flash timing, the number of flashes in one
exposure, or other light during one long exposure (see the
[checklist](../../docs/acquisition.md#checklist)). Report them separately from
multi-flash measurements.

### Recording

```bash
python python/examples/live_camera.py --mode average --source pylon --exposures 9 --frame-rate 80 --exposure-us 2000 --save data/measurements/2026-10-20_average --note "axis 20 mm/s"
```

```matlab
addpath('matlab', 'matlab/examples')       % from the repository root
live_basler('Mode', 'average', 'Exposures', 9, 'FrameRate', 80, 'ExposureUs', 2000, ...
    'SaveDir', 'data/measurements/2026-10-20_average', 'Note', 'axis 20 mm/s')
```

`--exposure-us` and `'ExposureUs'` are the exposure time in µs; the folder is
one session. Python needs pypylon for Basler cameras
(`python -m pip install pypylon`). In MATLAB the figure has a Start/Stop
button: recording starts with Start, Stop pauses it, and closing the figure
releases the camera.

Every burst of K frames that the script estimates is saved, as
`burst_0001.npz` (Python) or `burst_0001.mat` (MATLAB), numbered on from the
files already in the folder. A `.npz` holds `frames` (rows x columns x K, in
the data type of the camera), `t` (the timestamps in s), `clock`
(`camera` or `host`), and `exposure_s` and `frame_rate` when known; a `.mat`
holds `frames` as `getdata` returns them (rows x columns x 1 x K), `t`,
`exposure_s`, `frame_rate` and `clock` (`host`). A line for each file is
appended to `meta.csv` in the folder, with the mode, `frames`, `exposure_s`
and the note filled in. A `meta.csv` saved from Excel keeps its header and
separator; with semicolons the numbers are written with decimal commas, as
German Excel expects them. Fill in `gsd_mm_per_px`, `set_speed_mm_s`,
`set_angle_deg` and `static` before running the script. Record windows of
the stopped surface with the same settings as well and give them `static`
1 (a separate run with `--note "axis stopped"` makes them easy to find);
they must come out not valid. In Python, `--static` writes the 1 for you, and
`--count N` stops after N bursts: `--count 1 --static --note "patch 7"`
records one static burst per patch
([first recording](../../docs/acquisition.md#first-recording-step-by-step)).
In MATLAB, press Stop after the burst.

A burst in which a frame was lost is neither estimated nor saved; Python
prints `burst skipped` for it and `grab failed` for each failed grab. Count
these lines and give the number in the report, since every recorded burst
counts ([reporting](#reporting-the-result)).

In Python, `--roi` only sets the analysed square of the live estimate; the
saved files hold the frames as the camera delivers them, and
`real_measurements.py` uses its own `--roi`. A smaller live `--roi` makes the
estimate faster, so bursts are saved more often. To make the files smaller,
set a camera ROI (Width, Height and the offsets in pylon Viewer); `live_basler`
crops to its `'ROI'` itself. A full-sensor burst of 9 frames of 1936 x 1216 px
in Mono8 is 21 MB before compression (Python saves compressed `.npz` files),
a 512 x 512 px burst 2.4 MB.

Decide before committing a session what goes into the repository and gets
published: all bursts, or only `meta.csv`, `results.csv`, the figure and a
few bursts, with the rest stored elsewhere. Committed data fall under the
data licence ([data/README.md](../README.md#licence-of-the-data)). Folder names, file names and
notes are published with them.

Python takes the timestamps from the camera clock when the camera sends
them. MATLAB `getdata` gives times taken on the host when a frame arrives,
which jitter by the transfer; for these the line in `meta.csv` gets the
camera's frame period in `frame_period_s`. One estimate on a 512 px square
takes a few seconds, so the scripts save one burst every few seconds.

### meta.csv in average mode

| Column | Unit | Meaning |
|---|---|---|
| `mode` | | `average` |
| `file` | | a stack: `.npz`, `.mat` (from `getdata`, rows x columns x 1 x frames, or rows x columns x frames), `.npy` or a multi-page TIFF. The frames lie on the last axis (or on the first, if that is clearly the short one). `file#k` is not allowed here |
| `frames` | | number of consecutive frames averaged, at least 2 (also even). Empty: `exposures`, else all frames of the stack, if their number is odd. The averaged image is estimated with this number as K |
| `first_frame` | | first frame averaged, counted from 1 (empty: 1). A longer stack can be listed in several rows, for example with `first_frame` 1, 10, 19 and `frames` 9 |
| `exposures` | | optional; if given, it must equal `frames` |
| `frame_period_s` | s | time from one frame to the next. Empty: `flash_period_s` if given, else for host times the period of the `frame_rate` stored in the file, else the period measured from the timestamps in the file (see below). Empty without timestamps: speeds in px per frame only |
| `exposure_s` | s | exposure of one frame; gives duty = `exposure_s` / frame period unless `duty` is given. Empty: the `exposure_s` stored in the file, if any |
| `duty` | | exposure / frame period, in [0, 1); empty: from `exposure_s`, else 0. In `results.csv`: the duty used, also when it comes from `exposure_s` (and `exposure_s` the exposure used, also when it comes from the file) |

The other columns (`gsd_mm_per_px`, `set_speed_mm_s`, `set_angle_deg`,
`static`, `roi`, `note`) mean the same in both modes ([meta.csv](#metacsv)); `ref_frame_a` and
`ref_frame_b` are not used, because the frames themselves give the
reference (step check below). `flash_period_s` is accepted in place of
`frame_period_s`; if both are given they must be equal.

```
file,mode,frames,first_frame,frame_period_s,exposure_s,gsd_mm_per_px,set_speed_mm_s,set_angle_deg,static,note
burst_0001.npz,average,9,,,0.002,0.05,20,0,0,axis at 20 mm/s
long.npy,average,9,10,0.0125,0.002,0.05,20,0,0,frames 10 to 18 of a longer stack
burst_0031.npz,average,9,,,0.002,0.05,,,1,axis stopped
```

Timestamps: a `.npz` entry or `.mat` variable named `t`, `timestamps` or
`time` (the name in `[data, time] = getdata(vid)`), in s, one per frame of
the stack (`.npy` and TIFF files hold none). The intervals of the averaged
frames are checked as in the live scripts: their spread, max - min, must be
at most 10 % of their median. The measured period (`measured_period_s`) is
the median interval, which a dropped frame does not change.

Times taken on the host jitter by the transfer: those with `clock` `host` in
the file, and the times in a `.mat` file without `clock`, since `getdata`
gives host times. For them the measured period is the mean interval of the
window, and they are accepted when that mean agrees with the nominal period
within 10 % and no interval is off by half a period, as in `live_basler`.
The nominal period is `frame_period_s` (or `flash_period_s`), else the
period of the `frame_rate` stored in the file, which then also gives the
speeds in mm/s, as in the live scripts. Without a nominal period the spread
rule above applies to host times too, and their jitter alone can then flag
a window; give the camera's frame period in `frame_period_s`.

The script warns if the measured period differs from the period used by
more than 10 %. `period_source` in `results.csv` says which period the
speeds in mm/s use: `frame_period_s`, `flash_period_s`, `frame_rate` (host
times) or `timestamps`.

### Step check

Phase correlation (`fftvel.phase_correlation`) of every pair of consecutive
averaged frames, on the analysed square, gives the step (dx, dy) in px from
one frame to the next. Their mean is a second measurement of the motion per
frame from the same frames, independent of the estimator:
`step_dx_px`, `step_dy_px` and `step_speed_px` (its length) in px per frame.
It fills the reference columns as the reference frames do in flash mode:
`ref_speed_mm_s` is the mean step in mm/s, `ref_speed_diff` and
`ref_angle_diff_deg` compare the estimate with it, and the summary prints
their median, 95th percentile and maximum over the valid windows with even
steps, and how many of these are off by more than `--tol-speed` or
`--tol-angle`. This comparison depends neither on the axis nor on the GSD.

With a set speed, `step_vs_set` is the mean step divided by the set speed in
px per frame (`expected_px`), minus 1. The summary counts the windows with
even steps whose `step_vs_set` lies beyond `--tol-speed`, or whose frames do
not move although the set speed says they should (such a window also gets a
warning): the axis did not run at the set speed during these bursts, for
example while it paused or reversed. A valid window that is off against the
set motion but agrees with the mean step points to the reference (set speed
or angle, GSD or frame period), not to the estimate; the summary names their
number. Phase correlation
gives a signed direction (dy points down, so the direction is atan2(-dy,
dx)), while the estimate gives one modulo 180 deg; the difference is taken
with `fftvel.angle_difference`.

The steps are uneven (`uneven_steps` 1) when their lengths spread (max -
min) by more than 10 % of their median (`step_spread`) or their directions
by more than 5 deg (`step_angle_spread_deg`). Causes are a dropped frame
(one step twice as long), a speed that changes within the window, or
vibration. `uneven_steps` is also 1 when the timestamps show uneven frame
intervals; the `warnings` column says which check found it, and the summary
counts the windows with uneven steps. Steps shorter than 0.25 px are too
short to check; for them the steps are reported but not compared with the
estimate, and a static window whose frames move by more than that gets a
warning.

The estimate assumes equal steps. A window with uneven steps can still come
out valid: in `synthetic_average/`, the stack with one missing frame gives the
set speed (7 of its 8 steps are equal), while its mean step is 12.5 % too
long. Decide before looking at the estimates whether windows with uneven
steps are excluded, and give the reason (for example, a dropped frame) in
`note`.

## Multi-flash images (mode flash)

Rows without a mode, or with mode `flash`, hold one multi-flash image each:
a strobe fired K times during one camera exposure. The set-up is in
[docs/acquisition.md](../../docs/acquisition.md#multi-flash-images).

### meta.csv in flash mode

| Column | Unit | Meaning |
|---|---|---|
| `exposures` | | number of flashes K in the image, at least 2 (also even) |
| `flash_period_s` | s | time from the start of one flash to the start of the next. Empty: speeds in px per flash only |
| `duty` | | flash pulse length divided by the flash period, in [0, 1). 0 or empty if unknown |
| `exposure_s` | s | optional: camera exposure time; the script warns if it is shorter than (K - 1 + duty) flash periods |
| `ref_frame_a`, `ref_frame_b` | | optional: two single-flash frames for the [phase-correlation cross-check](#phase-correlation-cross-check) |
| `ref_dt_s` | s | optional: time from `ref_frame_a` to `ref_frame_b`; empty: only the directions are compared |

```
file,exposures,flash_period_s,duty,gsd_mm_per_px,set_speed_mm_s,set_angle_deg,static,roi,exposure_s,note
img_0001.png,9,0.0005,0.04,0.0512,500,0,0,,0.0045,axis at 500 mm/s
img_0002.png,9,0.0005,0.04,0.0512,250,0,0,,0.0045,axis at 250 mm/s
static_0001.png,9,0.0005,0.04,0.0512,,,1,,0.0045,axis stopped
```

### Phase-correlation cross-check

`ref_frame_a` and `ref_frame_b` name two single-flash frames (one flash per
frame) of the moving surface at the same speed, `ref_dt_s` apart. Phase
correlation of their central squares (the same size as the analysed square)
gives the displacement (dx, dy) in px. The script reports it, its direction,
the reference speed in px per flash, |d| * T / `ref_dt_s`, and in mm/s,
|d| * GSD / `ref_dt_s`, and the differences to the estimate. The comparison
in px per flash does not depend on the GSD or on the accuracy of the axis,
only on the timing. Keep |d| well below half of the analysed square; the
script warns above a quarter. The frames can be in a stack: `ref.mat#1` and
`ref.mat#2`.

## results.csv

One row per burst or image (per frame for stacks in flash mode), in the
order of `meta.csv`, with the columns of `meta.csv` and the ones below. The
first columns are `file`, `line`, `valid`, `reason`, `speed_mm_s`,
`set_speed_mm_s`, `speed_err`, `angle_deg`, `set_angle_deg` and
`angle_err_deg`; the settings, the cross-check, `warnings`, `error` and
`note` follow. The file is comma-separated with a decimal point, in UTF-8.
German Excel splits it into columns when it is opened through Data > From
Text/CSV (delimiter comma); a double click puts each row into one cell.

| Column | Meaning |
|---|---|
| `line` | line of the image in `meta.csv` |
| `roi` | side of the analysed square in px |
| `angle_deg`, `speed_px`, `vx_px`, `vy_px` | estimate in px per step: per frame in average mode, per flash in flash mode (vx along the column index, vy along the row index, pointing down). One image cannot tell a direction from its opposite, so the motion may be (vx, vy) or (-vx, -vy); `vx_px` can be negative for motion to the right |
| `speed_mm_s`, `vx_mm_s`, `vy_mm_s` | the same in mm/s, from `fftvel.to_physical` |
| `quality`, `quality_perp`, `threshold`, `geba` | fit quality, quality of the perpendicular fit, quality needed, ripple ratio |
| `in_range`, `isotropic`, `valid` | 1 or 0, as in the estimate |
| `quality_blur`, `blur_length`, `blur_like` | fit quality and length of a continuous blur along the found direction; `blur_like` is 1 when that fits at least as well and the speed is below about 2.2 px per step, where one blurred exposure and K exposures cannot be told apart |
| `saturated` | share of the analysed pixels that are saturated |
| `reason` | why the estimate is not valid |
| `expected_px` | set speed in px per step, `set_speed_mm_s` times the frame or flash period divided by `gsd_mm_per_px` |
| `speed_err` | `speed_mm_s / set_speed_mm_s - 1` (signed) |
| `angle_err_deg` | estimated minus set direction, modulo 180, in [-90, 90) |
| `mode` | `flash` or `average` |
| `frames`, `first_frame` | average mode: frames averaged (also in `exposures`) and the first of them |
| `frame_period_s` | average mode: as in `meta.csv` |
| `measured_period_s`, `period_source` | average mode: frame period measured from the timestamps in the file, and which period the speeds in mm/s use (`frame_period_s`, `flash_period_s`, `frame_rate` or `timestamps`) |
| `step_dx_px`, `step_dy_px`, `step_speed_px` | average mode: mean step between consecutive frames by phase correlation, in px per frame, and its length; the `ref_` columns compare the estimate with it |
| `step_vs_set` | average mode: `step_speed_px / expected_px - 1`, the mean step against the set speed (empty without a set speed, or with steps below 0.25 px) |
| `step_spread`, `step_angle_spread_deg`, `uneven_steps` | average mode: spread of the step lengths (relative to their median) and of their directions, and 1 if the steps or the frame intervals are uneven (see [Step check](#step-check)) |
| `ref_dx_px`, `ref_dy_px`, `ref_peak`, `ref_angle_deg` | flash mode: phase correlation of the reference frames (average mode: the mean step) |
| `ref_speed_px`, `ref_speed_mm_s` | reference speed in px per step and mm/s |
| `ref_speed_diff`, `ref_angle_diff_deg` | `speed_px / ref_speed_px - 1` and estimated minus reference direction |
| `warnings` | exposure too short, set speed outside the measurable range, saturation, JPEG, large reference displacement, a stack of frames read in flash mode; in average mode also uneven steps or frame intervals, frames that do not move although the set speed says they should, a duty above 0.25, a measured period that differs from the one used, windows that share frames |
| `error` | why the row was skipped (`row skipped: ...`, a problem in `meta.csv`) or the file could not be read or estimated |

## Reporting the result

* Report every burst or image that was recorded. Decide about exclusions
  before looking at the estimates and only for reasons that do not depend on
  them (for example, the axis was still accelerating, or a frame was
  dropped); comment those rows out with `#` and give the reason in `note`.
  The summary counts them.
* Give the share of moving images marked valid and the reasons for the
  others (column `reason`), the number of images behind each statistic, the
  median, 95th percentile and maximum of the errors of the valid images
  (with fewer than 20 images the 95th percentile lies between the largest
  values), the valid images that are off by more than 5 % or 2 deg, and the
  static images marked valid, each with its upper bound. The bound assumes
  independent images, so count static images of one surface patch as one
  ([static images](../../docs/acquisition.md#static-images)).
* In average mode, also give the comparison with the mean step between the
  frames. It tests the estimator alone, while the comparison with the set
  speed also contains the errors of the axis, the GSD and the timing.
* State the uncertainty of the reference: the set speed against the actual
  speed of the axis (encoder), the GSD and the frame interval or flash
  period. Their relative errors add to the speed error in mm/s; a smaller
  error than theirs cannot be claimed.
* Give the settings: K, frame rate or flash period, exposure or pulse
  length, ROI, camera, lens, working distance, illumination and surface.

## Example

`synthetic_average/` holds a synthetic example of average mode: three stacks
of 9 frames of 192 x 192 px, 8-bit, as `.npz` files in the format that
`live_camera.py --save` writes. The frames are the same texture shifted by
6 px per frame at 30 deg (5 ms frame period, 0.5 ms exposure, noise of 2
counts per frame): `moving_regular.npz` with regular timestamps and steps,
`moving_dropped.npz` with one frame missing (one step and one frame interval
twice as long, flagged by the step check), and `static.npz` of the stopped
surface. `synthetic_average/make_example.py` writes them. The script marks
both moving windows valid with speed errors below 0.1 %, flags the one with
the missing frame, and marks the static one not valid:

```bash
python experiments/real_measurements.py data/measurements/synthetic_average/meta.csv --out results.csv
```

`synthetic_flash/` holds a synthetic example of flash mode: three
multi-flash images (9 flashes at 4 px per flash, duty 0.02, noise 10 %; 8
flashes at 5.6 px per flash, duty 0.2, noise 30 %; a stopped surface) and two
single-flash reference frames, 256 x 256 px, 16-bit PNG.
`synthetic_flash/make_example.py` writes them with
`fftvel.simulate_pulsed_image` from `data/surface_texture.png` (a real Basler
texture) with exactly known motion, so the errors it shows (below 0.1 % and
0.2 deg) come from the estimator alone, without the errors of a real axis,
GSD or flash timing.

```bash
python experiments/real_measurements.py data/measurements/synthetic_flash/meta.csv --out results.csv
```
