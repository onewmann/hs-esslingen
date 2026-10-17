# Graph paper: a periodic texture

Two phone videos of graph paper with 5 mm squares, recorded in slow motion at
120 frames per second: one with the phone moved by hand, one with the phone
held still over the paper. A regular grid breaks assumption 1 of
[method.md](../../docs/method.md#assumptions) (a texture that does not repeat),
so this session shows how the estimator behaves on a real periodic texture.
The grid also gives the scale: the squares turn pixels into millimetres. Data
and figure are under CC BY 4.0
([licence of the data](../../README.md#licence-of-the-data)).

| File | What |
|---|---|
| `grid_moving_120fps.mp4` | phone moved by hand, 546 frames, 4.7 s; the video stream as recorded, without sound track and file metadata |
| `grid_still_120fps.mp4` | phone held still, 447 frames, 3.9 s; the same |
| `extract_frames.py` | writes `frames_moving.npz` and `frames_still.npz` (the central 512 x 512 px square of every frame, with the timestamps) and `meta.csv`, with the ground sample distance of each window of the moving video measured on the grid |
| `meta.csv` | 109 windows of 5 frames of the moving video, 89 of the still one (`static` 1) |
| `results.csv`, `results.png` | the output of `experiments/real_measurements.py` |

```bash
python -m pip install -e "./python[plot,data,video]"
python data/measurements/2026-10-17_phone_grid/extract_frames.py
python experiments/real_measurements.py data/measurements/2026-10-17_phone_grid/meta.csv --figure data/measurements/2026-10-17_phone_grid/results.png
```

## Scale

A square is 122 to 138 px wide in the moving video, measured on the 1024 px
field around the analysed square; the phone's height changed as the hand
moved. Each window of the moving video therefore gets its own ground sample
distance, 5 mm divided by the square's width, between 0.036 and 0.041 mm/px.
The last two windows, where the paper's edge enters the field, get none. In
the still video the squares grow from about 168 to 191 px: the phone came
about 13 % closer while held "still", a change of scale rather than a shift.

## Result

Moving video, against the mean step between the same frames by phase
correlation:

| | Windows | Valid |
|---|---:|---:|
| Even steps | 91 | 83 |
| Uneven steps (start of the recording, changing motion) | 18 | 5 |

| 83 valid windows with even steps | Median | 95th percentile | Max |
|---|---:|---:|---:|
| Speed, \|est/step - 1\| | 0.19 % | 0.60 % | 1.39 % |
| Direction | 0.25 deg | 3.63 deg | 5.22 deg |

The speed is right: 7.3 to 15.4 px per frame, 35 to 75 mm/s (median 51 mm/s). The direction is
not: 10 of the 83 windows are more than 2 deg off, and in 9 of them the
estimate lies closer to the direction of the grid lines (0 deg) than the step
does. The phone moved at -6 to 11 deg to the lines; the estimate snaps towards
them, by up to 5.2 deg.

Still video (`static` 1, all must come out not valid): 3 of 89 windows are
marked valid, at 4.1 to 4.4 px per frame along the grid lines (0 and 180 deg),
while their frames moved by 0.03 to 0.11 px (one-sided 95 % upper bound of the
rate: 8.5 %). On the speckled desk top no recording of the phone lying still
exists yet, so this rate belongs to the grid, not to a random texture.

## What this shows

A regular pattern pulls the direction towards its own axes and lets a still
image pass as motion now and then. The forward model predicted wrong results
marked valid for periods of 64 and 128 px, with far larger errors than here
(the period here is about 125 px; the paper has a fine texture between its
thin lines, and the motion ran nearly along them).
The speed along the motion stayed within 1.4 % of the step check. A
measurement surface must therefore have a random texture
([test surface](../../docs/acquisition.md#test-surface)); graph paper is
useful only for the scale. The mm/s values rest on the 5 mm squares and on
the frame times of the video; there is no independent reference speed.
