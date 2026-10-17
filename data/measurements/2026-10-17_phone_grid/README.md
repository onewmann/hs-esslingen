# Graph paper: a periodic texture

Two phone videos of graph paper with 5 mm squares, recorded in slow motion at
120 frames per second: one with the phone moved by hand, one with the phone
held still over the paper. A regular grid breaks assumption 1 of
[method.md](../../../docs/method.md#assumptions) (a texture that does not repeat),
so this session shows how the estimator behaves on a real periodic texture.
The grid also gives the scale: the squares turn pixels into millimetres. Data
and figure are under CC BY 4.0
([licence of the data](../../README.md#licence-of-the-data)).

| File | What |
|---|---|
| `grid_moving_120fps.mp4` | phone moved by hand, 546 frames, 4.7 s; the video stream as recorded, without sound track and file metadata |
| `grid_still_120fps.mp4` | phone held still, 447 frames, 3.9 s; the same |
| `extract_frames.py` | writes `frames_moving.npz` and `frames_still.npz` (the central 512 x 512 px square of every frame, with the timestamps), `meta.csv` and `grid_lines.csv` |
| `meta.csv` | 109 windows of 5 frames of the moving video, each with the ground sample distance measured on the grid, and 89 of the still one (`static` 1) |
| `grid_lines.csv` | for every window of both videos: the period of the grid in px, the ground sample distance and the direction of the grid lines in the image |
| `results.csv`, `results.png` | the output of `experiments/real_measurements.py` |
| `directions.py` | compares the estimated directions with the grid lines and with the steps |

```bash
python -m pip install -e "./python[plot,data,video]"
python data/measurements/2026-10-17_phone_grid/extract_frames.py
python experiments/real_measurements.py data/measurements/2026-10-17_phone_grid/meta.csv --figure data/measurements/2026-10-17_phone_grid/results.png
python data/measurements/2026-10-17_phone_grid/directions.py
```

The committed `results.csv` was made with OpenCV 5.0.0. Other versions
decode the video slightly differently (grey values up to 2 levels apart), and
the still video is sensitive to that (see below).

## Scale

A square is about 130 px wide in the moving video (122 to 138 px), measured
on the 1024 px field around the analysed square; the phone's height changed as
the hand moved. Each window of the moving video therefore gets its own ground
sample distance, 5 mm divided by the square's width, between 0.036 and
0.041 mm/px. Every window gets one except the last two, where the paper's edge
enters the field.

## Moving video

Against the mean step between the same frames by phase correlation:

| | Windows | Valid |
|---|---:|---:|
| Even steps | 91 | 83 |
| Uneven steps (start of the recording, changing motion, one dropped frame) | 18 | 5 |

| 83 valid windows with even steps | Median | 95th percentile | Max |
|---|---:|---:|---:|
| Speed, \|est/step - 1\| | 0.19 % | 0.60 % | 1.39 % |
| Direction | 0.25 deg | 3.63 deg | 5.22 deg |

The speed is right: 7.3 to 15.4 px per frame, 35 to 75 mm/s (median
51 mm/s). The direction is not. The hand turned the phone a little, so the
grid lines lie at -0.1 to 4.1 deg in the image (`grid_lines.csv`), and the
phone moved at -9.9 to 9.2 deg to them. 10 of the 83 windows are more than
2 deg off the step, by up to 5.2 deg (frames 31 to 35); in all 10 the
estimate lies within 0.36 deg of the grid lines, and in all 10 the motion ran
within 5.3 deg of them. Most windows with motion that close to the lines are
not pulled that far, but over all 83 windows the estimate lies closer to the
lines than the step does in 62.

The window of frames 226 to 230 holds the one dropped frame of this video (a
gap of 25 ms, three frame intervals). It is valid and 34 % below the mean step,
but only 1.0 % below the median step, the regular step between the other
frames, as in the [desk-top video](../2026-10-16_phone_120fps/README.md#result).
The step check flags it, and it is left out of the numbers above.

## Still video

All windows of the still video (`static` 1) must come out not valid. 3 of 89
are marked valid (frames 306 to 320): 4.1 to 4.4 px per frame along the grid
lines (0 and 180 deg), while their frames moved by 0.03 to 0.11 px. The
one-sided 95 % upper bound of the rate is 8.5 %.

This number is fragile. The quality of the three windows is 0.60 to 0.63,
just above the threshold of 0.60. With OpenCV 4.12 instead of 5.0.0, 13 of 89
windows come out valid (upper bound 22 %), all with a quality of 0.62 to
0.70. On the moving video the two versions agree closely: 81 instead of 83
valid windows with even steps, largest speed deviation 1.25 % instead of
1.39 %, again 10 windows more than 2 deg off, 9 of them the same.

The phone was not quite still. In 9 windows (frames 181 to 200, 396 to 410
and 426 to 435) the hand jolted it by 1.5 to 4.5 px per frame; the step check
flags them, and none is valid. Around frames 160 to 200 the squares grow from
168 to 191 px, as the phone came about 12 % closer to the paper, and stay at
184 to 191 px after that.

On the speckled desk top no recording of the phone lying still exists yet,
so this rate belongs to the grid, not to a random texture.

## What this shows

A regular pattern pulls the direction towards its own axes and lets a still
image pass as motion now and then; how often depends on details as small as
the decoding of the video, because the quality of a still grid sits at the
threshold. The forward model predicted wrong results marked valid for
periods of 64 and 128 px, with far larger errors than here (the period here is
about 130 px; the paper has a fine texture between its thin lines, and the
motion ran nearly along them). The speed along the motion stayed within 1.4 %
of the step check. A measurement surface must therefore have a random texture
([test surface](../../../docs/acquisition.md#test-surface)); graph paper is
useful only for the scale. The mm/s values rest on the 5 mm squares and on
the frame times of the video; there is no independent reference speed.
