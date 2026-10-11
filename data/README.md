# Data

Both files come from the original study project and are small enough for the
repository.

## `testsequences.mat`

Four synthetic sequences of 256 x 256 px, stored as `uint8` (divide by 255).
Each frame is the previous frame shifted down and to the right with zero fill,
so the black border grows from frame to frame. Analyse only the area that
holds image content in every frame of a window (see `tests/matlab/test_testsequences.m`).

| Variable | Frames | Motion per frame (x right, y down) |
|---|---|---|
| `model` | 20 | (5, 5) px, constant |
| `tyre` | 20 | (4, 4) px, constant |
| `vehicle` | 30 | about (5, 5) px, single steps with vx and vy between 1.9 and 7.9 px |
| `accelerating` | 30 | vy = 3 px, vx varies smoothly between 2 and 8 px |

Converted from `Testsequenzen1.mat` of the original project (variables
`modelSequence`, `reifenSequence`, `fahrzeugSequence`, `simSequence`; values in
[0, 1] scaled to 0..255).

## `surface_texture.png`

Central 1024 x 1024 px of the first frame of recording `vid2.mat`, taken with a
Basler camera (Mono8, 1936 x 1216 px) looking at a moving surface. Used as a
realistic texture for simulations: `simulate_pulsed_image(..., 'Texture', tex)`
moves a fixed window over it with exact, known shifts.

The full recordings (`vid1.mat` to `vid4.mat`, 10 frames each) stay in the
archive of the original project (`2D-FFT/data` in onewmann/hs-esslingen).
`vid1` shows a static room scene; `vid2` to `vid4` show the moving surface.
Each file holds the variable `recording`, 1216 x 1936 x 1 x 10 `uint8`, as
`getdata` of the Image Acquisition Toolbox returns it. `pulsed_from_stack`
(MATLAB and Python) and `python -m fftvel estimate vid2.mat --stack --L 4`
read it as it is; with 10 frames, an even number, they need `L` or
`Exposures` to know how many frames to average. The frame-to-frame motion in
`vid2` to `vid4` is not uniform, which
`benchmark/run_benchmark.py --only recordings --recordings PATH/TO/2D-FFT/data`
analyses.
