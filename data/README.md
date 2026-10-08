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
| `vehicle` | 30 | about (5, 5) px, single steps between 2.2 and 8 px |
| `accelerating` | 30 | vy = 3 px, vx varies smoothly between 1.8 and 8.2 px |

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
Their frame-to-frame motion is not uniform, which
`benchmark/run_benchmark.py --recordings PATH/TO/2D-FFT/data` analyses.
