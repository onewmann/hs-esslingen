# Status and how to resume

Saved on 11 October 2026, updated on 14 October 2026. This folder (`docs/dev/`) is a working note for
the rework and should be removed before the repository is published as
2D-FFT-Code. The plan this work follows is in
[restructure_plan.md](restructure_plan.md).

## Where things are

Work happens on branch `claude/pensive-gates-pzcxla` of `onewmann/hs-esslingen`.

On GitHub, the first seven commits after `8f0aec4` (`c3f77ab` to `ebe9cb6`)
still name "Claude" as author and carry `Co-Authored-By`/`Claude-Session`
lines; replacing them needs a force push, which the session's permission check
refused. Every later commit is authored by Oliver Neumann. Since 16 October
the session's clone follows the GitHub branch again (the earlier locally
rewritten history was only a copy with the same trees). publish.sh does the
cleaning for publishing (authors, trailers, tag), so no clean local history is
needed. Commit signing is off in the clone (`git config commit.gpgsign
false`); new commits are unsigned.

## Done

| Commit on GitHub | What |
|---|---|
| `4db98bc` to `4754d9c` | README: thesis credit, AI disclosure from commit `8f0aec4` on |
| `8fde67a` | Estimator safer on real multi-flash images: options `exposures`, `duty`, `roi`, `clip_outliers`; strict option checks; `pulsed_from_stack`; isotropy check (defocus, square repeating textures); threshold no longer lowered by padding beyond 2 |
| `66dc4ef` | CLI: `--L`/`--exposures` required, MATLAB v7.3 files, 16-bit and float images, frame axis of stacks, 512 px default ROI, JPEG warning |
| `5d2f735` | Live scripts: multi-flash mode, exposure set and read back, duty passed on, interval spread check, missing timestamps, failed grabs, pixel formats |
| `33ebf27` | Python 3.10+, matplotlib 3.6+, CI on 3.10 and 3.13, run_tests prints file and line, `--quick` no longer overwrites results, README updated |
| `1b3f653` | Recordings moved to `data/recordings`, archived 2025 scripts and the Basler driver removed from the tree (still at the tag) |
| `7706f29` | MATLAB code is now the package `fftvel` (`fftvel.estimate_velocity`, helpers in `+fftvel/+internal`) |
| `fcba12c` | `benchmark/` renamed to `experiments/`, one script per study, `fftvel.angle_difference`, version 0.1.0 |

Hashes are those on GitHub; the clean local history has other hashes for the
same changes.

The failure review of 11 October (46 findings, plus 2 from a final critic)
is closed: every finding is fixed or documented. Still open or only
documented:

- Defaults keep `duty = 0` and `clip_outliers = false`, so long flashes and
  glints still need the option (README Limitations says so).
- Rectangular repeating textures can still pass as motion; about 2 of 236
  simulated defocused images still pass.
- MATLAB live script uses host times from `getdata`, not camera timestamps
  (chunk data, R2024a+, not implemented).
- Nothing has run on camera hardware yet. The Python camera code was checked
  with the pylon emulator. (MATLAB: green in CI since run 12, see 14 October.)
- Host timing for webcams: alternating intervals pass when checked against
  the nominal period.

## Studies and tools added on 11 October (commits `93b60ec`, `38bc248`)

- `matlab/+fftvel/+internal/constants.m` and `tests/constants.json`: 28
  algorithm constants checked against one file in both languages.
- Keyword-only options in Python, smaller `fftvel.__all__`,
  `fftvel.to_physical(est, gsd, period)` in both languages, CLI `--gsd` and
  `--period`.
- `experiments/timing.py`: 256 px takes 0.84 s in Python and 1.30 s in
  Octave per call, 512 px 4.2 s and 6.2 s (4 cores, `results/timing.csv`).
- `experiments/real_measurements.py` and `data/real/README.md`: evaluation of
  the owner's multi-flash recordings from a `meta.csv`, with a synthetic
  example in `data/real/example/`.
- `experiments/calibrate_threshold.py`: largest k on the 1008 null images is
  8.41 (the docstring still says 8.6); 0 of 2268 simulated null images valid;
  held out 4 of 1512 valid, all crops of the static room scene vid1.
- `experiments/baselines.py`: share correct over 840 images, random / camera
  texture: fftvel 80 / 63 %, peak spacing 69 / 56 %, cepstrum 70 / 60 %,
  autocorrelation 72 / 55 %, phase correlation with two frames 95 / 88 %.
- `experiments/model_violations.py`: independent forward model; wrong but
  valid results appear with pulse width 0.25 to 0.5 without `duty`, step
  jitter of 10 %, continuous light during most of a long exposure, a strong
  static edge, periodic textures with a period below about 160 px, glints.
- run.py parts `calibration baselines violations timing` and `--all`; CI has
  a MATLAB job, the Octave demos and a smoke run.
- 284 Python tests and 9 Octave test files pass.

Follow-ups done on 12 October:

- Blur check: fields `quality_blur`, `blur_length`, `blur_like` and the option
  `reject_blur`. Not part of `valid` by default: below about 2.2 px per flash
  a blurred single exposure and a slow multi-flash image cannot be told
  apart, and as a default it would cost 3.6 % of correct results on the noise
  grid. It flags 41 of the 47 valid single frames of the moving surface.
- Field `saturated` with a warning above 1 %; `clip_outliers` also clears a
  3 px rim around each outlier.
- `phase_correlation` refines its peak by upsampling: median error 0.002 px
  instead of 0.12 px.
- All studies rerun; default results unchanged.
- New README (status first, no chapter structure), `docs/method.md`,
  `docs/results.md`, `docs/hardware.md`, `docs/history.md`, `hero.png`,
  `CITATION.cff`, `CONTRIBUTING.md`. Numbers checked against the CSV files,
  style checked (British spelling, no dashes, no bold in sentences).
- 306 Python tests and 9 Octave test files pass.
- Final review (commit `4a15cd2`): test inputs that depended on Octave's
  random stream are stored, so the first MATLAB run checks the same images;
  float comparisons in the baseline tests tolerate other CPUs and the
  minimum numpy; several docstrings and doc statements corrected.

Testing without a flash (added on 12 October): `experiments/real_measurements.py`
has a mode `average` (meta.csv columns `mode`, `frames`, `first_frame`,
`frame_period_s`): K short exposures at a fixed frame rate are averaged into
one pulsed image, the frame intervals (from timestamps) and the steps
between frames (phase correlation) are checked, and uneven windows are
flagged. `live_camera.py --save DIR` and `live_basler(..., 'SaveDir', DIR)`
store bursts with timestamps and append lines to `meta.csv`. Example in
`data/real/example_average/`, instructions in `data/real/README.md`
("Testing without a flash"). Results from this mode count as emulated
multi-flash.

Naming and structure review (13 October): `data/real` is now
`data/measurements` (`synthetic_flash/`, `synthetic_average/`),
`docs/hardware.md` is `docs/acquisition.md`, `experiments/testsequences.py`,
`docs/figures/violations.png`, new module `fftvel.io`, `fftvel/phasecorr.py`
(was `reference.py`), console command `fftvel`, MATLAB `Contents.m`,
live_basler options `'ExposureUs'` and `'UpdateInterval'`, README reordered
(hero first, one "Who did what", Quick start, Usage), data under CC BY 4.0.
The repository name stays 2D-FFT-Code. 344 Python tests and 10 Octave test
files pass.

Reframe (14 October): the project is now presented around averaged camera
frames (K frames at a fixed frame rate, averaged); the single multi-flash
image is the extension. The README states openly that phase correlation of
separate frames is more accurate and gives the sign. The owner's paragraph
in "Who did what" describes the AI-assisted way of working.

Audit and fixes (14 October, afternoon). A verified sweep (four auditors,
one skeptic each, a completeness critic) listed what was still open before
publishing. Done in this session:

- CI was red since run 17 only because the Octave demo step lacked a font;
  fixed. Run 24 is green in all five jobs: Python 3.10, 3.13 and 3.14,
  Octave (tests and demos) and MATLAB (tests, demos and the code analyser
  over all .m files, live_basler.m included). The actions use Node.js 24.
- Docs no longer say that the code has not run in MATLAB. `vid1` is now
  described as a room scene filmed with a moving camera (it is not static;
  only single frames of it are used, so no number changed). README: "also
  analyses" for multi-flash images, a status sentence in the intro, a header
  for the status table, caveats on whole-pixel test sequences, on the null
  images the threshold constant was checked on, and on averaging inside the
  camera.
- real_measurements.py: warning for a frame burst read in flash mode,
  statistics against the mean step (median, p95, max, number off with
  bound), column `step_vs_set`, a warning for a stopped axis labelled moving,
  valid windows that agree with the step but not with the set motion are
  blamed on the reference, and a figure against the mean step without a set
  speed. live_camera.py: `--count`, `--static`, compressed files, a message
  for skipped bursts. docs/acquisition.md: test surface rules and "First
  recording, step by step" with counts, exclusion rule and success criteria.
  Package metadata has project URLs.
- publish.sh also removes the GPS data of the 2025 photo, the Windows user
  folder in one 2025 script and the early README sentences about the time
  before the rework, and checks every blob for them (check_history.py).
  Dry run on 14 October from the local branch: 43 commits, all checks pass,
  final tree equal to the branch minus docs/dev, photo pixels unchanged.
- 367 Python tests and 10 Octave test files pass.

Phone video (16 October). The owner recorded the desk top with a phone: three
videos at 30 fps (76 to 93 px per frame, motion blur; 51 of 268 windows of 5
frames valid; not committed) and one at 120 fps, which is now the session
data/measurements/2026-10-16_phone_120fps (video stream unchanged, sound and
metadata removed; extract_frames.py with OpenCV; 68 of 68 windows with even
steps valid, median 0.21 % and 0.13 deg against the mean step). README Status,
results.md, acquisition.md (section Phone videos), CITATION.cff and
data/README.md updated; CI extracts and evaluates the session. Commit signing
is off in this clone (git config commit.gpgsign false), as the owner asked.
Still missing: a recording of the phone lying still, and a physical reference
(GSD from a ruler in the image, a known speed).

Graph paper and sound tracks (17 October). Two phone videos of 5 mm graph
paper at 120 fps (moved, held still) are the session
data/measurements/2026-10-17_phone_grid: speed within 1.4 % of the step check,
direction pulled towards the grid lines (10 of 83 windows more than 2 deg off),
3 of 89 still windows valid; the scale of each window comes from the grid
(0.036 to 0.041 mm/px, the hand changed height). All videos in the tree have
no sound track; publish.sh replaces the two 2025 phone videos at the tag by
copies without sound and metadata (clean_videos/), and check_history.py fails
on any video with a sound track. The first two commits get new hashes
(5f2dad8, bae6504). Still missing: the phone lying still on the desk top (a
random texture) and an independent reference speed.

## Next steps

1. Record the test without a flash: follow docs/acquisition.md, section
   "First recording, step by step", and data/measurements/README.md, section
   "Averaged frames (mode average)". Confirm the success criteria there
   (90 % valid, 5 % / 2 deg, none of the static bursts valid) before the
   session. Evaluate with experiments/real_measurements.py; then Claude
   updates README Status, acquisition.md and results.md.
2. Publishing:
   - Create the empty public repository `onewmann/2D-FFT-Code` on GitHub (no
     README, licence or .gitignore).
   - In Git Bash: `pip install git-filter-repo`, clone hs-esslingen, check out
     `claude/pensive-gates-pzcxla`, run `bash docs/dev/publish/publish.sh`. It
     must end with "clean history: ... (N commits, tag student-project-2025)";
     N was 43 on 14 October and grows with every later commit outside
     docs/dev. Then `bash docs/dev/publish/publish.sh --push`.
   - Check that all five CI jobs are green on 2D-FFT-Code, that the badge
     shows passing, and that the tag link in README "Who did what" opens.
   - Settings of 2D-FFT-Code: About "Surface velocity from averaged camera
     frames or a single multi-flash image. 2-D FFT and Radon transform, Python
     and MATLAB/Octave."; topics velocity-measurement, machine-vision,
     frame-averaging, image-processing, fft, radon-transform, python, matlab,
     octave; social preview `docs/dev/publish/social_preview.png` (1280 x
     640, in the clone used for publishing); Wiki and Projects off.
   - Release: on release day set `date-released` in CITATION.cff, tag v0.1.0
     and write release notes that state the status honestly (simulation and
     synthetic sequences; real recordings pending).
   - Then make hs-esslingen private (Settings > General > Danger Zone >
     Change visibility), after `git clone --mirror` as a backup. No README
     pointing to the new repository. Its MATLAB CI job will fail once it is
     private; that does not matter any more. Running publish.sh later then
     needs your GitHub login or `SRC=<local clone>`.
3. Check your public GitHub profile (bio, company, location, pinned
   repositories) for the university's name.
4. Future commits on 2D-FFT-Code: author and committer Oliver Neumann with
   the noreply address, no Co-Authored-By or Claude-Session lines, no
   `claude/` branch names on the public repository. Cloud sessions add these
   lines by default; say so at the start of each session (or in a CLAUDE.md,
   which would be public). For local Claude Code, `"attribution": {"commit":
   "", "pr": ""}` in ~/.claude/settings.json turns them off.

## Decisions

| | Decision |
|---|---|
| D1 data | May be published: data produced by the owner's code (test sequences) and the videos vid1 to vid4 |
| D2 history | Cleaned at publishing (no university name, owner as author, no AI trailers, no photo metadata or user folder); `docs/dev/publish/publish.sh` does it |
| D3 old repository | hs-esslingen is made private after publishing (14 October), no pointer README |
| D4 attribution | No Claude author and no AI trailers; the README section "Who did what" is the disclosure. "Claude wrote most of the new code" stays: the owner wrote and developed parts together with AI |
| D5 real experiment | Frame bursts at a fixed frame rate first (docs/acquisition.md, "First recording, step by step"), multi-flash images later; format in data/measurements/README.md |
| D6 Python | 3.10 and later; CI on 3.10, 3.13 and 3.14 |
| D7 names | Hochschule Esslingen is not named anywhere. Universität Karlsruhe (TH) stays in the citation of Jan Horn's thesis |
| D8 AI before the rework | No statement about the 2025 scripts either way; the early README sentences "my own work without AI" are removed from the published history |

Note on D7: the history still contains `2D-FFT/README.md` (a GitLab template
with a gitlab.hs-esslingen.de link) at the tag, and the archive repository
is called `hs-esslingen`; publish.sh removes the first, D3 handles the second.

## The owner's part

- Record and evaluate the frame bursts (next step 1), later multi-flash
  images. Record the stopped surface too; it must come out not valid.
- Check that the project terms allow publishing the 2025 scripts and
  recordings under MIT and CC BY 4.0, and that whoever is in the foreground of
  `vid1` (a dark shape close to the camera, no face) agrees, if it is a
  person other than you.
- Interview preparation:
  - one concrete 2025 error with numbers (docs/history.md): the fixed
    `findpeaks` prominence of 0.05 on summed Radon projections gave 32 to
    78 px per step whatever the true speed (median error 733 % in the Python
    port, coarse), and `imresize(I, [256 256])` of a 1936 x 1216 frame moved
    45 deg to 57.9 deg;
  - why averaging K frames gives the multi-flash image, why the sign of the
    motion is missing (phase correlation of two frames gives it), and why
    phase correlation is more accurate when the frames exist;
  - the validity threshold: k is a calibrated score, not a significance
    level, and it has not been checked on real bursts yet;
  - which parts you wrote or shaped yourself, why all published commits
    carry your name (D4), and two or three cases where you or a review
    rejected generated output (for example 7f870e8, the reframe 5ab74ef);
  - a walk through `estimate_velocity` and `fit_speed` in
    python/src/fftvel/core.py.
- Derive the image model, the speed template, the window autocorrelation and
  the threshold on paper, and rewrite those sections of docs/method.md in
  your own words.
- Read the relevant chapters of Horn's thesis and add section or page numbers
  to the [thesis] tags in README and method.md; move anything the thesis
  already has from [rework] to [thesis].
- Optional: run the 2025 MATLAB chain (finaltest1.m, lines 92 to 126) on the
  36 before/after images in MATLAB with the Image Processing and Signal
  Processing toolboxes, to confirm the port behind the 733 %.

## Resume in a new session

```
git clone https://github.com/onewmann/hs-esslingen && cd hs-esslingen
git checkout claude/pensive-gates-pzcxla
python -m pip install --upgrade pip
pip install -e "./python[test]"
python -m pytest -q tests/python
cd tests/matlab && octave --no-gui --eval run_tests
```

Then ask Claude to read `docs/dev/status.md` and continue with the next steps.
Say at the start that commits must have no Claude author or AI trailers.
