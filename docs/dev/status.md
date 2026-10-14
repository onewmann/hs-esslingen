# Status and how to resume

Saved on 11 October 2026, updated on 13 October 2026. This folder (`docs/dev/`) is a working note for
the rework and should be removed before the repository is published as
2D-FFT-Code. The plan this work follows is in
[restructure_plan.md](restructure_plan.md).

## Where things are

Work happens on branch `claude/pensive-gates-pzcxla` of `onewmann/hs-esslingen`.

There are two versions of the same history:

- On GitHub, the first eight commits after `8f0aec4` (`c3f77ab` to `ebe9cb6`)
  still name "Claude" as author and carry `Co-Authored-By`/`Claude-Session`
  lines. Replacing them needs a force push, which the session's permission
  check refused. Every later commit on GitHub is authored by Oliver Neumann.
- The clean history (all commits after `8f0aec4` authored by Oliver Neumann,
  no AI trailers) exists only in the local clone of the session. It can be
  rebuilt at any time from the GitHub branch. The trees are identical; this
  was checked on 11 October:

```
git clone https://github.com/onewmann/hs-esslingen && cd hs-esslingen
git checkout -b clean origin/claude/pensive-gates-pzcxla
git filter-branch -f \
  --env-filter 'export GIT_AUTHOR_NAME="Oliver Neumann" GIT_AUTHOR_EMAIL="149961666+onewmann@users.noreply.github.com" GIT_COMMITTER_NAME="Oliver Neumann" GIT_COMMITTER_EMAIL="149961666+onewmann@users.noreply.github.com"' \
  --msg-filter 'grep -v -e "^Co-Authored-By: Claude" -e "^Claude-Session:"' \
  8f0aec4..HEAD
git tag -a student-project-2025 8f0aec4 -m "Student project scripts as uploaded. The AI-assisted rework starts after this commit."
```

The tag `student-project-2025` (pointing to `8f0aec4`) exists only locally;
pushing tags was refused as well, hence the last command above.

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
- Nothing has run in real MATLAB or on camera hardware yet. The Python
  camera code was checked with the pylon emulator.
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

## Next steps

1. Done: the README section "Who did what" now says that the owner noticed
   errors in the 2025 results and rebuilt the project with Claude. Be ready
   to name a concrete error in an interview.
2. Record the test without a flash (`data/measurements/README.md`, section
   "Testing without a flash") and later multi-flash images; evaluate them
   with `experiments/real_measurements.py`; update Status and Results.
3. Publishing:
   - Create the empty repository `onewmann/2D-FFT-Code` on GitHub.
   - Run `bash docs/dev/publish/publish.sh` (Git Bash on Windows, needs
     `pip install git-filter-repo`). It clones the working branch, removes
     `docs/dev/` and the 2025 GitLab template README from every commit,
     replaces the remaining mentions of the university, sets every author
     to Oliver Neumann, removes the AI trailers and checks the result.
     Tested on 13 October from both the local and the GitHub branch: same
     final tree, 29 commits, nothing left to clean.
   - Then `bash docs/dev/publish/publish.sh --push`, or add the new
     repository to a Claude session and let it push.
   - Set About ("Surface velocity from a single multi-flash image. 2-D FFT
     and Radon transform, Python and MATLAB/Octave"), topics
     (velocity-measurement, machine-vision, fft, radon-transform, python,
     matlab, octave), `docs/figures/hero.png` as social preview, set
     `date-released` in CITATION.cff to the release day, tag v0.1.0 and make
     a release. Check that the MATLAB CI job is green, then update the
     MATLAB row of the README Status table.
   - Afterwards add a root README to hs-esslingen pointing to the new
     repository and archive hs-esslingen read-only.

## Decisions

| | Decision |
|---|---|
| D1 data | May be published: data produced by the owner's code (test sequences) and the videos vid1 to vid4 |
| D2 history | Cleaned at publishing (no university name, owner as author, no AI trailers); `docs/dev/publish/publish.sh` does it |
| D3 old branch | Keep the branch in hs-esslingen and archive the repository (recommended; not yet confirmed) |
| D4 attribution | No Claude author and no AI trailers; the README section "How this was built" is the disclosure |
| D5 real experiment | The owner records real multi-flash images in the coming days (`data/real/README.md` describes the format) |
| D6 Python | 3.10 and later (done) |
| D7 names | Hochschule Esslingen is not named anywhere. Universität Karlsruhe (TH) stays in the citation of Jan Horn's thesis |

Note on D7: the history still contains `2D-FFT/README.md` (a GitLab template
with a gitlab.hs-esslingen.de link) at the tag, and the archive repository
is called `hs-esslingen`. Renaming that repository on GitHub would remove
the name from links; GitHub redirects the old URL.

## The owner's part

- Record and evaluate real multi-flash images (`experiments/real_measurements.py`).
  Give the exact number of flashes and, unless the flashes are very short,
  the duty. Record the stopped surface too; it must come out not valid.
- Derive the speed template, the window autocorrelation and the threshold on
  paper and rewrite `docs/method.md` in your own words.
- Read the related work before anything is cited.
- After publishing, check that the MATLAB CI job is green.

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
