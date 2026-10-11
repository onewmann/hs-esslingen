# Status and how to resume

Saved on 11 October 2026. This folder (`docs/dev/`) is a working note for
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

## In progress when this was saved

A workflow was adding the studies and tools below. Its files are committed
as they were at that moment, under the commit "Save work in progress".
268 Python tests passed on that state; the reviews were still running.

- `matlab/+fftvel/+internal/constants.m` and `tests/constants.json`: all
  algorithm constants checked against one file in both languages.
- Keyword-only options in Python, smaller `fftvel.__all__`,
  `fftvel.to_physical(est, gsd, period)` in both languages, CLI `--gsd` and
  `--period`.
- `experiments/timing.py` (run times, `results/timing.csv`,
  `results/environment.json`).
- `experiments/real_measurements.py` and `data/real/README.md`: evaluation of
  the owner's own multi-flash recordings from a `meta.csv`, with a synthetic
  example in `data/real/example/`.
- `experiments/calibrate_threshold.py`: reproducible threshold calibration
  (1008 null images), held-out test, real null test on single Basler frames,
  power for constants 8 to 10, GEBA comparison, half-speed check
  (`results/threshold_*.csv`, `docs/figures/threshold.png`).
- `experiments/baselines.py`: cepstrum, autocorrelation, corrected peak
  spacing, phase correlation as a two-frame reference
  (`results/baselines.csv`, `docs/figures/baselines.png`).
- `experiments/model_violations.py`: independent forward model with pulse
  width, step and energy jitter, shot and read noise, quantisation, ambient
  light, defocus, static edge, periodic texture, glints, wrong K
  (`results/violations.csv`, `docs/figures/failures.png`).

Still to do for this step: finish the reviews of the three studies, register
the new studies in `experiments/run.py` (`--only calibration baselines
violations timing`, `--all`), extend CI (MATLAB job with
`matlab-actions/setup-matlab@v2`, Octave demos, a short experiments smoke
run and the real-data example), run both suites.

## Next steps

1. Finish the step above and commit it in logical commits.
2. Rewrite README.md following section 2 of the plan: no chapter structure,
   a status table first ("works in simulation; real multi-flash images
   not shown yet"), who did what, hero figure, try it, why one image, how it
   works with [thesis]/[rework] tags, results ordered by independence, where
   it fails, use it, development, how this was built, reference.
3. Docs pages: `docs/method.md`, `docs/results.md`, `docs/hardware.md`,
   `docs/history.md` (index of the 2025 scripts by function). `hero.png`,
   retitled figures, `CITATION.cff`, `CONTRIBUTING.md`.
   Style for everything published: British spelling, no em dashes, no bold
   inside sentences, no marketing words.
4. Final review: every number in the docs against `experiments/results`,
   both suites, a check that no text names Hochschule Esslingen.
5. Publishing: create the empty repository `onewmann/2D-FFT-Code`, add it to
   the Claude session with push access (or push it yourself), push the clean
   history as `main` plus the tag. Then set About, topics and the social
   preview, make release v0.1.0, add a root README to hs-esslingen pointing
   to the new repository and archive hs-esslingen read-only.

## Decisions

| | Decision |
|---|---|
| D1 data | May be published: data produced by the owner's code (test sequences) and the videos vid1 to vid4 |
| D2 history | Push the clean history unchanged into 2D-FFT-Code (recommended; not yet confirmed) |
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
