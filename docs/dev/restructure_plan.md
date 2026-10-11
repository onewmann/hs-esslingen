# 2D-FFT-Code: the restructuring plan, merged from three proposals

Written on 9 October 2026, before the work started. The current state and the
decisions taken since are in [status.md](status.md); where the two differ,
status.md is right.

## 0. The exact point where AI use starts

I checked this in the local repo (branch `claude/pensive-gates-pzcxla`, HEAD `73821f1`).

| Point | Commit | Time |
|---|---|---|
| 2025 files uploaded, part 1 | `54b90cc` "2D-FFT Studienprojekt importieren" | 8 Oct 2026, 10:41 CEST |
| 2025 files uploaded, part 2: **AI use starts after this commit** | `8f0aec4` "Aktuelle Fassung der Skripte ergänzen (Radon-GUI, finaltest1, Filterbank)" | 8 Oct 2026, 12:02 CEST |
| **First use of AI**: code review of the 2025 scripts | no commit | between 12:02 and 23:23 CEST |
| **First AI-assisted commit** | `680a45f` "Rewrite the 2D-FFT velocity estimator for MATLAB, Octave and Python" | 8 Oct 2026, 23:23 CEST (21:23 UTC) |
| Every later commit, up to `73821f1`, and all future restructure commits | `0d90144` … `73821f1` | from 8 Oct 2026 on |

Checks behind the table:
- `git diff 8f0aec4 HEAD -- 2D-FFT Matlab_BaslerCamDriver` is empty. The 2025 files have not changed since they were uploaded.
- `54b90cc` also contains `Matlab_BaslerCamDriver/`. That folder is third-party code (MIT, Hannes Badertscher, HSR 2015) and must be named as such wherever the upload is described.
- `gh api` reports onewmann/hs-esslingen as **public**, so two things are already visible:
  - The same AI commits exist there with "Claude" as author: `origin/claude/pensive-gates-pzcxla`, `c3f77ab`…`ebe9cb6`.
  - The local copies `680a45f`…`9d5d141` (plus `e9c70f9`/`4a62d6c`) have the student as author. Their trees are identical: `git diff` is empty for all 8 pairs.
- The README names the commit, not a date: both sides of the boundary fall on 8 October 2026.

**Text for the first screen** (adjusted to the student's decision: no claim about the time before AI):
> **Who did what.** Method: Jan Horn, doctoral thesis, Universität Karlsruhe (TH), 2006/2007. The 2025 student project files are uploaded unchanged up to and including commit [8f0aec4](https://github.com/onewmann/2D-FFT-Code/commit/8f0aec4) (tag `student-project-2025`); the Basler driver in that upload is third-party code. **Everything after 8f0aec4 was done with the assistance of AI (Claude, Anthropic):** first a code review of those scripts, then every commit from [680a45f](https://github.com/onewmann/2D-FFT-Code/commit/680a45f) "Rewrite the 2D-FFT velocity estimator…" on. [All AI-assisted changes as one diff](https://github.com/onewmann/2D-FFT-Code/compare/student-project-2025...main).

The full section "How this was built" adds the table above and, if the student wants, a paragraph "My part in the rework" with only true statements. No commit trailers (student's decision); the README section is the disclosure.

These links only work if the new repo contains the same hashes; see C1 below.

## 1. Conflicts between the proposals and how they are resolved

| # | Conflict | Decision | Reason |
|---|---|---|---|
| C1 | History: filter-repo into a fresh repo (reader-first) vs. push the existing history unchanged (engineering-first, science-first) | **Push unchanged** | Both `8f0aec4` and `680a45f` then exist in 2D-FFT-Code, and the compare link shows exactly the AI part. The history was already rewritten once (Claude to student); a second rewrite would create a third set of hashes and look like laundering. filter-repo also protects no data, because hs-esslingen is already public with the same files. Cost: about 110 MB of 2025 data in the history, so the README mentions `git clone --filter=blob:none`. Use reader-first's filter-repo only as a fallback if the supervisor refuses the data (D1). |
| C2 | 2D-FFT/ in the working tree: `git mv` to `original-2025/` (science-first) vs. remove it (the other two) | **Remove it; it stays reachable at the tag** | The student asked to get away from the chapter structure. `original-2025/` would keep about 50 chapter-named files and 112 MB in every checkout. At the tag the files are byte-identical by construction. Science-first's topic index moves to docs/history.md. |
| C3 | Where the README links the 2025 files | **`2D-FFT-Code/tree/student-project-2025/2D-FFT`** | One repo to check, and a permalink that never moves. The current link points to hs-esslingen `main`. |
| C4 | MATLAB helper folder: `+fftvel/private` (reader-first, science-first, critique) vs. `+fftvel/+internal` (engineering-first) | **`+internal` with fully qualified calls** | I retested this in Octave 8.4.0. A function in `+pk/private/` is "undefined" when called from `+pk/pub.m`. `pk.internal.h2` works. Using `private` would break the Octave CI. |
| C5 | Name of `benchmark/`: keep it, `experiments/` or `validation/` | **`experiments/`** | "benchmark" suggests speed. "validation" overclaims for numbers that come only from simulation, which is the critique's main point. |
| C6 | Commit trailers: rewrite all past messages, add from now on, or none | **None** (student's decision) | The README section and the compare link document the boundary. |
| C7 | Version: v0.1.0 or v1.0 | **0.1.0** (replaces the hard-coded `__version__ = "1.0.0"`) | Nothing has been validated on real data. |
| C8 | `pyproject.toml` at the root, or kept in `python/` | **Keep it in `python/`** | Keeps matlab/ and python/ symmetric. Moving it touches CI and the install commands for little gain. `python/README.md` shrinks to a pointer. `python/LICENSE` stays for the wheel. |
| C9 | Which docs pages | **method.md, results.md, hardware.md, history.md, related-work.md** | usage.md goes into the README "Use it" section and hardware.md. student-project-2025.md merges into history.md. |
| C10 | README numbers generated by `sync_docs` + `summary.json`, or written by hand | **By hand.** Optional later: a test that checks the README numbers against the CSVs | The generator is too heavy for a portfolio repo; a check gives the same safety. |
| C11 | Hero figure name | **`docs/figures/hero.png`**, also set as the GitHub social preview | |
| C12 | Limitations bullets, a failure table, or folded into "What is open" | **One section "Where it fails"**: a table with the silent failures first. The full assumptions list goes to method.md | The failure table is the most useful form for a machine-vision reader. |
| C13 | Related methods in the README | **A 3-row table** inside "Why one image", listing only methods that get a baseline. Everything else goes to related-work.md once the student has read the sources | No unread citations on the first page. |
| C14 | The archive repo hs-esslingen | **Do not touch anything inside 2D-FFT/ there.** Add a root README pointing to 2D-FFT-Code, keep the claude/* branch, archive the repo read-only | |
| C15 | Where the constants spec lives: `spec/constants.json` or elsewhere | **`tests/constants.json`** | Only the tests read it, so it needs no top-level folder. |

## 2. README outline, in final order (about 200 lines)

**First screen**
1. **Title:** "Surface velocity from a single multi-flash image".
   - Subtitle: "2-D FFT and Radon transform · Python and MATLAB/Octave".
   - One tests badge: a single workflow with python, octave and matlab jobs.
   - GitHub About: "Direction and speed of a moving textured surface from one multi-flash camera image". Topics: machine-vision, signal-processing, fft, radon-transform, matlab, octave, python.
2. **What it is for**, 3 plain lines:
   - A strobe flashes 2L+1 = 9 times during one exposure, so the surface appears 9 times, slightly shifted, in one image.
   - The shifts become evenly spaced lines in the 2-D spectrum. Their angle gives the direction, their spacing the distance moved per flash.
   - Use: contactless speed over ground and sideslip angle of a vehicle (Horn's application), or a conveyor or paper web. One image per measurement, no tracking between frames.
3. **Status**, a 3-row table:
   - **Works in simulation and on synthetic sequences.** Every number is labelled "in simulation".
   - **Not shown yet on a real multi-flash image.** The 2025 recordings break the model (unequal steps of 1.6–18 px, direction spread of about 30–40°). All 6 windows are rejected: quality 0.23–0.40, threshold 0.6.
   - **Tested in CI:** Python and GNU Octave 8. "Not yet run in MATLAB" until the MATLAB job is green.
4. **Who did what:** the block from §0.
5. **Hero figure, `hero.png`:** the top row of pipeline.png, without "quality 1.00".
   - Panels: one exposure of the real Basler texture | the 9-flash image | the log spectrum with the direction arrow and line spacing marked.
   - Caption: "Simulated from a real surface texture: estimate 31.0°, 5.83 px/flash; true 30.96°, 5.83".
6. **Try it, 2 lines:** `pip install -e "./python[examples]"` and `python -m fftvel simulate --vx 5 --vy -3 --plot estimate.png`.

**Five-minute reader**

7. **Why one image.**
   - The time base is the flash period (µs), not the frame rate (ms).
   - Example: 30 m/s with a 0.2 mm/px, 256 px field. At 100 fps consecutive frames are 300 mm apart and do not overlap. Flashes 133 µs apart give 20 px per flash.
   - Say plainly: "If you have overlapping separate frames, use phase correlation; it is more accurate and gives the sign." The live frame-averaging scripts are a test harness.
   - Small table: phase correlation (oracle), cepstrum and autocorrelation / Young's fringes (baselines). Link to related-work.md.
8. **How it works.**
   - One equation, $B_p = B\cdot D(\mathbf f\cdot\mathbf v)$, and four steps, each tagged [thesis] or [rework]: spectrum, direction (mean Radon projection, ripple profile WR), speed (Dirichlet fit), validity.
   - No "Chapter 3" in headings. Link to method.md.
9. **What I built on top of the thesis.** A provenance table that replaces the chapter structure.
   - Columns: Step | Thesis (Horn) | 2025 student project (up to `8f0aec4`) | 2026 rework (AI-assisted).
   - Rows: Image, Direction, Speed, Validity, Code (content as in reader-first's table).
10. **Results, ordered by how independent they are.** First a compact ladder table: rung | data | independent of the estimator's model? | result | next step. Then three short claims:
    - **(a) Synthetic test sequences** (`testsequences.png`, the only figure with an independent reference). Constant motion: 135.0 ± 0.07°, speed within 0.004 px. Accelerating sequence: 17/22 windows valid. Irregular steps: 4/22 valid.
    - **(b) The 2025 recordings break the model and are rejected** (`recordings.png`).
    - **(c) Simulation.**
      - Noise-free: ≤0.71° and ≤0.5 % from 0.4 to 64 px/flash.
      - With noise, per texture, camera texture first (verified from noise.csv, n = 60 each): 70 % / 50 % correct at noise = 0.5 / 1.0 × signal; random texture 92 % / 72 %.
      - 0 of 120 still images marked valid (≤2.5 % at 95 %).
      - 11 of 592 valid estimates wrong (1.9 %, 95 % CI 0.9–3.3 %), five of them at twice the speed.
      - Add one sentence: these images come from the same model the estimator fits, so they verify the implementation, not the model on a real camera.
    - Link to results.md.
11. **Where it fails.**
    - Silent failures first (wrong result marked valid):
      - regular uneven steps or a dropped flash: about 2× speed, can pass;
      - speed change during the burst: up to 11 % error, 17/22 windows pass;
      - a sharp static edge: false motion;
      - periodic texture, continuous blur, ambient light: untested;
      - wrong L: not detected, about 1 % effect.
    - Then the flagged ones: motion below the measurable range, heavy noise, no texture or no motion.
    - Direction is only known modulo 180°.
12. **Use it.**
    - The Python CLI and API.
    - MATLAB/Octave: `addpath('matlab'); est = fftvel.estimate_velocity(bp,'L',4); fftvel.plot_estimate(bp,est)`.
    - Units: v = s·GSD/T, with a link to hardware.md.
    - Live camera, 3 lines: "EXPERIMENTAL test harness, never run on camera hardware; the Python path was checked with the pylon emulator and video files".
    - Runtime: 0.7 s in Python and 1.4 s in Octave for 256 px; not real-time.
13. **Development.**
    - Tests: 2 commands.
    - What the tests prove: self-consistency and that the two languages agree, not that the model holds.
    - Python is the reference; to change the algorithm, change both. Link to CONTRIBUTING.md.
    - An 8-line layout and the `blob:none` clone tip.
14. **How this was built:** the full text from §0.
15. **Reference, citation, acknowledgements, licence.**
    - Full Horn reference; Chapter 3 is named only here.
    - CITATION.cff.
    - Acknowledgements only after asking the people named.
    - Code MIT; data terms in data/README.md.

## 3. Target repo layout

```
2D-FFT-Code/
├── README.md  LICENSE  CITATION.cff (new)  CONTRIBUTING.md (new: Python↔MATLAB map, change checklist)
├── .github/workflows/ci.yml   python 3.10/3.13 + pip cache; octave + demos (xvfb-run); matlab (new); experiments --quick smoke (new)
├── python/
│   ├── pyproject.toml  LICENSE  README.md (pointer)
│   ├── src/fftvel/  __init__ (smaller __all__, version via importlib.metadata), core, constants (new), simulate,
│   │                reference, plotting (+ hero figure), stream, __main__ (+ --gsd --period)
│   └── examples/    demo_simulation.py  demo_testsequences.py  live_camera.py (EXPERIMENTAL header)
├── matlab/
│   ├── +fftvel/     estimate_velocity simulate_pulsed_image plot_estimate phase_correlation make_texture
│   │   │            angle_difference (new) constants (new)
│   │   └── +internal/  crop_square dirichlet_template fit_speed hann_window log_spectrum mean_projection
│   │                   parabolic_offset parse_options ripple_profile null_quality (new, split out)
│   └── examples/    demo_simulation.m  demo_testsequences.m  live_basler.m (EXPERIMENTAL header)
├── tests/  python/  matlab/  parity/ (+ make_cases.py --check, simulate case, stage values)  constants.json (new)
├── experiments/                     (was benchmark/)
│   ├── run.py  accuracy.py  noise.py  test_sequences.py  recordings_2025.py
│   ├── calibrate_threshold.py  baselines.py  model_violations.py  timing.py   (new)
│   ├── history/  original_2025.py  before_after.py
│   └── results/  *.csv (+ threshold_calibration, baselines, violations, environment.json)
├── data/   README.md  testsequences.mat  surface_texture.png  (bench_2026/ after WP11)
├── tools/  make_data.py (new)  fetch_recordings.sh (new)
└── docs/   method.md  results.md  hardware.md  history.md  related-work.md   figures/ (+ hero.png)
```

| Action | From | To / note |
|---|---|---|
| Remove from tree (stays at tag) | `2D-FFT/` | linked as `/tree/student-project-2025/2D-FFT` |
| Remove from tree (stays at tag) | `Matlab_BaslerCamDriver/` | one line in history.md |
| Move (own commit; Octave suite and parity must stay green) | `matlab/{estimate_velocity,simulate_pulsed_image,plot_estimate,phase_correlation,make_texture}.m` | `matlab/+fftvel/` |
| Move | the 9 helper `.m` files | `matlab/+fftvel/+internal/` |
| Rename | `benchmark/` | `experiments/` |
| Rename and split | `run_benchmark.py` | `run.py` + one script per study; `--quick --out DIR` never writes `results/` |
| Rename | `benchmark/original_estimator.py` | `experiments/history/original_2025.py` (port of `aktuell/finaltest1.m` at the tag) |
| Merge the 3 copies | `tests/matlab/angle_error.m`, helper in `conftest.py`, benchmark copy | `fftvel.angle_difference` in both languages |
| Move text | derivation, constants, full Limitations | `docs/method.md` (+ `pipeline.png`) |
| Move text | "Background", before/after, 733 % | `docs/history.md` (+ `before_after.png`) |
| Move text | live camera details | `docs/hardware.md` |
| Retitle figures | `testsequences.png` | "Synthetic test sequences with known motion" |
| Retitle figures | `noise.png` | one panel per texture, camera texture first |
| Retitle figures | `quality.png` | full legend, sample sizes, "(256 px)" at the threshold |
| Delete | README sections "Background" and "Tests" as now; the label "Ground truth" | |
| Fix | `__version__ = "1.0.0"` | 0.1.0 |
| New | docs pages, `hero.png`, CITATION.cff, CONTRIBUTING.md, `tools/*`, new experiment scripts, `tests/constants.json`, MATLAB CI job | |

## 4. 2D-FFT/ and Matlab_BaslerCamDriver/

**2D-FFT/ leaves the working tree.** It stays byte-identical at tag `student-project-2025` (= `8f0aec4`) in 2D-FFT-Code itself, because the history is pushed unchanged.
- **docs/history.md replaces the chapter structure** with an index by function:
  - Image formation and simulation: Kapitel21_v3, Kapitel33_v1, bilddarstellung, simulateTexture, generateImageSequence, exponentialFilter.
  - Single-image analysis: Kapitel34_v1, main, untitled.
  - Pixel reduction: Kapitel35_v1.
  - GEBA: Kapitel34_v2, Kapitel36_v1, video_test.
  - Camera and live: `aktuell/finaltest1.m` (source of `original_2025.py`), kamera_live, Chapter3_4_live*, camera3_feed_sample, Treiber_install.
  - Video GUIs: radon_video_gui_*, test.m.
  - Unrelated: `aktuell/polyphasen_filterbank.m`.
- **Also recorded in history.md:**
  - The file names contradict their own headers: `…Kapitel33_v1` is headed "3.1.3", `…35_v1` "3.3", `…36_v1` "3.4" (checked).
  - 11 scripts exist twice (root vs. `aktuell/Grundlagen_Skripte`).
  - The review findings: an absolute prominence threshold, anisotropic resize to 256 px, unequal steps.
  - The 733 % before/after, presented as code history.
- **GEBA:** the scripts describe it as "Gütemaß … Selbstbeurteilung in der Einzelbild-Analyse (Kap. 3.4)". Take the exact expansion from the thesis before writing it out.
- **Recordings and data files:**
  - `tools/fetch_recordings.sh` extracts vid1–4 from the tag (`git archive student-project-2025 2D-FFT/data/…`) for `recordings_2025.py` and the real null test.
  - `tools/make_data.py` rebuilds `testsequences.mat` and `surface_texture.png` from the tag and checks their SHA-256.

**Matlab_BaslerCamDriver/ leaves the working tree** and stays at the tag.
- It is third-party code, and no script calls it; the 2025 code uses `videoinput('gentl')`.
- It is mentioned only in one line of history.md and in the "Who did what" exception. Nowhere else.

**onewmann/hs-esslingen:**
- Do not rewrite it and do not touch anything inside 2D-FFT/.
- Keep the `claude/pensive-gates-pzcxla` branch.
- Add a root README: "Unchanged since 8f0aec4; the rework is in 2D-FFT-Code". That commit is AI-era, so say so, or the student writes it themselves.
- Then archive the repo read-only.

## 5. Work packages, in priority order

**Phase 0: before anything else is pushed**
- **WP0: Data rights.** Effort S. Needs the student.
  - Ask the supervisor about:
    - `Testsequenzen1.mat`: who generated it?
    - vid1–4 and `surface_texture.png` (cut from vid2);
    - publishing the code at all, and naming the supervisor and HS Esslingen.
  - The data are already public in hs-esslingen. If the answer is no, both repos need the fallback (filter-repo out `2D-FFT/data`).
  - Answers: data rights, acknowledgements.

**Phase 1: restructure without any algorithm change.** CSVs, parity cases and numbers stay identical. The assistant can do this alone in this session; Octave 8.4 is available here.
- **WP1: AI disclosure.** Effort S. The assistant writes the text from §0 and the command `git tag -a student-project-2025 8f0aec4 -m "Student project files as uploaded; AI-assisted rework starts after this commit"`. The student writes "My part" and decides D3. Answers: authorship and the user's request.
- **WP2: Tree restructure.** Effort M.
  - `git rm -r 2D-FFT Matlab_BaslerCamDriver` in one commit whose message names the tag.
  - `benchmark/` → `experiments/`.
  - The `+fftvel/+internal` move as its own commit.
  - Update the test paths, demos, `live_basler.m` and CI.
  - Merge the three angle helpers; set version 0.1.0.
  - Answers: generic MATLAB names, chapter structure, third-party code.
- **WP3: README, docs split, figures.** Effort M.
  - `hero.png`; retitled and per-texture figures.
  - The README as in §2.
  - method.md with a constants table and GEBA; results.md; history.md; hardware.md with acquisition requirements, the v = s·GSD/T formula and the live-script status.
  - Answers: the 30-second test, paper-like README, quick start too late, simulation numbers unlabelled, hidden failure on real data, numbers per texture, before/after strawman, MATLAB wording, GEBA, use case (text only). The student reviews.
- **WP4: CI.** Effort S.
  - Python 3.10 and 3.13 with pip cache.
  - Octave job plus demos under `xvfb-run`.
  - MATLAB job: `matlab-actions/setup-matlab` + `cd tests/matlab; run_tests`.
  - The assistant writes it. It only runs once 2D-FFT-Code exists on GitHub.
  - Answers: MATLAB never run (once green), packaging.
- **WP5: Publish.** Effort S. The student creates the empty repo onewmann/2D-FFT-Code. The assistant can push only after `add_repo` with push access; this session has no access to 2D-FFT-Code yet.
  - Push only `HEAD:main` of the local branch (the student-authored chain including `680a45f`) and the tag. Not `--all`, so `backup-before-author-rewrite` stays local.
  - Set About, topics and the social preview; create release v0.1.0.

**Phase 2: new computation, no hardware.** The assistant can do this alone.
- **WP6: Reproducible threshold calibration.** Effort M.
  - Commit the scratch `calib/calib.py` as `experiments/calibrate_threshold.py`. It already rebuilds the 1008-image null set: 7 sizes × 4 L × 3 null kinds × 12.
  - A scratch run from before the last fixes gave a maximum k of 8.5; the docstring says 8.6. Report whatever the committed script gives; do not tune it towards 8.6.
  - Add held-out seeds and textures.
  - Add a real null test on vid1 and single frames of vid2–4, with a Clopper–Pearson bound.
  - Add an ROC of fit quality vs. GEBA and a score check at v/2 vs. 2v.
  - Answers: magic number, a valid flag that lets 2× through, GEBA unused.
- **WP7: Baselines.** Effort M.
  - `experiments/baselines.py`, kept out of the fftvel API.
  - Cepstrum peak; autocorrelation peak; the thesis's peak-spacing rule done correctly (robust prominence, sub-bin refinement); phase-correlation oracle.
  - Same images and metrics as the noise study; `baselines.png`.
  - Answers: no baselines, strawman comparison.
- **WP8: Independent forward model and failure modes.** Effort M–L.
  - `model_violations.py`: `scipy.ndimage.shift`, finite pulse width, step and energy jitter, Poisson plus read noise, 8-bit quantisation, ambient light, defocus.
  - Tests for periodic texture and a static edge; `failures.png`; confidence intervals per texture.
  - Answers: circular validation, static and periodic structures, thin statistics.
- **WP9: Keeping the two implementations in sync.** Effort M.
  - `constants.py` / `constants.m` plus `tests/constants.json`, checked on both sides.
  - Stage-level parity values; a parity case for `simulate_pulsed_image` with a fixed texture; `make_cases.py --check` in CI.
  - CONTRIBUTING.md.
  - API: L keyword-only, frozen `Estimate`, smaller `__all__`.
  - Answers: hand-synced implementations, API width.
- **WP10: Units and runtime.** Effort M.
  - `--gsd` / `--period` options.
  - An error budget in hardware.md: GSD, flash jitter, and the height sensitivity of an entocentric lens.
  - A precomputed sparse Radon operator and coarse-to-fine angle search in both languages; `timing.py` and a timing table.
  - This changes the algorithm, so regenerate the parity cases, rerun the CSVs and update the README numbers in the same commit.
  - Answers: no physical units, no real-time.

**Phase 3: needs the student. Can run in parallel with phase 2.**
- **WP11: Real multi-flash experiment.** Effort L. Needs hardware and the lab.
  - **Set-up:**
    - Basler camera.
    - Power LED switched by a MOSFET.
    - A microcontroller as the timing master: it triggers the camera and fires 2L+1 flashes of duration τ at period T inside one exposure.
    - A constant-speed linear axis or turntable with sandpaper or asphalt.
    - GSD from a checkerboard image.
  - **Shots:**
    - L = 4: 6 speeds × 10 shots;
    - L = 10: 2 low speeds × 10;
    - about 100 static shots with the flashes on;
    - 10 shots with one flash dropped;
    - 10 shots with ambient light.
  - **Two references per shot:** v_set·T/GSD, and phase correlation of two single-flash frames.
  - **Assistant's part, in advance:** firmware, acquisition script, a `meta.csv` schema, `bench_2026.py` and the figure.
  - **Afterwards:** the hero image and the status line switch to real data.
  - **Fallback:** averaged hardware-triggered short exposures, labelled "emulated".
  - Answers: no real data, use case, units, real false-alarm rate.
- **WP12: The student's own understanding.** Effort M. Needs the student.
  - Derive P(j), the Hann autocorrelation and the m−5 / Fisher-z threshold on paper.
  - Rewrite the method.md derivations in the student's own words.
  - Read and check every related-work citation; write "My part" and the acknowledgements.
  - Answers: authorship and the competence gap, related work.
- **WP13: One MATLAB run, only if the CI job cannot run.** Effort S. Needs the student and the campus licence.

## 6. Decisions for the student, with recommendations

- **D1: Data rights.** Ask the supervisor now. Recommendation: push nothing new until the supervisor answers.
- **D2: History.** Recommendation: push it unchanged (C1). filter-repo only if D1 is a no.
- **D3: hs-esslingen and the Claude-authored branch.** Recommendation: keep both public, name the branch in "How this was built", and archive the repo read-only. Making it private or deleting the branch after publication would look like hiding.
- **D4: Trailer from now on.** Decided by the student: no Claude trailer or author; the README section "How this was built" is the disclosure.
- **D5: Real experiment (WP11).** Recommendation: yes, strongly. Without it, README and CV must keep saying "validated in simulation only".
- **D6: Drop Python 3.9.** Recommendation: yes, test 3.10 and 3.13; 3.9 is end-of-life.
- **D7: Name the supervisor and HS Esslingen.** Recommendation: only with explicit consent.

## 7. What not to do

- Do not describe the AI use by a date. Always give the commit, its link and the compare link.
- Do not rewrite or squash the history, do not use filter-repo (unless D1 forces it), do not add trailers to old commits, and do not push with `--all`.
- Do not delete the claude/* branch silently, and do not make hs-esslingen private after the README mentions it.
- Do not copy, move or rename the 2025 files, do not add a submodule, and do not edit inside 2D-FFT/ in the archive.
- Do not put the Basler driver in the new repo, and name it as third-party code wherever the 2025 upload is described.
- In MATLAB code, do not use `+fftvel/private/`, `import fftvel.*`, or `exist()` on package functions. Each fails or misbehaves in Octave 8.4.
- Do not change the algorithm in phase 1.
- Do not lead with `pipeline.png` "quality 1.00" or with simulation numbers that lack "in simulation".
- Do not use "Ground truth" for the test sequences, and do not keep the 733 % as a headline.
- Do not use "Chapter", "Kapitel" or "Background" as section titles or file names.
- Do not tune the threshold on data you then report the false-alarm rate on, and do not tune the calibration script to reproduce 8.6.
- Do not put the baselines into the fftvel public API.
- Do not present frame averaging as the use case.
- Do not cite related work nobody has read, and do not guess the expansion of GEBA.
- Do not let `--quick` or CI write into `experiments/results` or `docs/figures`, and do not compare PNG pixels in CI.
- Do not write "MATLAB-tested", "real-time" or "measurement system" in the README, CV or LinkedIn before the matching CI job or measurement exists. CV line: "Single-image surface velocity estimation (2-D FFT, Radon): MATLAB implementation of a doctoral-thesis method (student project 2025); reworked with AI assistance into a tested Python/MATLAB package (2026). Validated in simulation; real-data validation open."
