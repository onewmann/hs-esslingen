"""experiments/run.py and the study modules: argument checks and output folders (no simulation is run)."""

import re
import shutil
import subprocess
import sys

import pytest

STUDY_MODULES = ["run", "common", "accuracy", "noise", "pipeline_figure", "recordings_2025", "testsequences",
                 "calibrate_threshold", "baselines", "model_violations", "timing",
                 "history", "history.before_after", "history.original_2025"]
ALL_VIDEOS = "vid1.mat, vid2.mat, vid3.mat, vid4.mat not found"


def experiment(repo, *args, script="run.py"):
    """Run an experiments script with the given arguments; return the finished process."""
    return subprocess.run([sys.executable, str(repo / "experiments" / script), *map(str, args)],
                          capture_output=True, text=True, timeout=300)


@pytest.fixture
def experiments(repo, monkeypatch):
    """The modules run and common, imported from experiments/ and removed again afterwards."""
    monkeypatch.setattr(sys, "path", [str(repo / "experiments"), *sys.path])
    for name in STUDY_MODULES:
        monkeypatch.delitem(sys.modules, name, raising=False)
    import common
    import run
    yield run, common
    for name in STUDY_MODULES:
        sys.modules.pop(name, None)


@pytest.mark.parametrize("args, message", [
    (["--only", "accuracy", "--recordings", "{tmp}"], "--only leaves out recordings"),
    (["--recordings", "{tmp}"], "vid2.mat, vid3.mat, vid4.mat not found"),
    (["--only", "recordings", "--recordings", "{tmp}"], "vid2.mat, vid3.mat, vid4.mat not found"),
    (["--quick", "--out", "{committed}"], "--quick must not write to experiments/results"),
    (["--quick", "--out", "{docs}"], "--quick must not write to docs/figures"),
    (["--quick", "--all", "--out", "{committed}"], "--quick must not write to experiments/results"),
    (["--only", "baselines", "timing", "--recordings", "{tmp}"], "--only leaves out recordings and calibration"),
    (["--only", "calibration", "--recordings", "{tmp}"], ALL_VIDEOS),
    (["--all", "--recordings", "{tmp}"], ALL_VIDEOS),
    (["--all", "--only", "noise"], "argument --only: not allowed with argument --all"),
    (["--only", "real_measurements"], "invalid choice: 'real_measurements'"),
    (["--only", "violations", "--jobs", "0"], "--jobs must be at least 1"),
])
def test_argument_errors(repo, tmp_path, args, message):
    args = [a.format(tmp=tmp_path, committed=repo / "experiments" / "results", docs=repo / "docs") for a in args]
    r = experiment(repo, *args)
    assert r.returncode == 2 and message in r.stderr, r.stderr


@pytest.mark.parametrize("script, args, message", [
    ("recordings_2025.py", ["--recordings", "{tmp}"], "vid2.mat, vid3.mat, vid4.mat not found"),
    ("noise.py", ["--only", "accuracy"], "unrecognized arguments: --only"),
    ("noise.py", ["--recordings", "{tmp}"], "unrecognized arguments: --recordings"),
    ("history/before_after.py", ["--quick", "--out", "{committed}"], "--quick must not write"),
])
def test_study_module_argument_errors(repo, tmp_path, script, args, message):
    args = [a.format(tmp=tmp_path, committed=repo / "experiments" / "results") for a in args]
    r = experiment(repo, *args, script=script)
    assert r.returncode == 2 and message in r.stderr, r.stderr


def test_replot_uses_only_the_output_folder(repo, tmp_path):
    shutil.copy(repo / "experiments" / "results" / "testsequences.csv", tmp_path)
    figures = sorted((repo / "docs" / "figures").glob("*.png"))
    before = [p.stat().st_mtime_ns for p in figures]
    for script, args in (("run.py", ["--only", "testsequences"]), ("testsequences.py", [])):
        r = experiment(repo, "--replot", *args, "--out", tmp_path, script=script)
        assert r.returncode == 0, r.stderr
        assert "windows valid" in r.stdout
        assert [p.name for p in (tmp_path / "figures").iterdir()] == ["testsequences.png"]
    assert [p.stat().st_mtime_ns for p in figures] == before       # README figures untouched


def test_replot_needs_the_csv_file(repo, tmp_path):
    for script, args in (("run.py", ["--only", "recordings"]), ("recordings_2025.py", [])):
        r = experiment(repo, "--replot", *args, "--out", tmp_path, script=script)
        assert r.returncode == 1 and "recordings.csv not found" in r.stderr, r.stderr


@pytest.mark.parametrize("part", ["calibration", "baselines", "violations", "timing"])
def test_replot_of_a_part_on_request_needs_its_csv_files(tmp_path, experiments, part):
    run, common = experiments
    with pytest.raises(SystemExit, match=f"--replot {part}: .*{re.escape(run.STUDIES[part].CSV[0])}.* not found"):
        run.main(["--replot", "--only", part, "--out", str(tmp_path)])
    assert common.RESULTS == tmp_path.resolve()


def test_replot_all_skips_parts_without_csv_files(repo, tmp_path):
    for name in ("testsequences.csv", "timing.csv"):
        shutil.copy(repo / "experiments" / "results" / name, tmp_path)
    r = experiment(repo, "--replot", "--all", "--out", tmp_path)
    assert r.returncode == 0, r.stderr
    for part in ("accuracy", "noise", "before_after", "recordings", "calibration", "baselines", "violations"):
        assert f"\n{part}: skipped, " in r.stdout, part
    assert "windows valid" in r.stdout and "Run time of estimate_velocity" in r.stdout
    # the pipeline figures have no CSV file and are computed again; timing has no figure
    assert sorted(p.name for p in (tmp_path / "figures").iterdir()) == ["hero.png", "pipeline.png",
                                                                         "testsequences.png"]


def test_default_folders(repo, experiments):
    run, common = experiments
    assert common.COMMITTED == repo / "experiments" / "results"
    assert common.QUICK == repo / "experiments" / "results_quick"
    assert common.RECORDINGS == repo / "data" / "recordings"
    assert run.PARTS == ["accuracy", "noise", "before_after", "pipeline", "testsequences", "recordings"]
    assert run.EXTRA == ["calibration", "baselines", "violations", "timing"]
    assert run.ALL_PARTS == run.PARTS + run.EXTRA and list(run.STUDIES) == run.ALL_PARTS
    for part in run.ALL_PARTS:
        study = run.STUDIES[part]
        assert isinstance(study.CSV, list) and callable(study.main) and callable(study.replot), part
    assert set(run.SIMULATIONS) <= set(run.ALL_PARTS) and set(run.SUMMARISED) <= set(run.SIMULATIONS)
    assert run.RECORDING_FILES["calibration"] == ["vid1.mat", "vid2.mat", "vid3.mat", "vid4.mat"]
    assert run.RECORDING_FILES["recordings"] == ["vid2.mat", "vid3.mat", "vid4.mat"]


def test_quick_defaults_to_results_quick(tmp_path, monkeypatch, experiments):
    run, common = experiments
    quick = (tmp_path / "quick").resolve()
    monkeypatch.setattr(common, "QUICK", quick)
    with pytest.raises(SystemExit, match="recordings.csv not found"):
        run.main(["--quick", "--replot", "--only", "recordings"])
    assert common.RESULTS == quick and common.FIGURES == quick / "figures"


def test_ci_runs_known_parts(repo, experiments):
    # the experiments smoke step of .github/workflows/ci.yml names parts of run.py and a synthetic measurement session
    run, _ = experiments
    text = (repo / ".github" / "workflows" / "ci.yml").read_text()
    commands = re.findall(r"python experiments/run\.py ([^\n]*)", text)
    assert commands
    for args in commands:
        words = args.split()
        assert "--quick" in words and "--out" in words, args
        if "--only" in words:
            named = words[words.index("--only") + 1:]
            named = named[:next((i for i, w in enumerate(named) if w.startswith("--")), len(named))]
            assert named and set(named) <= set(run.ALL_PARTS), args
    metas = re.findall(r"python experiments/real_measurements\.py (\S+)", text)
    assert metas and all((repo / m).is_file() for m in metas), metas


def test_no_script_looks_like_a_test_module(repo):
    # a bare pytest run from the repository root would import such a script as a test module
    scripts = sorted((repo / "experiments").rglob("*.py"))
    assert scripts
    assert [p.name for p in scripts if p.name.startswith("test_") or p.stem.endswith("_test")] == []
