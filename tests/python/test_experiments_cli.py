"""experiments/run.py and the study modules: argument checks and output folders (no simulation is run)."""

import shutil
import subprocess
import sys

import pytest

STUDY_MODULES = ["run", "common", "accuracy", "noise", "pipeline_figure", "recordings_2025", "test_sequences",
                 "history", "history.before_after", "history.original_2025"]


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
    for script, args in (("run.py", ["--only", "testsequences"]), ("test_sequences.py", [])):
        r = experiment(repo, "--replot", *args, "--out", tmp_path, script=script)
        assert r.returncode == 0, r.stderr
        assert "windows valid" in r.stdout
        assert [p.name for p in (tmp_path / "figures").iterdir()] == ["testsequences.png"]
    assert [p.stat().st_mtime_ns for p in figures] == before       # README figures untouched


def test_replot_needs_the_csv_file(repo, tmp_path):
    for script, args in (("run.py", ["--only", "recordings"]), ("recordings_2025.py", [])):
        r = experiment(repo, "--replot", *args, "--out", tmp_path, script=script)
        assert r.returncode == 1 and "recordings.csv not found" in r.stderr, r.stderr


def test_default_folders(repo, experiments):
    run, common = experiments
    assert common.COMMITTED == repo / "experiments" / "results"
    assert common.QUICK == repo / "experiments" / "results_quick"
    assert common.RECORDINGS == repo / "data" / "recordings"
    assert run.PARTS == ["accuracy", "noise", "before_after", "pipeline", "testsequences", "recordings"]
    for part in run.PARTS:
        assert isinstance(run.STUDIES[part].CSV, list) and callable(run.STUDIES[part].main), part


def test_quick_defaults_to_results_quick(tmp_path, monkeypatch, experiments):
    run, common = experiments
    quick = (tmp_path / "quick").resolve()
    monkeypatch.setattr(common, "QUICK", quick)
    with pytest.raises(SystemExit, match="recordings.csv not found"):
        run.main(["--quick", "--replot", "--only", "recordings"])
    assert common.RESULTS == quick and common.FIGURES == quick / "figures"


def test_study_modules_have_no_test_functions(repo):
    # experiments/test_sequences.py matches pytest's file pattern; nothing in it may look like a test
    text = (repo / "experiments" / "test_sequences.py").read_text()
    assert "\ndef test" not in text and "\nclass Test" not in text
