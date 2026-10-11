"""benchmark/run_benchmark.py: argument checks and output folders (no simulation is run)."""

import importlib.util
import shutil
import subprocess
import sys

import pytest


def bench(repo, *args):
    """Run the benchmark script with the given arguments; return the finished process."""
    return subprocess.run([sys.executable, str(repo / "benchmark" / "run_benchmark.py"), *map(str, args)],
                          capture_output=True, text=True, timeout=300)


@pytest.mark.parametrize("args, message", [
    (["--only", "recordings"], "--only recordings needs --recordings"),
    (["--only", "accuracy", "--recordings", "{tmp}"], "--only leaves out recordings"),
    (["--recordings", "{tmp}"], "vid2.mat, vid3.mat, vid4.mat not found"),
    (["--quick", "--out", "{committed}"], "--quick must not write to benchmark/results"),
])
def test_argument_errors(repo, tmp_path, args, message):
    args = [a.format(tmp=tmp_path, committed=repo / "benchmark" / "results") for a in args]
    r = bench(repo, *args)
    assert r.returncode == 2 and message in r.stderr, r.stderr


def test_replot_uses_only_the_output_folder(repo, tmp_path):
    shutil.copy(repo / "benchmark" / "results" / "testsequences.csv", tmp_path)
    figures = sorted((repo / "docs" / "figures").glob("*.png"))
    before = [p.stat().st_mtime_ns for p in figures]
    r = bench(repo, "--replot", "--only", "testsequences", "--out", tmp_path)
    assert r.returncode == 0, r.stderr
    assert "windows valid" in r.stdout
    assert [p.name for p in (tmp_path / "figures").iterdir()] == ["testsequences.png"]
    assert [p.stat().st_mtime_ns for p in figures] == before       # README figures untouched


def test_replot_needs_the_csv_file(repo, tmp_path):
    r = bench(repo, "--replot", "--only", "recordings", "--out", tmp_path)
    assert r.returncode == 1 and "recordings.csv not found" in r.stderr, r.stderr


def test_quick_defaults_to_results_quick(repo, tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "path", list(sys.path))                  # the script extends sys.path
    spec = importlib.util.spec_from_file_location("run_benchmark", repo / "benchmark" / "run_benchmark.py")
    rb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rb)
    assert rb.QUICK == repo / "benchmark" / "results_quick"
    quick = (tmp_path / "quick").resolve()
    monkeypatch.setattr(rb, "QUICK", quick)
    monkeypatch.setattr(sys, "argv", ["run_benchmark.py", "--quick", "--replot", "--only", "recordings"])
    with pytest.raises(SystemExit, match="recordings.csv not found"):
        rb.main()
    assert rb.RESULTS == quick and rb.FIGURES == quick / "figures"
