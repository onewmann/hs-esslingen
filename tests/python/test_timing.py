"""experiments/timing.py: a short run on small images, the table and the argument checks."""

import csv
import json
import subprocess
import sys

import pytest


def timing(repo, *args):
    return subprocess.run([sys.executable, str(repo / "experiments" / "timing.py"), *map(str, args)],
                          capture_output=True, text=True, timeout=300)


def test_short_run(repo, tmp_path):
    r = timing(repo, "--sizes", 32, 48, "--first-runs", 2, "--repeats", 2, "--no-octave", "--out", tmp_path)
    assert r.returncode == 0, r.stderr
    with open(tmp_path / "timing.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    assert [(row["implementation"], int(row["n"]), int(row["npad"])) for row in rows] == \
        [("python", 32, 64), ("python", 48, 96)]
    for row in rows:
        assert int(row["first_runs"]) == 2 and int(row["repeated_runs"]) == 2 and int(row["L"]) == 4
        assert 0 < float(row["repeated_min_s"]) <= float(row["repeated_median_s"]) < 60
        assert 0 < float(row["first_min_s"]) <= float(row["first_median_s"]) < 60
        assert float(row["import_median_s"]) > 0
    env = json.loads((tmp_path / "environment.json").read_text())
    for key in ("platform", "cpu_count", "python", "numpy", "fftvel", "date_utc"):
        assert env[key], key
    assert "median in s" in r.stdout and "python     32    64" in r.stdout
    r = timing(repo, "--replot", "--out", tmp_path)
    assert r.returncode == 0 and "python     48    96" in r.stdout, r.stderr


@pytest.mark.parametrize("args, message", [
    (["--quick", "--out", "{committed}"], "--quick must not write to experiments/results"),
    (["--sizes", "8"], "--sizes must be at least 16"),
    (["--repeats", "0"], "--repeats must be at least 1"),
    (["--first-runs", "0"], "--first-runs must be at least 1"),
])
def test_argument_errors(repo, args, message):
    args = [a.format(committed=repo / "experiments" / "results") for a in args]
    r = timing(repo, *args)
    assert r.returncode == 2 and message in r.stderr, r.stderr


def test_replot_needs_the_csv_file(repo, tmp_path):
    r = timing(repo, "--replot", "--out", tmp_path)
    assert r.returncode == 1 and "timing.csv not found" in r.stderr, r.stderr
