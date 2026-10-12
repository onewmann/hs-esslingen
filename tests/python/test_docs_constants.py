"""The Constants table of docs/method.md is a copy of tests/constants.json.

test_constants.py checks the package against the JSON file; this file checks
the table in the documentation against the same file, so that a changed
constant cannot leave the documented value behind.
"""

import html
import json
import re

import pytest

NAME = re.compile(r"`([A-Z][A-Z0-9_]*)`")


def cells(line):
    """Cells of a Markdown table row; an escaped pipe (\\|) stays inside its cell."""
    line = line.strip()
    if not (line.startswith("|") and line.endswith("|") and len(line) > 1):
        return None
    return re.split(r"(?<!\\)\|", line[1:-1])


def parse_number(text):
    text = text.strip()
    return int(text) if re.fullmatch(r"-?\d+", text) else float(text)


def unescape(cell):
    """Markdown cell text to plain text: HTML entities and backslash escapes."""
    return re.sub(r"\\(.)", r"\1", html.unescape(cell.strip()))


@pytest.fixture(scope="module")
def spec(repo):
    with open(repo / "tests" / "constants.json", encoding="ascii") as f:
        return json.load(f)["constants"]


@pytest.fixture(scope="module")
def table(repo):
    """Rows of the table under '## Constants' as {name: (value cell, meaning cell)}."""
    text = (repo / "docs" / "method.md").read_text(encoding="utf-8")
    match = re.search(r"^## Constants\n(.*?)(?=^## )", text, flags=re.M | re.S)
    assert match, "docs/method.md has no '## Constants' section"
    rows = {}
    for line in match.group(1).splitlines():
        row = cells(line)
        if row is None or len(row) != 3:
            continue
        m = NAME.fullmatch(row[0].strip())
        if m:
            name = m.group(1)
            assert name not in rows, f"{name} appears twice in the table"
            rows[name] = (row[1], row[2])
    assert rows, "no table rows found under '## Constants'"
    return rows


def test_same_names_in_the_same_order(spec, table):
    assert list(table) == list(spec)


def test_same_values(spec, table):
    for name, entry in spec.items():
        cell = table[name][0]
        expected = entry["value"] if isinstance(entry["value"], list) else [entry["value"]]
        got = [parse_number(t) for t in cell.split(",")]
        assert got == expected, f"{name}: table {cell.strip()!r}, constants.json {entry['value']!r}"
        assert [type(v) for v in got] == [type(v) for v in expected], f"{name}: int against float"


def test_same_meanings(spec, table):
    for name, entry in spec.items():
        assert unescape(table[name][1]) == entry["meaning"], name
