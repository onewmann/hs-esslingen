from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="session")
def repo() -> Path:
    return REPO


def angle_error(a: float, b: float) -> float:
    """Difference of two directions modulo 180 degrees, in [-90, 90)."""
    return (a - b + 90) % 180 - 90
