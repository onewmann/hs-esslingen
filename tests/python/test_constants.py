"""The constants in python/src/fftvel/core.py equal those of tests/constants.json.

tests/matlab/test_constants.m checks matlab/+fftvel/+internal/constants.m
against the same file, so both implementations use the same numbers.
"""

import inspect
import json

import pytest

from fftvel import core, estimate_velocity, simulate_pulsed_image
from fftvel.simulate import make_texture


@pytest.fixture(scope="module")
def spec(repo):
    with open(repo / "tests" / "constants.json", encoding="ascii") as f:
        return json.load(f)["constants"]


def module_constants():
    """Upper-case numbers and tuples of numbers defined in fftvel.core."""
    return {k: v for k, v in vars(core).items()
            if k.isupper() and isinstance(v, (int, float, tuple)) and not isinstance(v, bool)}


def test_same_names(spec):
    # a new constant must be added to the JSON file and to constants.m as well
    assert sorted(module_constants()) == sorted(spec)


def test_same_values(spec):
    for name, entry in spec.items():
        value = getattr(core, name)
        if isinstance(value, tuple):
            assert list(value) == entry["value"], name
        else:
            assert value == entry["value"] and type(value) is type(entry["value"]), name
        assert entry["meaning"] and entry["used_in"], name


def test_constants_fit_together():
    # the null-quality offset counts the points taken by the polynomial and by the fit
    assert core.NULL_OFFSET == core.TREND_DEGREE + 2
    assert (180 / core.ANGLE_STEP).is_integer()
    assert abs(sum(core.GREY_WEIGHTS) - 1) < 1e-3


def test_defaults_come_from_the_constants():
    p = inspect.signature(estimate_velocity).parameters
    assert p["L"].default == core.DEFAULT_L
    assert p["pad_factor"].default == core.DEFAULT_PAD_FACTOR
    assert p["min_quality"].default == core.DEFAULT_MIN_QUALITY
    assert inspect.signature(core.log_spectrum).parameters["pad_factor"].default == core.DEFAULT_PAD_FACTOR
    assert inspect.signature(core.ripple_profile).parameters["rfrac"].default == core.RFRAC
    assert inspect.signature(simulate_pulsed_image).parameters["corr_len"].default == core.DEFAULT_CORR_LEN
    assert inspect.signature(make_texture).parameters["corr_len"].default == core.DEFAULT_CORR_LEN


def test_options_after_the_image_are_keyword_only():
    for fn, first in ((estimate_velocity, ["bp"]), (core.pulsed_from_stack, ["frames"]),
                      (simulate_pulsed_image, ["n", "v", "L"])):
        params = inspect.signature(fn).parameters.values()
        positional = [p.name for p in params if p.kind is p.POSITIONAL_OR_KEYWORD]
        assert positional == first, fn.__name__
    with pytest.raises(TypeError):
        estimate_velocity(simulate_pulsed_image(64, (2, 1), 4, rng=0), 4)
