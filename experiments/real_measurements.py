"""Real multi-flash images against the set motion: estimate every image listed
in a meta.csv, convert to mm/s and summarise the errors.

    python experiments/real_measurements.py data/real/RUN/meta.csv [--figure PNG] [--out FILE]
    python experiments/real_measurements.py data/real/example/meta.csv --out results.csv

meta.csv has one row per image: the file, the number of flashes, the flash
period, the duty, the gsd, the set speed and direction, whether the surface
stood still, and optionally two single-flash frames for a phase-correlation
cross-check (format and advice: data/real/README.md). Writes results.csv
next to meta.csv (or to --out) with the estimate of every image, and prints
the share of moving images marked valid, the median and 95th percentile of
the speed and angle errors of the valid ones, the static images marked valid
with a Clopper-Pearson upper bound, and the worst images. --figure adds a
plot of estimated against set speed, the speed errors and the angle errors.

All rows are checked first. A row with a problem (a file that does not
exist, a value out of range) is reported with its line, skipped and listed
in results.csv with the reason; the other rows are still analysed. The exit
status is 0 when every row was analysed, 1 when some were not, and 2 when
meta.csv or the options cannot be used at all.

A file that holds a stack (a .mat from getdata, a .npy or a multi-page TIFF)
gives one image per frame, each estimated on its own; ``file.mat#3`` selects
frame 3 (counted from 1). Frames are never averaged here: every frame must be
a multi-flash image. data/real/example holds a small synthetic example.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import math
import re
import sys
import time
from functools import lru_cache
from multiprocessing import Pool
from pathlib import Path
from types import SimpleNamespace

import numpy as np

import common
from common import AQUA, BLUE, INK2, ORANGE, SURFACE
from fftvel import angle_difference, estimate_velocity, phase_correlation, to_physical
from fftvel.__main__ import _frames_axis, _load
from fftvel.core import DEFAULT_PAD_FACTOR, RFRAC, VMIN_FACTOR

COLUMNS = {
    "file": "image file, relative to meta.csv (.png, .tif, .npy, .mat); file.mat#k selects frame k",
    "exposures": "number of flashes K in the image (at least 2)",
    "flash_period_s": "time from one flash to the next in s (empty: speed in px per flash only)",
    "duty": "flash pulse length / flash period, in [0, 1); 0 or empty if unknown",
    "gsd_mm_per_px": "ground sample distance in mm per px (empty: speed in px per flash only)",
    "set_speed_mm_s": "speed set on the axis or turntable in mm/s (empty if unknown)",
    "set_angle_deg": "set direction, counter-clockwise from +x with y up (empty if unknown)",
    "static": "1 for an image that must come out not valid (stopped surface), else 0 (empty: 0)",
    "roi": "side of the analysed central square in px (empty: --roi, default 512; 0: largest square)",
    "note": "free text",
    "exposure_s": "camera exposure time in s, checked against (K - 1 + duty) flash periods (optional)",
    "ref_frame_a": "single-flash frame for the phase-correlation cross-check (optional)",
    "ref_frame_b": "second single-flash frame, taken ref_dt_s after ref_frame_a (optional)",
    "ref_dt_s": "time from ref_frame_a to ref_frame_b in s (empty: compare directions only)",
}
REQUIRED = ("file", "exposures")
META = ("exposures", "flash_period_s", "duty", "gsd_mm_per_px", "set_speed_mm_s", "set_angle_deg", "static",
        "exposure_s", "ref_frame_a", "ref_frame_b", "ref_dt_s", "note")    # meta.csv columns copied to results
# columns of results.csv in this order: the image, the outcome, the comparison with the set motion, the
# estimate, the settings from meta.csv, the cross-check, then warnings, errors and the note
FIELDS = ("file", "line", "valid", "reason", "speed_mm_s", "set_speed_mm_s", "speed_err", "angle_deg",
          "set_angle_deg", "angle_err_deg", "speed_px", "expected_px", "vx_px", "vy_px", "vx_mm_s", "vy_mm_s",
          "quality", "threshold", "quality_perp", "geba", "in_range", "isotropic", "roi", "static", "exposures",
          "flash_period_s", "duty", "gsd_mm_per_px", "exposure_s", "ref_frame_a", "ref_frame_b", "ref_dt_s",
          "ref_dx_px", "ref_dy_px", "ref_peak", "ref_angle_deg", "ref_angle_diff_deg", "ref_speed_px",
          "ref_speed_diff", "ref_speed_mm_s", "warnings", "error", "note")
DEFAULT_ROI = 512
TOL_SPEED = 0.05        # an estimate counts as correct within 5 % of the set speed ...
TOL_ANGLE = 2.0         # ... and 2 deg of the set direction, as in the simulated studies
CONFIDENCE = 0.95       # of the Clopper-Pearson upper bounds
SATURATED = 0.005       # warn when more than this share of the pixels sits at the maximum
FRAME = re.compile(r"^(.*)#(\d+)$")


class MetaError(ValueError):
    """meta.csv cannot be used at all (no header, a required column missing)."""


# meta.csv

def _decode(path: Path) -> str:
    """Text of a CSV file: UTF-8 (with or without BOM), else Windows-1252 as Excel writes it."""
    raw = path.read_bytes()
    for enc in ("utf-8-sig", "cp1252"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            pass
    return raw.decode("latin-1")


def _split_frame(text: str, folder: Path):
    """``(path, frame)`` of a file entry; ``a.mat#3`` gives frame 3 (from 1), no ``#k`` gives None."""
    t = text.strip().replace("\\", "/")          # a meta.csv written on Windows also works elsewhere
    path = folder / t
    m = FRAME.match(t)
    if m and not path.is_file():
        return folder / m.group(1), int(m.group(2))
    return path, None


def read_meta(path: Path, default_roi: int = DEFAULT_ROI):
    """Checked rows of meta.csv, the rows commented out, and the problems found.

    Comma-, semicolon- or tab-separated; with semicolons or tabs (German
    Excel) decimal commas are accepted. Empty rows are skipped, and so are
    rows whose file starts with ``#`` (they are counted as commented out).
    Returns ``(rows, commented, problems, unknown_columns)``; each row is a
    dict with the values of ``COLUMNS`` (None where empty), its line and
    ``problems``, the list of what is wrong with it (empty for a usable
    row). ``problems`` holds those of all rows, each with its line. Raises
    MetaError if the file is empty or lacks a required column.
    """
    text = _decode(path)
    header = next((ln for ln in text.splitlines() if ln.strip()), "")
    delim = max((",", ";", "\t"), key=header.count)
    comma = delim != ","
    reader = csv.reader(io.StringIO(text), delimiter=delim)
    names = None
    rows, commented, problems = [], [], []
    folder = path.parent
    for cells in reader:
        if not any(c.strip() for c in cells):
            continue
        if names is None:
            names = [c.strip().lower() for c in cells]
            missing = [c for c in REQUIRED if c not in names]
            if missing:
                raise MetaError(f"{path.name} has no column {', '.join(missing)} (header: {', '.join(names)}; "
                                f"separated by '{delim}')")
            continue
        line = reader.line_num
        raw = {k: (cells[i].strip() if i < len(cells) else "") for i, k in enumerate(names)}
        if raw["file"].startswith("#"):
            commented.append(f"line {line}: {raw['file']}")
            continue
        row, bad = _check_row(raw, folder, comma, default_roi)
        row["line"] = line
        if len(cells) > len(names):             # Excel writes as many cells in every row as in the header
            bad.insert(0, f"{len(cells)} cells but the header has {len(names)} columns; a decimal comma in "
                          f"a comma-separated file, or a comma in a note without quotes, shifts the values "
                          f"into the wrong columns (a comma at the end of the line also counts)")
        row["problems"] = bad
        rows.append(row)
        problems += [f"line {line}: {b}" for b in bad]
    if names is None:
        raise MetaError(f"{path.name} is empty")
    unknown = [c for c in names if c and c not in COLUMNS]
    return rows, commented, problems, unknown


def _check_row(raw: dict, folder: Path, comma: bool, default_roi: int):
    """One row of meta.csv as checked values, and the problems found in it."""
    bad = []

    def number(key, lowest=None, above=None, below=None, integer=False, default=None, choices=None):
        text = raw.get(key, "")
        if text == "":
            return default
        try:
            x = float(text.replace(",", ".") if comma else text)
        except ValueError:
            x = math.nan
        ok = math.isfinite(x) and (not integer or x == int(x))
        ok = ok and (lowest is None or x >= lowest) and (above is None or x > above)
        ok = ok and (below is None or x < below) and (choices is None or x in choices)
        if not ok:
            bad.append(f"{key} '{text}' must be {_rule(lowest, above, below, integer, choices)}")
            return default
        return int(x) if integer else x

    row = dict(line=None, file=raw["file"], note=raw.get("note", ""))
    row["exposures"] = number("exposures", lowest=2, integer=True)
    if raw["exposures"] == "":
        bad.append("exposures is empty; give the number of flashes in the image")
    row["flash_period_s"] = number("flash_period_s", above=0)
    row["duty"] = number("duty", lowest=0, below=1, default=0.0)
    row["gsd_mm_per_px"] = number("gsd_mm_per_px", above=0)
    row["set_speed_mm_s"] = number("set_speed_mm_s", lowest=0)
    row["set_angle_deg"] = number("set_angle_deg")
    row["static"] = number("static", integer=True, default=0, choices=(0, 1))
    row["roi"] = number("roi", integer=True, lowest=0, default=default_roi)
    if row["roi"] is not None and 0 < row["roi"] < 16:
        bad.append(f"roi '{raw['roi']}' must be 0 (largest square) or at least 16")
    row["exposure_s"] = number("exposure_s", above=0)
    row["ref_dt_s"] = number("ref_dt_s", above=0)
    if not row["static"] and row["set_speed_mm_s"] == 0:
        bad.append("set_speed_mm_s is 0 but static is 0; mark images of the stopped surface with static 1")

    row["path"], row["frame"] = None, None
    if raw["file"] == "":
        bad.append("file is empty")
    else:
        row["path"], row["frame"] = _split_frame(raw["file"], folder)
        bad += _check_path(row["path"], row["frame"], raw["file"])
    refs = [raw.get("ref_frame_a", ""), raw.get("ref_frame_b", "")]
    row["ref_frame_a"], row["ref_frame_b"] = refs
    row["ref"] = None
    if any(refs) and not all(refs):
        bad.append("give both ref_frame_a and ref_frame_b, or neither")
    elif all(refs):
        row["ref"] = [_split_frame(r, folder) for r in refs]
        for (p, k), r in zip(row["ref"], refs):
            bad += _check_path(p, k, r)
    return row, bad


def _rule(lowest, above, below, integer, choices):
    """What a value must be, for the message about a bad one."""
    if choices is not None:
        return " or ".join(str(c) for c in choices)
    kind = "an integer" if integer else "a number"
    if lowest is not None and below is not None:
        return f"{kind} in [{lowest:g}, {below:g})"
    if lowest is not None:
        return f"{kind} >= {lowest:g}"
    if above is not None:
        return f"{kind} > {above:g}"
    return kind


def _check_path(path: Path, frame, text: str):
    if not path.is_file():
        return [f"file not found: {text} ({path})"]
    if frame is not None and frame < 1:
        return [f"{text}: frames are counted from 1"]
    return []


# Images

def central(img, side: int):
    """Central ``side``-by-``side`` square, the same one that estimate_velocity analyses."""
    h, w = img.shape
    r0, c0 = (h - side) // 2, (w - side) // 2
    return img[r0:r0 + side, c0:c0 + side]


@lru_cache(maxsize=4)
def _read(path: str):
    """Array of an image or array file (fftvel CLI loader); its notes on stderr are dropped."""
    with contextlib.redirect_stderr(io.StringIO()):
        return _load(path)


def frames_of(path: Path, frame=None):
    """List of ``(label, 2-D array)``: the image, or one entry per frame of a stack.

    ``frame`` (from 1) selects one frame of a stack. Raises ValueError with a
    message if the file cannot be read or holds no grey image.
    """
    try:
        a = _read(str(path))
    except SystemExit as err:                   # the CLI loader exits with a message
        msg = str(err.code).replace("choose one with --var", "save only the image or the stack in the "
                                    "file, e.g. save('rec1.mat', 'recording')")
        raise ValueError(msg) from None
    shape = a.shape
    a = np.squeeze(a)
    if a.ndim == 2:
        if frame not in (None, 1):
            raise ValueError(f"holds one image ({' x '.join(map(str, shape))}), not frame {frame}")
        return [("", a)]
    if a.ndim != 3:
        raise ValueError(f"expected an image or a stack of images, got shape {shape}")
    try:
        axis = _frames_axis(a.shape, "auto")
    except SystemExit as err:
        raise ValueError(str(err.code)) from None
    nf = a.shape[axis]
    if frame is not None:
        if frame > nf:
            raise ValueError(f"has {nf} frames, not frame {frame}")
        return [(f"#{frame}", np.take(a, frame - 1, axis=axis))]
    return [(f"#{k + 1}", np.take(a, k, axis=axis)) for k in range(nf)]


def build_tasks(rows, clip_outliers: bool = False):
    """One task per image (a stack gives one per frame), and the result rows of files that failed.

    Each task holds a copy of the analysed square only, so that a long list
    of full camera frames does not stay in memory.
    """
    tasks, failed = [], []
    seen = {}                                   # (file, frame) -> line of its first listing
    for row in rows:
        try:
            images = frames_of(row["path"], row["frame"])
            ref = None
            if row["ref"] is not None:
                ref = [frames_of(p, k) for p, k in row["ref"]]
                if any(len(r) != 1 for r in ref):
                    raise ValueError("a reference frame is a stack; select one frame with file#k")
                ref = [r[0][1] for r in ref]
        except ValueError as err:
            failed.append(dict(_meta_columns(row, row["file"]), error=str(err)))
            continue
        for label, img in images:
            name = row["file"] + (label if row["frame"] is None else "")
            h, w = img.shape
            side = min(h, w) if row["roi"] == 0 else min(row["roi"], h, w)
            task = dict(row=row, name=name, image=central(img, side).copy(), clip=clip_outliers,
                        warnings=_image_warnings(row, img, side))
            key = (str(row["path"].resolve()), label or "#1")
            if key in seen:
                task["warnings"].append(f"the same image as in line {seen[key]}, so it counts twice in "
                                        f"the summary")
            seen.setdefault(key, row["line"])
            if ref is not None:
                s = min(side, *ref[0].shape, *ref[1].shape)
                task["ref"] = (central(ref[0], s).copy(), central(ref[1], s).copy())
            tasks.append(task)
    return tasks, failed


def _image_warnings(row, img, side):
    out = []
    K, T, duty = row["exposures"], row["flash_period_s"], row["duty"]
    if str(row["path"]).lower().endswith((".jpg", ".jpeg", ".jfif")):
        out.append("JPEG: its 8 x 8 block grid can give wrong results that are marked valid")
    if row["exposure_s"] is not None and T is not None and row["exposure_s"] < (K - 1 + duty) * T:
        out.append(f"exposure {row['exposure_s']:g} s is shorter than (K - 1 + duty) * T = "
                   f"{(K - 1 + duty) * T:g} s, so the image holds fewer than {K} flashes")
    crop = central(img, side)
    if crop.dtype.kind in "iu":
        share = float(np.mean(crop == crop.max()))
        if share > SATURATED:
            out.append(f"{100 * share:.1f} % of the analysed pixels are at the maximum {crop.max()} "
                       f"(saturation or glints; try --clip-outliers)")
    return out


# Estimates

def physical(est, gsd, period):
    """fftvel.to_physical of ``est`` (mm/s for gsd in mm/px), or None without gsd or period."""
    if gsd is None or period is None:
        return None
    return to_physical(est, gsd, period)


def lowest_speed(K: int, n: int) -> float:
    """Lowest speed in px per flash that estimate_velocity measures with K flashes on an n px square.

    The ``vmin`` of fftvel.core.fit_speed for the default padding, about
    3.35/K; below it the estimate comes out not in range.
    """
    npad = DEFAULT_PAD_FACTOR * n
    return VMIN_FACTOR * npad / (K * math.floor(RFRAC * (npad // 2 - 1)))


def _meta_columns(row, name):
    """Result row with the columns of meta.csv and empty result columns (every row has all FIELDS)."""
    out = dict.fromkeys(FIELDS)
    out.update(file=name, line=row["line"])
    out.update({k: row.get(k) for k in META})
    return out


def analyse(task):
    """Result row of one image. Top-level, so that it also runs in a worker process on Windows."""
    row = task["row"]
    K, T, gsd, duty = row["exposures"], row["flash_period_s"], row["gsd_mm_per_px"], row["duty"]
    res = _meta_columns(row, task["name"])
    warnings = list(task["warnings"])
    try:
        e = estimate_velocity(task["image"], L=max(1, (K - 1) // 2), exposures=K, duty=duty,
                              clip_outliers=task["clip"])
    except ValueError as err:
        res["error"], res["warnings"] = str(err), "; ".join(warnings)
        return res
    reason = []
    if e.quality < e.threshold:
        reason.append(f"quality {e.quality:.2f} < {e.threshold:.2f}")
    if not e.in_range:
        reason.append("out of range")
    if e.isotropic:
        reason.append("isotropic")
    res.update(roi=e.n, angle_deg=e.angle_deg, speed_px=e.speed, vx_px=e.vx, vy_px=e.vy,
               quality=e.quality, quality_perp=e.quality_perp, threshold=e.threshold, geba=e.geba,
               in_range=int(e.in_range), isotropic=int(e.isotropic), valid=int(e.valid),
               reason="; ".join(reason))
    mm = physical(e, gsd, T)
    res["speed_mm_s"] = None if mm is None else mm.speed
    if mm is not None:
        res["vx_mm_s"], res["vy_mm_s"] = mm.vx, mm.vy
    set_speed, set_angle = row["set_speed_mm_s"], row["set_angle_deg"]
    if not row["static"]:
        if set_speed and gsd is not None and T is not None:
            res["expected_px"] = set_speed * T / gsd
            lo, hi = lowest_speed(K, e.n), e.n / 4
            if not lo <= res["expected_px"] <= hi:
                warnings.append(f"the set speed is {res['expected_px']:.3g} px per flash, outside "
                                f"{lo:.2g} to {hi:g} px, the range measured for K = {K} on a {e.n} px "
                                f"square (about 3.35/K to n/4)")
        if set_speed and res["speed_mm_s"] is not None:
            res["speed_err"] = res["speed_mm_s"] / set_speed - 1
        if set_angle is not None:
            res["angle_err_deg"] = angle_difference(e.angle_deg, set_angle)
    if "ref" in task:
        dx, dy, peak = phase_correlation(*task["ref"])
        res.update(ref_dx_px=dx, ref_dy_px=dy, ref_peak=peak,
                   ref_angle_deg=float(np.rad2deg(np.arctan2(-dy, dx)) % 180))
        res["ref_angle_diff_deg"] = angle_difference(e.angle_deg, res["ref_angle_deg"])
        if np.hypot(dx, dy) > task["ref"][0].shape[0] / 4:
            warnings.append(f"the reference frames are {np.hypot(dx, dy):.0f} px apart, more than a "
                            f"quarter of the {task['ref'][0].shape[0]} px square; phase correlation "
                            f"may fail")
        if row["ref_dt_s"] is not None:
            disp = float(np.hypot(dx, dy))
            if T is not None:                   # displacement per flash period
                res["ref_speed_px"] = disp * T / row["ref_dt_s"]
                res["ref_speed_diff"] = e.speed / res["ref_speed_px"] - 1 if disp > 0 else None
            ref = physical(SimpleNamespace(speed=disp, vx=dx, vy=dy, angle_deg=res["ref_angle_deg"]),
                           gsd, row["ref_dt_s"])
            if ref is not None:
                res["ref_speed_mm_s"] = ref.speed
    res["warnings"] = "; ".join(warnings)
    return res


def run(tasks, jobs: int = 1, log=print):
    """Result rows of all tasks in their order; ``jobs > 1`` uses a process pool."""
    jobs = max(1, min(jobs, len(tasks)))
    if jobs == 1:
        results = []
        for t in tasks:
            results.append(analyse(t))
            log(line(results[-1]))
        return results
    results = []
    with Pool(jobs) as pool:
        for r in pool.imap(analyse, tasks, chunksize=1):     # analyse is top-level (Windows spawn)
            results.append(r)
            log(line(r))
    return results


def line(r) -> str:
    """One line of progress for a result row."""
    if r.get("error"):
        return f"{r['file']}: error: {r['error']}"
    speed = (f"{r['speed_mm_s']:.4g} mm/s ({r['speed_px']:.3f} px/flash)" if r.get("speed_mm_s") is not None
             else f"{r['speed_px']:.3f} px/flash")
    state = "valid" if r["valid"] else f"not valid: {r['reason']}"
    tag = " [static]" if r["static"] else ""
    return f"{r['file']}{tag}: {speed} at {r['angle_deg']:.1f} deg, quality {r['quality']:.2f}, {state}"


# Summary

def upper_bound(k: int, n: int, confidence: float = CONFIDENCE) -> float:
    """One-sided Clopper-Pearson upper confidence bound of a rate with ``k`` events in ``n`` trials.

    For ``k = 0`` it is ``1 - (1 - confidence)**(1/n)``, about ``3/n`` at 95 %.
    """
    if n <= 0:
        return math.nan
    if k >= n:
        return 1.0
    from scipy.stats import beta
    return float(beta.ppf(confidence, k + 1, n - k))


def _off(r, tol_speed, tol_angle) -> bool:
    return ((r.get("speed_err") is not None and abs(r["speed_err"]) > tol_speed)
            or (r.get("angle_err_deg") is not None and abs(r["angle_err_deg"]) > tol_angle))


def _badness(r, tol_speed, tol_angle) -> float:
    s = abs(r["speed_err"]) / tol_speed if r.get("speed_err") is not None else 0.0
    a = abs(r["angle_err_deg"]) / tol_angle if r.get("angle_err_deg") is not None else 0.0
    return max(s, a)


def summary(results, commented=(), tol_speed=TOL_SPEED, tol_angle=TOL_ANGLE, worst=5, log=print):
    """Print the summary of the result rows and return its numbers as a dict."""
    ok = [r for r in results if not r.get("error")]
    failed = [r for r in results if r.get("error")]
    moving = [r for r in ok if not r["static"]]
    static = [r for r in ok if r["static"]]
    valid = [r for r in moving if r["valid"]]
    checked = [r for r in valid if r.get("speed_err") is not None or r.get("angle_err_deg") is not None]
    off = [r for r in checked if _off(r, tol_speed, tol_angle)]
    static_valid = [r for r in static if r["valid"]]
    out = dict(images=len(results), analysed=len(ok), failed=len(failed), commented=len(commented),
               moving=len(moving), moving_valid=len(valid), static=len(static),
               static_valid=len(static_valid), static_bound=upper_bound(len(static_valid), len(static)),
               checked=len(checked), off=len(off), off_bound=upper_bound(len(off), len(checked)))
    conf = f"{100 * CONFIDENCE:.0f} % upper bound of the rate"
    log("")
    log(f"Images: {len(ok)} analysed ({len(moving)} moving, {len(static)} static), {len(failed)} not "
        f"analysed, {len(commented)} rows commented out")
    if moving:
        log(f"Moving images marked valid: {len(valid)} of {len(moving)} ({100 * len(valid) / len(moving):.0f} %)")
    else:
        log("Moving images: none")
    few = False
    for key, name, unit, scale, need in (
            ("speed_err", "Speed error |est/set - 1|", " %", 100, "a set speed, gsd and flash period"),
            ("angle_err_deg", "Angle error |est - set|", " deg", 1, "a set angle")):
        x = np.array([abs(r[key]) for r in valid if r.get(key) is not None]) * scale
        if x.size:
            p50, p95 = np.percentile(x, [50, 95])
            k = key.split("_")[0]
            out.update({f"{k}_median": float(p50), f"{k}_p95": float(p95), f"{k}_max": float(x.max()),
                        f"{k}_n": x.size})
            log(f"{name} of the valid images ({x.size}): median {p50:.3g}{unit}, 95th percentile "
                f"{p95:.3g}{unit}, max {x.max():.3g}{unit}" + (" *" if x.size < 20 else ""))
            few = few or x.size < 20
        else:
            log(f"{name}: no valid moving image with {need}")
    if checked:
        log(f"Valid but off by more than {100 * tol_speed:g} % or {tol_angle:g} deg: {len(off)} of "
            f"{len(checked)} ({conf}: {100 * out['off_bound']:.1f} %)")
    if static:
        log(f"Static images marked valid: {len(static_valid)} of {len(static)} (should be 0; "
            f"{conf}: {100 * out['static_bound']:.1f} %)")
    else:
        log("Static images: none. Record some of the stopped surface with the same settings; "
            "they must come out not valid")
    refs = [r for r in valid if r.get("ref_angle_diff_deg") is not None]
    if refs:
        d = np.array([abs(r["ref_angle_diff_deg"]) for r in refs])
        s = np.array([r["ref_speed_diff"] for r in refs if r.get("ref_speed_diff") is not None])
        out["ref_angle_median"] = float(np.median(d))
        text = (f"Phase correlation of the reference frames ({len(refs)} valid "
                f"image{'s' if len(refs) > 1 else ''}): angle difference median {np.median(d):.2g} deg")
        if s.size:
            out["ref_speed_median"] = 100 * float(np.median(s))
            text += f", speed est/ref - 1 median {out['ref_speed_median']:+.2g} %"
        log(text)
    if few:
        log("  * fewer than 20 values: the 95th percentile lies between the largest ones")
    bad = sorted(checked, key=lambda r: -_badness(r, tol_speed, tol_angle))[:worst]
    if static_valid or bad:
        log(f"Worst images (static ones marked valid, then valid ones by error relative to "
            f"{100 * tol_speed:g} % / {tol_angle:g} deg):")
        for r in static_valid:
            log(f"  {r['file']}  static but valid: {r['speed_px']:.3f} px/flash at {r['angle_deg']:.1f} deg, "
                f"quality {r['quality']:.2f}")
        for r in bad:
            parts = []
            if r.get("speed_err") is not None:
                parts.append(f"speed {100 * r['speed_err']:+.2f} %")
            if r.get("angle_err_deg") is not None:
                parts.append(f"angle {r['angle_err_deg']:+.2f} deg")
            log(f"  {r['file']}  {', '.join(parts)}, quality {r['quality']:.2f}")
    rejected = [r for r in moving if not r["valid"]]
    if rejected:
        log(f"Moving images not valid ({len(rejected)}" + (f", the first {worst}" if len(rejected) > worst else "")
            + "):")
        for r in rejected[:worst]:
            log(f"  {r['file']}  {r['reason']}")
    if failed:
        log(f"Not analysed ({len(failed)}; problems in meta.csv, files that cannot be read, estimate errors):")
    for r in failed:
        log(f"  {r['file']} (line {r['line']}): {r['error']}")
    for c in commented:
        log(f"  commented out: {c}")
    warned = [r for r in ok if r.get("warnings")]
    if warned:
        log(f"Warnings for {len(warned)} image{'s' if len(warned) > 1 else ''} (column warnings of the "
            f"results), e.g. {warned[0]['file']}: {warned[0]['warnings']}")
    return out


# Figure

def plot(results, path: Path, tol_speed=TOL_SPEED, tol_angle=TOL_ANGLE, label=""):
    """Estimated against set speed, the speed errors and the angle errors, written to ``path``.

    The axes are scaled to the set values and the valid estimates; not valid
    estimates outside them are left out and counted in the subtitle.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ok = [r for r in results if not r.get("error")]
    moving = [r for r in ok if not r["static"]]
    static = [r for r in ok if r["static"]]
    fig, (ax1, ax2, ax3) = common.figure(3, 13.5)
    hidden = set()
    pts = [r for r in moving if r.get("speed_err") is not None]
    if pts:
        set_v = np.array([r["set_speed_mm_s"] for r in pts])
        est_v = np.array([r["speed_mm_s"] for r in pts])
        err = 100 * np.array([r["speed_err"] for r in pts])
        valid = np.array([bool(r["valid"]) for r in pts])
        top = 1.15 * max(set_v.max(), est_v[valid].max() if valid.any() else 0)
        x = np.array([0, top])
        ax1.fill_between(x, x * (1 - tol_speed), x * (1 + tol_speed), color=AQUA, alpha=0.2, lw=0)
        ax1.plot(x, x, color=INK2, lw=1)
        _scatter(ax1, set_v, est_v, valid)
        ax1.set_xlim(0, top)
        ax1.set_ylim(0, top)
        lim = _limit(err[valid], 100 * tol_speed, 100)
        ax2.axhspan(-100 * tol_speed, 100 * tol_speed, color=AQUA, alpha=0.2, lw=0)
        ax2.axhline(0, color=INK2, lw=1)
        _scatter(ax2, set_v, err, valid)
        ax2.set_xlim(0, top)
        ax2.set_ylim(-lim, lim)
        hidden |= {r["file"] for r, e, y in zip(pts, err, est_v) if abs(e) > lim or y > top}
    else:
        for ax in (ax1, ax2):
            ax.text(0.5, 0.5, "no set speed with gsd and flash period", transform=ax.transAxes,
                    ha="center", color=INK2)
    common.style(ax1, "set speed (mm/s)", "estimated speed (mm/s)")
    common.style(ax2, "set speed (mm/s)", "speed error (%)")
    if pts:
        ax1.legend(frameon=False, fontsize=8, labelcolor=INK2, loc="upper left")

    ang = [r for r in moving if r.get("angle_err_deg") is not None]
    if ang:
        err = np.array([r["angle_err_deg"] for r in ang])
        valid = np.array([bool(r["valid"]) for r in ang])
        lim = _limit(err[valid], tol_angle, 90)
        ax3.axhspan(-tol_angle, tol_angle, color=AQUA, alpha=0.2, lw=0)
        ax3.axhline(0, color=INK2, lw=1)
        _scatter(ax3, [r["set_angle_deg"] % 180 for r in ang], err, valid)
        ax3.set_xlim(0, 180)
        ax3.set_ylim(-lim, lim)
        hidden |= {r["file"] for r, e in zip(ang, err) if abs(e) > lim}
    else:
        ax3.text(0.5, 0.5, "no set angle", transform=ax3.transAxes, ha="center", color=INK2)
    common.style(ax3, "set direction (deg, mod 180)", "angle error (deg)")
    nv = sum(r["valid"] for r in moving)
    sub = (f"{label}{nv} of {len(moving)} moving images valid, {sum(r['valid'] for r in static)} of "
           f"{len(static)} static ones (should be 0); filled = valid, open = not valid, band = "
           f"{100 * tol_speed:g} % / {tol_angle:g} deg")
    if hidden:
        sub += f"; {len(hidden)} not valid outside the axes"
    common.title(fig, "Real measurements: estimate against set motion", sub)
    path.parent.mkdir(parents=True, exist_ok=True)
    common.FIGURES = path.parent
    common.save(fig, path.name)
    plt.close(fig)


def _limit(values, tol, cap):
    """Symmetric axis limit: at least 3 tolerances, room for the given values, at most ``cap``."""
    big = 1.15 * float(np.max(np.abs(values))) if np.size(values) else 0.0
    return min(cap, max(3 * tol, big))


def _scatter(ax, x, y, valid):
    x, y, v = np.asarray(x, float), np.asarray(y, float), np.asarray(valid, bool)
    if v.any():
        ax.plot(x[v], y[v], "o", color=BLUE, ms=6, mec=SURFACE, mew=1.2, label="valid")
    if (~v).any():
        ax.plot(x[~v], y[~v], "o", color=ORANGE, ms=6, mfc="none", label="not valid")


# Command line

def write_results(path: Path, results) -> None:
    """results.csv: one row per image, the columns in the order of FIELDS, UTF-8.

    Not common.write_csv, which sorts the columns by name and writes in the
    locale encoding; on Windows (cp1252) that fails at the end of the run on
    a note or file name with a character such as a Greek letter.
    """
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
        w.writeheader()
        w.writerows(results)


def _cannot_write(path: Path):
    """None if ``path`` can be written, else the reason (on Windows: the file is open in Excel)."""
    existed = path.exists()
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "a", encoding="utf-8"):
            pass
    except OSError as err:
        return err.strerror or str(err)
    if not existed:
        path.unlink()
    return None


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("meta", type=Path, help="meta.csv, one row per image (see data/real/README.md)")
    p.add_argument("--out", type=Path, metavar="FILE",
                   help="results CSV (default: results.csv next to meta.csv; a folder: results.csv in it)")
    p.add_argument("--figure", type=Path, metavar="PNG", help="also plot estimated against set motion")
    p.add_argument("--roi", type=int, default=DEFAULT_ROI,
                   help=f"side of the analysed central square for rows without roi (default {DEFAULT_ROI}; "
                        f"0: the largest square). The run time grows with its square")
    p.add_argument("--clip-outliers", action="store_true",
                   help="set pixels more than 6 robust sigma from the median to the median (glints)")
    p.add_argument("--tol-speed", type=float, default=TOL_SPEED,
                   help=f"relative speed error up to which a valid estimate counts as correct (default {TOL_SPEED})")
    p.add_argument("--tol-angle", type=float, default=TOL_ANGLE,
                   help=f"angle error in deg up to which a valid estimate counts as correct (default {TOL_ANGLE:g})")
    p.add_argument("--worst", type=int, default=5, help="number of worst images listed (default 5)")
    p.add_argument("--jobs", type=int, default=1, help="worker processes (default 1)")
    a = p.parse_args(argv)
    if a.roi < 0 or 0 < a.roi < 16:
        p.error("--roi must be 0 (the largest square) or at least 16")
    if a.jobs < 1:
        p.error("--jobs must be at least 1")
    if not (a.tol_speed > 0 and a.tol_angle > 0) or a.worst < 0:
        p.error("--tol-speed and --tol-angle must be greater than 0, --worst at least 0")
    meta = a.meta.resolve()
    if not meta.is_file():
        p.error(f"no such file: {a.meta}")
    out = (a.out.resolve() if a.out is not None else meta.parent / "results.csv")
    if out.is_dir():
        out = out / "results.csv"
    if out == meta:
        p.error("--out must not be meta.csv itself")
    try:
        rows, commented, problems, unknown = read_meta(meta, a.roi)
    except MetaError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2
    if unknown:
        print(f"note: columns not used: {', '.join(unknown)} (known: {', '.join(COLUMNS)})", file=sys.stderr)
    if not rows:
        print(f"error: {a.meta} lists no image", file=sys.stderr)
        return 2
    reason = _cannot_write(out)
    if reason:
        print(f"error: cannot write {out} ({reason}). Close it if it is open in Excel, or pass "
              f"--out OTHER.csv", file=sys.stderr)
        return 2
    skipped = [r for r in rows if r["problems"]]
    if skipped:
        print(f"{len(skipped)} of {len(rows)} rows of {a.meta.name} skipped because of problems"
              + ("; the other rows are analysed:" if len(skipped) < len(rows) else ":"))
        for msg in problems:
            print(f"  {msg}")
    tasks, failed = build_tasks([r for r in rows if not r["problems"]], a.clip_outliers)
    failed += [dict(_meta_columns(r, r["file"]), error="row skipped: " + "; ".join(r["problems"]))
               for r in skipped]
    for r in failed:
        if not r["error"].startswith("row skipped"):
            print(line(r))
    if not tasks:
        print(f"error: no image of {a.meta} could be analysed; results.csv not written", file=sys.stderr)
        summary(failed, commented, a.tol_speed, a.tol_angle, a.worst)
        return 1
    results = run(tasks, a.jobs) + failed
    results.sort(key=lambda r: r["line"])
    try:
        write_results(out, results)
    except OSError as err:                      # e.g. opened in Excel during the run
        alt = out.with_name(f"{out.stem}_{time.strftime('%Y%m%d_%H%M%S')}{out.suffix}")
        print(f"cannot write {out} ({err.strerror or err}); writing {alt.name} instead")
        out = alt
        write_results(out, results)
    print(f"results: {out}")
    summary(results, commented, a.tol_speed, a.tol_angle, a.worst)
    if a.figure is not None:
        plot(results, a.figure.resolve(), a.tol_speed, a.tol_angle, f"{meta.parent.name}: ")
        print(f"figure:  {a.figure.resolve()}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
