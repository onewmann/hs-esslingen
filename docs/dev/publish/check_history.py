"""Fail if a blob of the rewritten history still holds photo or video metadata, a local user path or a claim about 2025.

Run inside the rewritten repository: python check_history.py. Exit status 0: clean, 3: a blob needs
cleaning (listed on stderr); anything else means the check itself failed.
"""
import re
import struct
import subprocess
import sys

CONTAINERS = {b"moov", b"trak", b"mdia", b"minf", b"stbl", b"edts", b"udta"}


def boxes(data, start=0, end=None):
    """(type, start of the payload, end) of the MP4/MOV boxes in data[start:end]; the rest that is no box: b"tail"."""
    end = len(data) if end is None else end
    p = start
    while p + 8 <= end:
        n, kind = struct.unpack(">I4s", data[p:p + 8])
        h = 8
        if n == 1 and p + 16 <= end:
            n, h = struct.unpack(">Q", data[p + 8:p + 16])[0], 16
        elif n == 0:
            n = end - p
        if n < h or p + n > end:
            yield b"tail", p, end
            return
        yield kind, p + h, p + n
        p += n


def video_problems(data):
    """Sound tracks and location, device and capture-time metadata of an MP4/MOV file; mdat (the frames) is skipped."""
    meta, times = b"", []

    def visit(start, end):
        nonlocal meta
        for kind, a, b in boxes(data, start, end):
            if kind == b"mdat":
                continue
            if kind in CONTAINERS:
                visit(a, b)
                continue
            meta += kind + data[a:b]
            if kind in (b"mvhd", b"tkhd", b"mdhd") and b - a >= 20:
                wide = data[a] == 1                  # version 1: 64-bit times
                times.extend(struct.unpack(">QQ", data[a + 4:a + 20]) if wide else
                             struct.unpack(">II", data[a + 4:a + 12]))

    visit(0, len(data))
    out = []
    if re.search(rb"hdlr.{8}soun", meta, re.S):
        out.append("a sound track")
    for atom, what in ((b"\xa9xyz", "a location atom (\xa9xyz)"), (b"loci", "a location atom (loci)"),
                       (b"com.apple.quicktime.location", "a location key"), (b"location.ISO6709", "a location key"),
                       (b"com.android.", "an Android key"), (b"com.samsung.", "a Samsung key"),
                       (b"\xa9mak", "the maker of the device (\xa9mak)"), (b"\xa9mod", "the device model (\xa9mod)"),
                       (b"com.apple.quicktime.make", "the maker of the device"),
                       (b"com.apple.quicktime.model", "the device model")):
        if atom in meta:
            out.append(what)
    if any(times):
        out.append("a capture time in mvhd, tkhd or mdhd")
    return out


objs = subprocess.run(["git", "cat-file", "--batch-all-objects", "--batch-check=%(objectname) %(objecttype)"],
                      capture_output=True, text=True, check=True).stdout.split("\n")
blobs = [line.split()[0] for line in objs if line.endswith(" blob")]
cat = subprocess.Popen(["git", "cat-file", "--batch"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
bad = []
for sha in blobs:
    cat.stdin.write(sha.encode() + b"\n")
    cat.stdin.flush()
    size = int(cat.stdout.readline().split()[2])
    data = cat.stdout.read(size)
    cat.stdout.read(1)
    if data[:8] == b"\x89PNG\r\n\x1a\n":
        p = 8
        while p + 12 <= len(data):
            n = struct.unpack(">I", data[p:p + 4])[0]
            kind = data[p + 4:p + 8]
            if kind == b"eXIf" or (kind == b"iTXt" and data[p + 8:p + 25] == b"XML:com.adobe.xmp"):
                bad.append(f"{sha}: PNG with an {kind.decode()} chunk (photo metadata)")
            p += 12 + n
    if data[4:8] == b"ftyp":                     # MP4/MOV of a phone
        bad += [f"{sha}: video with {what}" for what in video_problems(data)]
    low = data.lower()
    for needle, what in ((b"gpslatitude", "GPS position"), (b"\\users\\", "a Windows user folder"),
                         (b"without ai", "a claim about the time before the rework"),
                         (b"my own work", "a claim about the time before the rework"),
                         (b"rework only", "a claim about the time before the rework")):
        if needle in low:
            bad.append(f"{sha}: {what} ({needle.decode()})")
cat.stdin.close()
cat.wait()
for b in bad:
    print(f"error: {b}", file=sys.stderr)
sys.exit(3 if bad else 0)
