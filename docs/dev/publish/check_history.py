"""Fail if any blob of the rewritten history still holds photo metadata, a local user path or a claim about 2025.

Run inside the rewritten repository: python check_history.py. Exit status 0: clean, 3: a blob needs
cleaning (listed on stderr); anything else means the check itself failed.
"""
import re
import struct
import subprocess
import sys

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
    if data[4:8] == b"ftyp":                     # MP4/MOV: sound tracks, and location atoms (GPS) of phones
        if re.search(rb"hdlr.{8}soun", data, re.S):
            bad.append(f"{sha}: video with a sound track")
        for atom, what in ((b"\xa9xyz", "a location atom (\xa9xyz)"), (b"loci", "a location atom (loci)"),
                           (b"com.apple.quicktime.location", "a location key"), (b"location.ISO6709", "a location key")):
            if atom in data:
                bad.append(f"{sha}: video with {what}")
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
