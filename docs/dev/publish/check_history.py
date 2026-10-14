"""Fail if any blob of the rewritten history still holds photo metadata, a local user path or a claim about 2025.

Run inside the rewritten repository: python check_history.py
"""
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
    low = data.lower()
    for needle, what in ((b"gpslatitude", "GPS position"), (b"\\users\\", "a Windows user folder"),
                         (b"without ai", "a claim about the time before the rework"),
                         (b"my own work", "a claim about the time before the rework")):
        if needle in low:
            bad.append(f"{sha}: {what} ({needle.decode()})")
cat.stdin.close()
cat.wait()
for b in bad:
    print(f"error: {b}", file=sys.stderr)
sys.exit(1 if bad else 0)
