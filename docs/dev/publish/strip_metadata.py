# git filter-repo --blob-callback body.
# 1. Videos: a blob whose id names a file in clean_videos/ (next to the rewritten repository) is replaced by
#    it: the same video stream without the sound track and the file metadata (made with
#    ffmpeg -map 0:v:0 -c copy -map_metadata -1; the frames are identical). These are the two phone videos
#    of 2025 at the tag.
# 2. PNG files: drop the EXIF and XMP chunks (location, camera, time); the image data stay byte for byte.
#    Only 2D-FFT/data/image001.png at the tag has such chunks.
import os
import struct
clean = os.path.join("..", "clean_videos", blob.original_id.decode() + ".mp4")
if os.path.isfile(clean):
    with open(clean, "rb") as f:
        blob.data = f.read()
d = blob.data
if d[:8] == b"\x89PNG\r\n\x1a\n" and (b"eXIf" in d or b"XML:com.adobe.xmp" in d):
    out, p = [d[:8]], 8
    while p + 12 <= len(d):
        n = struct.unpack(">I", d[p:p + 4])[0]
        kind, chunk = d[p + 4:p + 8], d[p:p + 12 + n]
        if not (kind == b"eXIf" or (kind == b"iTXt" and chunk[8:25] == b"XML:com.adobe.xmp")):
            out.append(chunk)
        p += 12 + n
    blob.data = b"".join(out)
