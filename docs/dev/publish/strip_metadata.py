# git filter-repo --blob-callback body: drop the EXIF and XMP chunks (location, camera, time) from PNG
# files; the image data stay byte for byte. Only 2D-FFT/data/image001.png at the tag has such chunks.
import struct
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
