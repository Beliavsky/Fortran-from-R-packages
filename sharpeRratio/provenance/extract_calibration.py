#!/usr/bin/env python3
"""Document the numeric-vector extraction from sharpeRratio R data files.

R serialization uses big-endian XDR. At the offsets below, each object is a
REALSXP encoded as a four-byte flags word, a four-byte length, and binary64
values. These are the x, y, b, c, and d arrays in the three splinefun closure
environments. The production Fortran source stores their exact bit patterns.
"""
from pathlib import Path
import bz2
import gzip
import struct

ROOT = Path(__file__).resolve().parents[1]
UPSTREAM = ROOT / "original" / "sharpeRratio-master"

def vector(data: bytes, offset: int) -> list[float]:
    flags = struct.unpack(">I", data[offset:offset + 4])[0]
    if flags & 0xFF != 14:
        raise ValueError(f"offset {offset} is not a REALSXP")
    length = struct.unpack(">i", data[offset + 4:offset + 8])[0]
    end = offset + 8 + 8 * length
    return list(struct.unpack(f">{length}d", data[offset + 8:end]))

a_data = gzip.open(UPSTREAM / "data" / "a.rda", "rb").read()
f_data = gzip.open(UPSTREAM / "data" / "f.rda", "rb").read()
sysdata = bz2.open(UPSTREAM / "R" / "sysdata.rda", "rb").read()

OFFSETS = {
    "a": (a_data, [135, 807, 1479, 2151, 2823]),
    "a_medium": (sysdata, [11776, 14992, 18208, 21424, 24640]),
    "f": (f_data, [135, 383, 631, 879, 1127]),
}

for name, (data, offsets) in OFFSETS.items():
    arrays = [vector(data, offset) for offset in offsets]
    print(name, [len(array) for array in arrays])
    print("  x range:", arrays[0][0], arrays[0][-1])
    print("  y range:", min(arrays[1]), max(arrays[1]))
