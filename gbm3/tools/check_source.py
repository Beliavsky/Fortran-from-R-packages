#!/usr/bin/env python3
from pathlib import Path
import re
import sys

root = Path(__file__).resolve().parents[1]
files = sorted(root.glob("src/*.f90")) + sorted(root.glob("test/*.f90")) + sorted(root.glob("example/*.f90"))
text = "\n".join(path.read_text(encoding="utf-8") for path in files)
errors = []

for pattern, label in [
    (r"\bdouble\s+precision\b", "double precision"),
    (r"\breal\s*\*\s*8\b", "real*8"),
    (r"kind\s*\(\s*0\.0[dD]0\s*\)", "kind(0.0d0)"),
    (r"(?<![A-Za-z0-9_])(?:\d+(?:\.\d*)?|\.\d+)[dD][+-]?\d+", "D exponent literal"),
]:
    if re.search(pattern, text, flags=re.IGNORECASE):
        errors.append(f"forbidden maintained-source construct: {label}")

real64_files = [path for path in files if re.search(r"\breal64\b", path.read_text(encoding="utf-8"))]
expected = root / "src" / "gbm3_kinds.f90"
if real64_files != [expected]:
    errors.append("real64 must appear only in src/gbm3_kinds.f90")

for path in files:
    data = path.read_text(encoding="utf-8")
    if not data.isascii():
        errors.append(f"non-ASCII maintained source: {path.relative_to(root)}")

if errors:
    print("source hygiene check failed:")
    for error in errors:
        print(f"- {error}")
    sys.exit(1)

print(f"source hygiene check passed ({len(files)} maintained Fortran files)")
