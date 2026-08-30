#!/usr/bin/env python3
from pathlib import Path
import re
import sys

roots = [Path("src"), Path("test"), Path("example")]
files = sorted(p for root in roots for p in root.rglob("*.f90"))
problems = []
texts = {p: p.read_text(encoding="utf-8") for p in files}

iso_uses = sum(len(re.findall(r"^\s*use(?:\s*,\s*intrinsic\s*::|\s+)\s*iso_fortran_env\b", t, re.I | re.M)) for t in texts.values())
if iso_uses != 1:
    problems.append("iso_fortran_env should be imported in exactly one maintained source file")
if sum(len(re.findall(r"\bdp\s*=\s*real64\b", t, re.I)) for t in texts.values()) != 1:
    problems.append("dp = real64 should be defined exactly once")

patterns = {
    "double precision": re.compile(r"\bdouble\s+precision\b", re.I),
    "real*8": re.compile(r"\breal\s*\*\s*8\b", re.I),
    "kind(0.0d0)": re.compile(r"kind\s*\(\s*0\.0d0\s*\)", re.I),
    "D exponent": re.compile(r"(?<![A-Za-z0-9_])(?:\d+(?:\.\d*)?|\.\d+)[dD][+-]?\d+"),
}
for p, text in texts.items():
    for label, rx in patterns.items():
        if rx.search(text):
            problems.append(f"{p}: forbidden {label}")
    for n, line in enumerate(text.splitlines(), 1):
        if len(line) > 132:
            problems.append(f"{p}:{n}: line exceeds 132 columns")

if problems:
    print("source hygiene: FAILED")
    print("\n".join(problems))
    sys.exit(1)
print(f"source hygiene: ok ({len(files)} maintained Fortran files)")
