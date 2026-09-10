#!/usr/bin/env python3
from __future__ import annotations
import hashlib
import re
import sys
import tomllib
from pathlib import Path

root = Path(__file__).resolve().parents[1]
errors: list[str] = []

with (root / "fpm.toml").open("rb") as fh:
    manifest = tomllib.load(fh)
tr = manifest["extra"]["translation"]
entries = tr.get("function", [])
if tr["r_functions_total"] != 122:
    errors.append("r_functions_total is not 122")
if tr["r_functions_translated"] != 122 or len(entries) != 122:
    errors.append("translation mapping count is not 122")
if len(tr["untranslated_r_functions"]) != 0:
    errors.append("untranslated function list is not empty")
if abs(tr["frac_functions_translated"] - 1.0) > 1.0e-14:
    errors.append("translation fraction is inconsistent")

fortran = sorted(list((root / "src").glob("*.f90")) + list((root / "test").glob("*.f90")) +
                 list((root / "app").glob("*.f90")) + list((root / "example").glob("*.f90")))
seen: dict[str, Path] = {}
for path in fortran:
    data = path.read_bytes()
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError:
        errors.append(f"non-ASCII Fortran: {path.relative_to(root)}")
        continue
    digest = hashlib.sha256(data).hexdigest()
    if digest in seen:
        errors.append(f"duplicate Fortran content: {path.relative_to(root)} and {seen[digest].relative_to(root)}")
    seen[digest] = path
    for lineno, line in enumerate(text.splitlines(), 1):
        if len(line) > 132:
            errors.append(f"line >132: {path.relative_to(root)}:{lineno}")
        code = line.split("!", 1)[0]
        if ";" in code:
            errors.append(f"semicolon statement: {path.relative_to(root)}:{lineno}")
        lower = code.lower()
        if re.search(r"\bdouble\s+precision\b|\breal\s*\*\s*8\b|kind\s*\(\s*0\.0d0\s*\)", lower):
            errors.append(f"disallowed real kind: {path.relative_to(root)}:{lineno}")
        if re.search(r"\d(?:\.\d*)?d[+-]?\d", lower):
            errors.append(f"D-exponent literal: {path.relative_to(root)}:{lineno}")
        if re.search(r"\b([a-z_]\w*)\s*/=\s*\1\b", lower):
            errors.append(f"self-comparison NaN idiom: {path.relative_to(root)}:{lineno}")
    if "implicit none" not in text.lower():
        errors.append(f"missing implicit none: {path.relative_to(root)}")

bad_suffixes = {".o", ".obj", ".mod", ".smod", ".exe", ".a", ".so", ".dll", ".zip"}
for path in root.rglob("*"):
    if path.is_file() and path.suffix.lower() in bad_suffixes:
        errors.append(f"build/archive artifact in package: {path.relative_to(root)}")
    if path.is_dir() and path.name.lower() in {"build", "__pycache__", ".pytest_cache"}:
        errors.append(f"cache/build directory in package: {path.relative_to(root)}")

readme = (root / "README.md").read_text(encoding="utf-8")
if "122 of 122 (100.0%)" not in readme:
    errors.append("README coverage summary is inconsistent")

if errors:
    for error in errors:
        print("audit: ERROR:", error)
    sys.exit(1)
print(f"audit: PASS ({len(fortran)} Fortran files, {len(entries)} mapping entries)")
