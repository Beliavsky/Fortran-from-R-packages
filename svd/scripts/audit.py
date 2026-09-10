#!/usr/bin/env python3
"""Static release audit for the svd Fortran translation."""
from __future__ import annotations

import hashlib
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ERRORS: list[str] = []


def fail(message: str) -> None:
    ERRORS.append(message)


fortran_files = sorted(ROOT.rglob("*.f90"))
for path in fortran_files:
    raw = path.read_bytes()
    try:
        text = raw.decode("ascii")
    except UnicodeDecodeError:
        fail(f"non-ASCII Fortran: {path.relative_to(ROOT)}")
        continue
    if "implicit none" not in text.lower():
        fail(f"missing implicit none: {path.relative_to(ROOT)}")
    for lineno, line in enumerate(text.splitlines(), 1):
        if len(line) > 132:
            fail(f"line >132: {path.relative_to(ROOT)}:{lineno}")
        code = line.split("!", 1)[0]
        if ";" in code:
            fail(f"semicolon-packed statement: {path.relative_to(ROOT)}:{lineno}")
        if re.search(r"\b([A-Za-z_]\w*)\s*/=\s*\1\b", code):
            fail(f"self-comparison NaN idiom: {path.relative_to(ROOT)}:{lineno}")
    for pattern in (
        r"\bdouble\s+precision\b",
        r"\breal\s*\*\s*8\b",
        r"kind\s*\(\s*0\.0d0\s*\)",
        r"\d(?:\.\d*)?[dD][+-]?\d+",
    ):
        if re.search(pattern, text, re.IGNORECASE):
            fail(f"disallowed real-kind form: {path.relative_to(ROOT)}")

# Every dummy in maintained Fortran must have one declaration, INTENT/VALUE, and FORD comment.
procedure_header = re.compile(
    r"\s*(?:pure\s+|elemental\s+|pure\s+elemental\s+)?"
    r"(?:[\w()=*,: ]+\s+)?(subroutine|function)\s+(\w+)\s*\((.*)",
    re.IGNORECASE,
)
for path in fortran_files:
    lines = path.read_text(encoding="ascii").splitlines()
    i = 0
    while i < len(lines):
        match = procedure_header.match(lines[i])
        if not match:
            i += 1
            continue
        kind, name, tail = match.groups()
        header = tail
        j = i
        while ")" not in header and j + 1 < len(lines):
            j += 1
            header += " " + lines[j].strip().lstrip("&").rstrip("&")
        args = [part.strip() for part in header.split(")", 1)[0].split(",") if part.strip()]
        end = j + 1
        end_pattern = re.compile(rf"\s*end\s+{kind}(?:\s+{re.escape(name)})?\s*$", re.IGNORECASE)
        while end < len(lines) and not end_pattern.match(lines[end]):
            end += 1
        block = lines[j + 1 : end]
        lowered_args = [arg.lower() for arg in args]
        for arg in args:
            if not re.fullmatch(r"\w+", arg):
                fail(f"cannot parse dummy {arg!r}: {path.relative_to(ROOT)}:{i + 1}")
                continue
            declarations: list[tuple[int, str]] = []
            for lineno, declaration in enumerate(block, j + 2):
                if "::" not in declaration:
                    continue
                rhs = declaration.split("::", 1)[1].split("!!", 1)[0]
                names = [re.sub(r"\(.*", "", item.strip()).strip().lower() for item in rhs.split(",")]
                if arg.lower() in names:
                    declarations.append((lineno, declaration))
            if len(declarations) != 1:
                fail(f"dummy declaration count {len(declarations)} for {name}:{arg} in {path.relative_to(ROOT)}")
                continue
            lineno, declaration = declarations[0]
            lhs = declaration.split("::", 1)[0].lower()
            if "intent(" not in lhs and "value" not in lhs:
                fail(f"dummy lacks INTENT/VALUE: {path.relative_to(ROOT)}:{lineno}")
            if "!!" not in declaration:
                fail(f"dummy lacks FORD comment: {path.relative_to(ROOT)}:{lineno}")
            rhs = declaration.split("::", 1)[1].split("!!", 1)[0]
            declared = [re.sub(r"\(.*", "", item.strip()).strip().lower() for item in rhs.split(",")]
            if sum(item in lowered_args for item in declared) != 1:
                fail(f"multiple dummies on one declaration: {path.relative_to(ROOT)}:{lineno}")
        i = max(i + 1, end + 1)

# Duplicate maintained Fortran sources by byte content.
seen: dict[str, Path] = {}
for path in fortran_files:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest in seen:
        fail(f"duplicate Fortran source: {seen[digest].relative_to(ROOT)} and {path.relative_to(ROOT)}")
    seen[digest] = path

# Coverage accounting.
with (ROOT / "fpm.toml").open("rb") as handle:
    manifest = tomllib.load(handle)
translation = manifest["extra"]["translation"]
entries = translation.get("function", [])
total = translation["r_functions_total"]
mapped = translation["r_functions_translated"]
fraction = translation["frac_functions_translated"]
untranslated = translation["untranslated_r_functions"]
if mapped != len(entries):
    fail("coverage mapping-entry count disagrees with r_functions_translated")
if mapped + len(untranslated) != total:
    fail("coverage mapped + untranslated does not equal total")
if abs(fraction - mapped / total) > 1.0e-12:
    fail("coverage fraction is inconsistent")
readme = (ROOT / "README.md").read_text(encoding="utf-8")
if f"**{mapped} of {total} ({100.0 * mapped / total:.0f}%)**" not in readme:
    fail("README coverage headline disagrees with fpm.toml")
for entry in entries:
    if f"`{entry['r_name']}`" not in readme:
        fail(f"README coverage table missing {entry['r_name']}")

# Dependency sources and transient products must not be shipped.
forbidden_suffixes = {".o", ".obj", ".mod", ".smod", ".exe", ".a", ".so", ".dll", ".pyc", ".zip"}
for path in ROOT.rglob("*"):
    if not path.is_file():
        continue
    if path.suffix.lower() in forbidden_suffixes:
        fail(f"forbidden release product: {path.relative_to(ROOT)}")
for forbidden in ("rfortran_arpack", "rspectra_external", "la_lapack", "propack", "trlan"):
    for path in (ROOT / "src").glob("*.f90"):
        if forbidden in path.read_text(encoding="ascii").lower():
            if forbidden in {"propack", "trlan"}:
                continue
            fail(f"possible vendored dependency implementation in {path.relative_to(ROOT)}: {forbidden}")

if ERRORS:
    for error in ERRORS:
        print(f"ERROR: {error}")
    print(f"audit: FAIL ({len(ERRORS)} errors)")
    sys.exit(1)
print(f"audit: PASS ({len(fortran_files)} Fortran files, {mapped}/{total} mappings)")
