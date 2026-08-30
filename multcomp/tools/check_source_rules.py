#!/usr/bin/env python3
"""Check maintained multcomp Fortran source conventions."""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES = sorted((ROOT / "src").glob("*.f90")) + sorted((ROOT / "test").glob("*.f90")) + sorted((ROOT / "example").glob("*.f90"))

PROC_RE = re.compile(r"\b(?:subroutine|function)\s+([a-z]\w*)\s*\((.*?)\)", re.I)
DECL_RE = re.compile(r"^\s*[^!]*\b(?:intent\s*\([^)]*\)|\bvalue\b)[^!]*::\s*([^!]+?)(?:\s*!!(.*))?$", re.I)
DUMMY_DECL_RE = re.compile(r"^\s*[^!]*::\s*([^!]+?)(?:\s*!!(.*))?$", re.I)


def logical_lines(lines: list[str]) -> list[tuple[int, str]]:
    out: list[tuple[int, str]] = []
    buf = ""
    start = 0
    continuing = False
    for no, raw in enumerate(lines, 1):
        code = raw.split("!", 1)[0].rstrip()
        if not continuing:
            buf = code
            start = no
        else:
            piece = code.lstrip()
            if piece.startswith("&"):
                piece = piece[1:].lstrip()
            buf += " " + piece
        continuing = code.endswith("&")
        if continuing:
            buf = buf[:-1].rstrip()
        elif buf.strip():
            out.append((start, buf))
            buf = ""
    if buf.strip():
        out.append((start, buf))
    return out


def split_names(text: str) -> list[str]:
    names: list[str] = []
    depth = 0
    current = []
    for char in text:
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        if char == "," and depth == 0:
            token = "".join(current).strip()
            if token:
                names.append(token)
            current = []
        else:
            current.append(char)
    token = "".join(current).strip()
    if token:
        names.append(token)
    return names


def base_name(decl: str) -> str:
    decl = decl.split("=", 1)[0].strip()
    return re.split(r"\s*\(", decl, maxsplit=1)[0].strip().lower()


errors: list[str] = []
for path in FILES:
    lines = path.read_text(encoding="utf-8").splitlines()
    for no, line in enumerate(lines, 1):
        if len(line) > 132:
            errors.append(f"{path.relative_to(ROOT)}:{no}: line exceeds 132 columns ({len(line)})")
        lower = line.lower()
        if "double precision" in lower or re.search(r"\breal\s*\*\s*8\b", lower):
            errors.append(f"{path.relative_to(ROOT)}:{no}: forbidden real declaration")
        if re.search(r"\b\d+(?:\.\d*)?[dD][+-]?\d+\b", line):
            errors.append(f"{path.relative_to(ROOT)}:{no}: D-exponent literal is forbidden")
        if "real64" in lower and path.name != "multcomp_kinds.f90":
            errors.append(f"{path.relative_to(ROOT)}:{no}: real64 may only appear in multcomp_kinds")

    dummy_sets: list[tuple[int, set[str], set[str]]] = []
    for no, logical in logical_lines(lines):
        match = PROC_RE.search(logical)
        if match:
            args = {x.strip().lower() for x in split_names(match.group(2)) if x.strip()}
            dummy_sets.append((no, args, set()))

    # Assign physical declaration lines to the innermost procedure whose declaration precedes them.
    for no, line in enumerate(lines, 1):
        decl = DUMMY_DECL_RE.match(line)
        if not decl:
            continue
        names = split_names(decl.group(1))
        if not names:
            continue
        active = None
        for item in dummy_sets:
            if item[0] <= no:
                active = item
            else:
                break
        if active is None:
            continue
        dummy_names = active[1]
        declared_dummies = [base_name(name) for name in names if base_name(name) in dummy_names]
        if not declared_dummies:
            continue
        if len(names) != 1:
            errors.append(f"{path.relative_to(ROOT)}:{no}: dummy arguments must be declared one per line")
        intent_decl = DECL_RE.match(line)
        if intent_decl is None:
            errors.append(f"{path.relative_to(ROOT)}:{no}: dummy declaration lacks explicit INTENT or VALUE")
        comment = decl.group(2)
        if comment is None or len(comment.strip()) < 12:
            errors.append(f"{path.relative_to(ROOT)}:{no}: dummy declaration lacks meaningful trailing FORD !! comment")
        active[2].update(declared_dummies)

    for no, dummies, declared in dummy_sets:
        missing = sorted(dummies - declared)
        if missing:
            errors.append(f"{path.relative_to(ROOT)}:{no}: undeclared/documentation-missing dummies: {', '.join(missing)}")

if errors:
    print("Fortran source-rule audit FAILED", file=sys.stderr)
    for error in errors:
        print(error, file=sys.stderr)
    raise SystemExit(1)

print(f"Fortran source-rule audit passed for {len(FILES)} maintained files")
