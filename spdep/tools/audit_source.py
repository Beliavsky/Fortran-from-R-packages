#!/usr/bin/env python3
"""Static source audit for the maintained spdep Fortran translation."""
from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FORTRAN = sorted(ROOT.rglob("*.f90"))
errors: list[str] = []

# Duplicate file contents.
seen_hash: dict[str, Path] = {}
for path in FORTRAN:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest in seen_hash:
        errors.append(f"duplicate Fortran source: {path} == {seen_hash[digest]}")
    else:
        seen_hash[digest] = path

# Content rules.
self_compare = re.compile(r"\b([A-Za-z_]\w*)\s*/=\s*\1\b", re.IGNORECASE)
d_exp = re.compile(r"(?<![A-Za-z_])(?:\d+(?:\.\d*)?|\.\d+)[dD][+-]?\d+")
for path in FORTRAN:
    for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        code = raw.split("!", 1)[0]
        if ";" in code:
            errors.append(f"{path}:{lineno}: semicolon-separated statement")
        if self_compare.search(code):
            errors.append(f"{path}:{lineno}: self-comparison NaN idiom")
        lower = code.lower()
        if "double precision" in lower or "real*8" in lower or "kind(0.0d0)" in lower:
            errors.append(f"{path}:{lineno}: disallowed real-kind spelling")
        if d_exp.search(code):
            errors.append(f"{path}:{lineno}: D-exponent real literal")
        if "intent(" in lower and "!!" not in raw:
            errors.append(f"{path}:{lineno}: dummy declaration lacks trailing FORD !! comment")

# Every dummy name in a routine header must appear on its own INTENT/VALUE declaration.
# Parse physical source lines so comments and continuation syntax do not confuse the check.
routine_re = re.compile(
    r"(?i)\\b(?:subroutine|function)\\s+(\\w+)\\s*\\((.*?)\\)"
)
end_re_template = r"(?i)^\\s*end\\s+(?:subroutine|function)\\s+{name}\\b"

for path in FORTRAN:
    lines = path.read_text(encoding="utf-8").splitlines()
    i = 0
    while i < len(lines):
        raw = lines[i]
        code = raw.split("!", 1)[0].strip()
        if re.search(r"(?i)\\b(?:subroutine|function)\\b", code) and not re.match(r"(?i)^\\s*end\\b", code):
            header = code
            j = i
            while header.rstrip().endswith("&") and j + 1 < len(lines):
                header = header.rstrip()[:-1] + " "
                j += 1
                continuation = lines[j].split("!", 1)[0].strip()
                if continuation.startswith("&"):
                    continuation = continuation[1:].lstrip()
                header += continuation
            match = routine_re.search(header)
            if match:
                name = match.group(1)
                args = [a.strip().lower() for a in match.group(2).split(",") if a.strip()]
                end_re = re.compile(end_re_template.format(name=re.escape(name)))
                k = j + 1
                while k < len(lines) and not end_re.match(lines[k].split("!", 1)[0]):
                    k += 1
                body = lines[j + 1:k]
                for arg in args:
                    decls = []
                    arg_re = re.compile(rf"(?i)::\\s*{re.escape(arg)}\\b")
                    for body_line in body:
                        if arg_re.search(body_line.split("!", 1)[0]):
                            decls.append(body_line)
                    if len(decls) != 1:
                        errors.append(
                            f"{path}: {name}: dummy {arg!r} is not declared exactly once on its own line"
                        )
                        continue
                    decl = decls[0]
                    low = decl.lower()
                    if "intent(" not in low and not re.search(r"\\bvalue\\b", low):
                        errors.append(f"{path}: {name}: dummy {arg!r} lacks INTENT or VALUE")
                    if "!!" not in decl:
                        errors.append(f"{path}: {name}: dummy {arg!r} lacks FORD documentation")

                    # A dummy declaration may contain array bounds but must not declare another entity.
                    entity_text = decl.split("::", 1)[1].split("!", 1)[0].strip()
                    depth = 0
                    has_top_level_comma = False
                    for ch in entity_text:
                        if ch == "(":
                            depth += 1
                        elif ch == ")":
                            depth = max(0, depth - 1)
                        elif ch == "," and depth == 0:
                            has_top_level_comma = True
                            break
                    if has_top_level_comma:
                        errors.append(f"{path}: {name}: dummy {arg!r} shares a declaration line")
                i = k
        i += 1

# Build products or nested archives must not be shipped.
for path in ROOT.rglob("*"):
    if not path.is_file():
        continue
    if path.suffix.lower() in {".o", ".obj", ".mod", ".smod", ".exe", ".zip", ".pyc"}:
        errors.append(f"disallowed build/archive artifact: {path}")
    if "__pycache__" in path.parts or "build" in path.parts:
        errors.append(f"disallowed cache/build path: {path}")

if errors:
    print("SOURCE AUDIT FAILED", file=sys.stderr)
    for error in errors:
        print(f"- {error}", file=sys.stderr)
    raise SystemExit(1)

print(f"SOURCE AUDIT PASS: {len(FORTRAN)} Fortran files checked")
