#!/usr/bin/env python3
"""Audit maintained classInt Fortran sources against repository translation rules."""
from pathlib import Path
import re
import sys

roots = [Path("src"), Path("test"), Path("example")]
files = sorted(p for root in roots if root.exists() for p in root.rglob("*.f90"))
errors = []

for path in files:
    lines = path.read_text(encoding="utf-8").splitlines()
    for lineno, line in enumerate(lines, 1):
        low = line.lower()
        if len(line) > 132:
            errors.append(f"{path}:{lineno}: line longer than 132 columns")
        if "double precision" in low or "real*8" in low or re.search(r"\b\d+(?:\.\d*)?[dD][+-]?\d+\b", line):
            errors.append(f"{path}:{lineno}: forbidden real-kind spelling")
        attrs = low.split("::", 1)[0] if "::" in low else low
        if "intent(" in attrs or re.search(r"\bvalue\b", attrs):
            if "::" in line and "!!" not in line:
                errors.append(f"{path}:{lineno}: dummy declaration lacks trailing FORD !! comment")
            if "::" in line and "!!" in line:
                comment = line.split("!!", 1)[1].strip()
                normalized = re.sub(r"[^a-z0-9]+", " ", comment.lower()).strip()
                placeholders = {"", "input", "output", "input value", "output value", "value", "argument"}
                if normalized in placeholders or len(normalized.split()) < 2:
                    errors.append(f"{path}:{lineno}: dummy FORD comment is not meaningful enough")

    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        if re.match(r"(?i)^(pure\s+|elemental\s+|recursive\s+|impure\s+)*(subroutine|function)\s+\w+\s*\(", stripped):
            header = stripped
            start = i
            while ")" not in header.split("!", 1)[0] and i + 1 < len(lines):
                i += 1
                header += " " + lines[i].strip().lstrip("&").strip()
            code = header.split("!", 1)[0]
            match = re.search(r"(?i)\b(?:subroutine|function)\s+\w+\s*\((.*?)\)", code)
            if match:
                args = [a.strip().lstrip("&").strip() for a in match.group(1).split(",") if a.strip()]
                body = lines[i + 1:]
                for arg in args:
                    declarations = []
                    pattern = re.compile(rf"(?i)::\s*{re.escape(arg)}(?:\s*\(|\s*!|\s*$)")
                    for offset, candidate in enumerate(body, i + 2):
                        state = candidate.strip().lower()
                        if state == "contains" or re.match(r"^end\s+(subroutine|function)\b", state):
                            break
                        if pattern.search(candidate):
                            declarations.append((offset, candidate))
                    if not declarations:
                        errors.append(f"{path}:{start+1}: dummy '{arg}' has no separate declaration")
                    else:
                        offset, declaration = declarations[0]
                        dlow = declaration.lower()
                        if "intent(" not in dlow and not re.search(r"\bvalue\b", dlow):
                            errors.append(f"{path}:{offset}: dummy '{arg}' lacks INTENT or VALUE")
                        if "!!" not in declaration:
                            errors.append(f"{path}:{offset}: dummy '{arg}' lacks trailing FORD comment")
            i = max(i, start)
        i += 1

real64_occurrences = []
dp_defs = []
for path in files:
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        low = line.lower()
        if "real64" in low:
            real64_occurrences.append((path, lineno, line))
        if re.search(r"(?i)integer\s*,[^:]*\bparameter\b[^:]*::\s*dp\s*=", line):
            dp_defs.append((path, lineno))
if any(path.name != "classint_kinds.f90" for path, _, _ in real64_occurrences):
    for path, lineno, _ in real64_occurrences:
        if path.name != "classint_kinds.f90":
            errors.append(f"{path}:{lineno}: real64 may only appear in classint_kinds.f90")
if len(dp_defs) != 1 or dp_defs[0][0].name != "classint_kinds.f90":
    errors.append(f"expected exactly one dp definition in classint_kinds.f90, found {dp_defs}")

if errors:
    print("source-rule audit FAILED")
    for error in errors:
        print(error)
    sys.exit(1)
print(f"source-rule audit passed for {len(files)} maintained Fortran files")
