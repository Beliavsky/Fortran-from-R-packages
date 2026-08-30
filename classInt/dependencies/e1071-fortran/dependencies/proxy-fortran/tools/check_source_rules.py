#!/usr/bin/env python3
"""Audit maintained Fortran sources against the repository translation rules."""
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
                if len(comment.split()) < 4:
                    errors.append(f"{path}:{lineno}: dummy FORD comment is not meaningful enough")

    # Build logical procedure headers, including continuation lines.
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
            m = re.search(r"(?i)\b(?:subroutine|function)\s+\w+\s*\((.*?)\)", code)
            if m:
                args = [a.strip().lstrip("&").strip() for a in m.group(1).split(",") if a.strip()]
                # Search until the first executable-looking line or nested contains/end.
                body = lines[i + 1:]
                for arg in args:
                    decls = []
                    pat = re.compile(rf"(?i)::\s*{re.escape(arg)}(?:\s*\(|\s*!|\s*$)")
                    for off, candidate in enumerate(body, i + 2):
                        st = candidate.strip().lower()
                        if st == "contains" or re.match(r"^end\s+(subroutine|function)\b", st):
                            break
                        if pat.search(candidate):
                            decls.append((off, candidate))
                    if not decls:
                        errors.append(f"{path}:{start+1}: dummy '{arg}' has no separate declaration")
                    else:
                        off, decl = decls[0]
                        dl = decl.lower()
                        if "intent(" not in dl and not re.search(r"\bvalue\b", dl):
                            errors.append(f"{path}:{off}: dummy '{arg}' lacks INTENT or VALUE")
                        if "!!" not in decl:
                            errors.append(f"{path}:{off}: dummy '{arg}' lacks trailing FORD comment")
            i = max(i, start)
        i += 1

# dp/real64 rules.
real64_occurrences = []
dp_defs = []
for path in files:
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        low = line.lower()
        if "real64" in low:
            real64_occurrences.append((path, lineno, line))
        if re.search(r"(?i)integer\s*,\s*parameter[^:]*::\s*dp\s*=", line):
            dp_defs.append((path, lineno))
if any(path.name != "proxy_kinds.f90" for path, _, _ in real64_occurrences):
    for path, lineno, _ in real64_occurrences:
        if path.name != "proxy_kinds.f90":
            errors.append(f"{path}:{lineno}: real64 may only appear in proxy_kinds.f90")
if len(dp_defs) != 1 or dp_defs[0][0].name != "proxy_kinds.f90":
    errors.append(f"expected exactly one dp definition in proxy_kinds.f90, found {dp_defs}")

if errors:
    print("source-rule audit FAILED")
    for error in errors:
        print(error)
    sys.exit(1)
print(f"source-rule audit passed for {len(files)} maintained Fortran files")
