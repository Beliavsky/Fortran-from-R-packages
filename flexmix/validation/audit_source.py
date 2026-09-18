#!/usr/bin/env python3
"""Static release-policy audit for the flexmix Fortran translation."""
from __future__ import annotations
from pathlib import Path
import hashlib
import re
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[1]
FORTRAN_DIRS = [ROOT / "src", ROOT / "test", ROOT / "example", ROOT / "validation"]
SOURCES = sorted(p for d in FORTRAN_DIRS for p in d.glob("*.f90"))
errors: list[str] = []


def code_before_comment(line: str) -> str:
    out = []
    quote = None
    i = 0
    while i < len(line):
        c = line[i]
        if quote:
            out.append(c)
            if c == quote:
                if i + 1 < len(line) and line[i + 1] == quote:
                    out.append(line[i + 1])
                    i += 2
                    continue
                quote = None
        else:
            if c in "'\"":
                quote = c
                out.append(c)
            elif c == "!":
                break
            else:
                out.append(c)
        i += 1
    return "".join(out)


def split_top_level(text: str) -> list[str]:
    parts, start, depth = [], 0, 0
    quote = None
    for i, c in enumerate(text):
        if quote:
            if c == quote:
                quote = None
            continue
        if c in "'\"":
            quote = c
        elif c == "(":
            depth += 1
        elif c == ")":
            depth = max(0, depth - 1)
        elif c == "," and depth == 0:
            parts.append(text[start:i].strip())
            start = i + 1
    parts.append(text[start:].strip())
    return [p for p in parts if p]


def logical_statements(lines: list[str]) -> list[str]:
    statements: list[str] = []
    current = ""
    continuing = False
    for raw in lines:
        code = code_before_comment(raw).strip()
        if not code:
            continue
        if continuing and code.startswith("&"):
            code = code[1:].lstrip()
        has_cont = code.endswith("&")
        if has_cont:
            code = code[:-1].rstrip()
        if current:
            current += " " + code
        else:
            current = code
        continuing = has_cont
        if not continuing:
            statements.append(current)
            current = ""
    if current:
        statements.append(current)
    return statements

# Disallowed syntax and source hygiene.
for path in SOURCES:
    text = path.read_text(encoding="utf-8")
    rel = path.relative_to(ROOT)
    for lineno, line in enumerate(text.splitlines(), 1):
        code = code_before_comment(line)
        if ";" in code:
            errors.append(f"{rel}:{lineno}: semicolon-separated statement")
        if re.search(r"\bdouble\s+precision\b|\breal\s*\*\s*8\b", code, re.I):
            errors.append(f"{rel}:{lineno}: legacy real-kind declaration")
        if re.search(r"(?:\d(?:\.\d*)?|\.\d+)[dD][+-]?\d+", code):
            errors.append(f"{rel}:{lineno}: D-exponent literal")
        if re.search(r"\bkind\s*\(\s*[-+]?\d*\.\d+[dD]0\s*\)", code, re.I):
            errors.append(f"{rel}:{lineno}: kind(0.0d0)-style declaration")
        m = re.search(r"\b([A-Za-z_]\w*)\s*/=\s*\1\b", code)
        if m:
            errors.append(f"{rel}:{lineno}: self-comparison NaN test")
        if re.search(r"-ffast-math|-Ofast|-ffinite-math-only", code):
            errors.append(f"{rel}:{lineno}: unsafe finite-math compiler flag")

    # Every dummy must have an INTENT/VALUE declaration; dummy declarations must be one per line with FORD docs.
    statements = logical_statements(text.splitlines())
    procedures: dict[str, list[str]] = {}
    for stmt in statements:
        m = re.search(r"\b(?:subroutine|function)\s+([A-Za-z_]\w*)\s*\(([^)]*)\)", stmt, re.I)
        if m:
            dummies = [x.strip() for x in m.group(2).split(",") if x.strip()]
            procedures[m.group(1).lower()] = [x.lower() for x in dummies]
    declarations: dict[str, tuple[int, str, bool, int]] = {}
    for lineno, line in enumerate(text.splitlines(), 1):
        code = code_before_comment(line)
        if "::" not in code:
            continue
        lhs, rhs = code.split("::", 1)
        if not (re.search(r"\bintent\s*\(", lhs, re.I) or re.search(r"\bvalue\b", lhs, re.I)):
            continue
        names = split_top_level(rhs)
        clean_names = []
        for name in names:
            name = re.sub(r"\s*=.*$", "", name).strip()
            name = re.sub(r"\(.*\)$", "", name).strip()
            if re.match(r"^[A-Za-z_]\w*$", name):
                clean_names.append(name.lower())
        has_ford = "!!" in line
        for name in clean_names:
            declarations[name] = (lineno, line, has_ford, len(clean_names))
    for proc, dummies in procedures.items():
        for dummy in dummies:
            if dummy == "*":
                continue
            dec = declarations.get(dummy)
            if dec is None:
                errors.append(f"{rel}: procedure {proc}: dummy {dummy} lacks explicit INTENT or VALUE")
                continue
            lineno, line, has_ford, ndecl = dec
            if ndecl != 1:
                errors.append(f"{rel}:{lineno}: dummy {dummy} is not declared on its own line")
            if not has_ford:
                errors.append(f"{rel}:{lineno}: dummy {dummy} lacks trailing FORD '!!' documentation")

# One canonical real kind, and every maintained source that declares reals imports/uses dp.
kinds = (ROOT / "src" / "flexmix_kinds.f90").read_text()
if not re.search(r"integer\s*,\s*parameter\s*,\s*public\s*::\s*dp\s*=\s*real64", kinds, re.I):
    errors.append("src/flexmix_kinds.f90: canonical dp=real64 definition missing")
for path in sorted((ROOT / "src").glob("*.f90")):
    text = path.read_text()
    if path.name == "flexmix_kinds.f90":
        continue
    if re.search(r"\breal\s*\(", text, re.I) and "real(dp)" not in text.lower():
        errors.append(f"{path.relative_to(ROOT)}: real declarations do not use dp")

# Duplicate Fortran content.
seen: dict[str, Path] = {}
for path in SOURCES:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest in seen:
        errors.append(f"duplicate Fortran source content: {path.relative_to(ROOT)} and {seen[digest].relative_to(ROOT)}")
    seen[digest] = path

# No build products, archives, caches, or vendored dependency directories.
bad_suffixes = {".o", ".obj", ".mod", ".smod", ".exe", ".a", ".so", ".dll", ".zip", ".pyc"}
for path in ROOT.rglob("*"):
    if path.is_file() and path.suffix.lower() in bad_suffixes:
        errors.append(f"build/archive artifact present: {path.relative_to(ROOT)}")
    if path.is_dir() and path.name == "__pycache__":
        errors.append(f"cache directory present: {path.relative_to(ROOT)}")
for dep in ["rfortran-core", "rfortran-linalg", "rfortran-arpack", "rfortran-compat", "nnet", "mclust", "fortran-lapack"]:
    if (ROOT / dep).exists():
        errors.append(f"vendored dependency directory present: {dep}")

# Manifest/README/source consistency.
data = tomllib.loads((ROOT / "fpm.toml").read_text())
translation = data["extra"]["translation"]
mapped = translation["function"]
if translation["coverage_basis"] != "exported computational R functions":
    errors.append("fpm.toml: unexpected coverage_basis")
if translation["r_functions_total"] != 68 or translation["r_functions_translated"] != 66:
    errors.append("fpm.toml: coverage counts are not 66 of 68")
if len(mapped) != 66:
    errors.append(f"fpm.toml: {len(mapped)} function mappings, expected 66")
expected_fraction = 66.0 / 68.0
if abs(translation["frac_functions_translated"] - expected_fraction) > 1e-15:
    errors.append("fpm.toml: coverage fraction does not equal 66/68")
mapped_names = [m["r_name"] for m in mapped]
untranslated = translation["untranslated_r_functions"]
if len(untranslated) != 2 or set(mapped_names) & set(untranslated):
    errors.append("fpm.toml: mapped/untranslated partition is inconsistent")
if len(set(mapped_names)) != len(mapped_names):
    errors.append("fpm.toml: duplicate r_name mapping")
for m in mapped:
    rsource = ROOT / m["r_source"]
    if not rsource.is_file():
        errors.append(f"fpm.toml: missing R source {m['r_source']}")
    source_names = m.get("fortran_sources")
    if source_names is None:
        source_names = [m.get("fortran_source", "")]
    ftexts = []
    for source_name in source_names:
        fsource = ROOT / source_name
        if not source_name or not fsource.is_file():
            errors.append(f"fpm.toml: missing Fortran source {source_name}")
        else:
            ftexts.append(fsource.read_text())
    combined = "\n".join(ftexts)
    for name in m["fortran_names"]:
        if not re.search(rf"\b(?:subroutine|function)\s+{re.escape(name)}\b", combined, re.I):
            errors.append(f"fpm.toml: {name} not found in mapped Fortran source(s) for {m['r_name']}")
readme = (ROOT / "README.md").read_text()
if "66 of 68 (97.1%)" not in readme:
    errors.append("README.md: coverage headline disagrees with manifest")

# Recompute the 68-name coverage basis from the preserved NAMESPACE.
ns = (ROOT / "upstream" / "NAMESPACE").read_text()
export_block = re.search(r"export\((.*?)\)\s*\n\n\nexportClasses", ns, re.S)
method_block = re.search(r"exportMethods\((.*?)\)\s*\n\nS3method", ns, re.S)
if export_block and method_block:
    exports = re.findall(r'"([^"]+)"', export_block.group(1))
    methods = re.findall(r'"([^"]+)"', method_block.group(1))
    excluded = {"flxglht", "plotEll", "coerce", "initialize", "plot", "show", "summary"}
    basis = []
    for name in exports + methods:
        if name not in excluded and name not in basis:
            basis.append(name)
    if len(basis) != 68:
        errors.append(f"coverage basis recomputation yielded {len(basis)}, expected 68")
    if set(basis) != set(mapped_names) | set(untranslated):
        errors.append("manifest mapped+untranslated names do not equal audited NAMESPACE coverage basis")
else:
    errors.append("could not parse preserved upstream NAMESPACE")

if errors:
    print("SOURCE AUDIT FAILED")
    for error in errors:
        print("-", error)
    sys.exit(1)
print(f"SOURCE AUDIT PASSED: {len(SOURCES)} Fortran files; 66/68 coverage mapping consistent.")
