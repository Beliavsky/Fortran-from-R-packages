#!/usr/bin/env python3
"""Static release-policy checks for the pbs Fortran translation."""
from __future__ import annotations

import hashlib
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FORTRAN_DIRS = [ROOT / 'src', ROOT / 'test', ROOT / 'example', ROOT / 'tools']
fortran = sorted(p for d in FORTRAN_DIRS for p in d.glob('*.f90'))
errors: list[str] = []

if not fortran:
    errors.append('no Fortran source files found')

# Legacy kinds, disallowed D exponents, self-comparison NaN tests, and line length.
legacy_patterns = [
    (re.compile(r'\bdouble\s+precision\b', re.I), 'double precision'),
    (re.compile(r'\breal\s*\*\s*8\b', re.I), 'real*8'),
    (re.compile(r'\bkind\s*\(\s*0?\.0d0\s*\)', re.I), 'kind(0.0d0)'),
    (re.compile(r'(?<![A-Za-z0-9_])(?:\d+(?:\.\d*)?|\.\d+)[dD][+-]?\d+'), 'D exponent'),
]
for p in fortran:
    text = p.read_text()
    for lineno, line in enumerate(text.splitlines(), 1):
        if len(line) > 132:
            errors.append(f'{p.relative_to(ROOT)}:{lineno}: line exceeds 132 columns')
        code = line.split('!')[0]
        if ';' in code:
            errors.append(f'{p.relative_to(ROOT)}:{lineno}: executable semicolon found')
        if re.search(r'\b([A-Za-z][A-Za-z0-9_]*)\s*/=\s*\1\b', code):
            errors.append(f'{p.relative_to(ROOT)}:{lineno}: self-comparison NaN idiom')
        for pat, label in legacy_patterns:
            if pat.search(code):
                errors.append(f'{p.relative_to(ROOT)}:{lineno}: disallowed {label}')

# Exactly one real64-backed dp definition in maintained source.
dp_defs = []
for p in fortran:
    for lineno, line in enumerate(p.read_text().splitlines(), 1):
        if re.search(r'\binteger\s*,\s*parameter.*\bdp\s*=\s*real64\b', line, re.I):
            dp_defs.append((p, lineno))
if len(dp_defs) != 1:
    errors.append(f'expected exactly one dp=real64 definition, found {len(dp_defs)}')

# No duplicate maintained Fortran file contents.
seen: dict[str, Path] = {}
for p in fortran:
    digest = hashlib.sha256(p.read_bytes()).hexdigest()
    if digest in seen:
        errors.append(f'duplicate Fortran contents: {seen[digest].relative_to(ROOT)} and {p.relative_to(ROOT)}')
    seen[digest] = p

# Dummy arguments: explicit intent/value and trailing FORD comment.
proc_re = re.compile(r'^\s*(?:pure\s+)?(?:elemental\s+)?(?:function|subroutine)\s+([A-Za-z0-9_]+)\s*\(([^)]*)\)', re.I)
for p in fortran:
    lines = p.read_text().splitlines()
    for idx, line in enumerate(lines):
        m = proc_re.match(line)
        if not m:
            continue
        args = [a.strip() for a in m.group(2).split(',') if a.strip()]
        if not args:
            continue
        block = '\n'.join(lines[idx + 1 : min(len(lines), idx + 30)])
        for arg in args:
            decl_line = None
            for candidate in block.splitlines():
                if '::' not in candidate:
                    continue
                if not re.search(rf'::\s*{re.escape(arg)}(?:\s*\([^)]*\))?\s*(?:!!|$)', candidate, re.I):
                    continue
                if not re.search(r'\bintent\s*\([^)]*\)|\bvalue\b', candidate, re.I):
                    continue
                decl_line = candidate
                break
            if decl_line is None:
                errors.append(f'{p.relative_to(ROOT)}:{idx+1}: dummy {arg} lacks explicit INTENT/VALUE declaration')
                continue
            if '!!' not in decl_line:
                errors.append(f'{p.relative_to(ROOT)}:{idx+1}: dummy {arg} lacks trailing FORD comment')

# Coverage manifest consistency.
manifest = tomllib.loads((ROOT / 'fpm.toml').read_text())
translation = manifest['extra']['translation']
mappings = translation['function']
if translation['r_functions_total'] != 2:
    errors.append('manifest r_functions_total is not 2')
if translation['r_functions_translated'] != len(mappings):
    errors.append('manifest translated count does not equal mapping entries')
if len(mappings) != 2:
    errors.append('expected two function mapping entries')
if translation['untranslated_r_functions']:
    errors.append('untranslated_r_functions should be empty')
if abs(float(translation['frac_functions_translated']) - 1.0) > 1e-15:
    errors.append('coverage fraction is not 1.0')
readme = (ROOT / 'README.md').read_text()
if '2 of 2 (100%)' not in readme:
    errors.append('README coverage statement does not match manifest')

# No generated build products in the package tree.
for p in ROOT.rglob('*'):
    if p.is_file() and (p.suffix.lower() in {'.o', '.mod', '.smod', '.exe', '.zip'} or p.name in {'a.out'}):
        errors.append(f'build/archive product present: {p.relative_to(ROOT)}')
for name in ['build', '.fpm', '__pycache__', '.pytest_cache']:
    for p in ROOT.rglob(name):
        if p.is_dir():
            errors.append(f'cache/build directory present: {p.relative_to(ROOT)}')

if errors:
    print('SOURCE AUDIT FAILED')
    for e in errors:
        print('-', e)
    sys.exit(1)
print(f'SOURCE AUDIT PASSED: {len(fortran)} Fortran files; 2/2 coverage mapping consistent.')
