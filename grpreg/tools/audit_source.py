#!/usr/bin/env python3
"""Static release-policy audit for the grpreg Fortran translation."""
from __future__ import annotations
from pathlib import Path
import hashlib
import re
import tomllib
import sys

root = Path(__file__).resolve().parents[1]
fortran = sorted(list((root/'src').glob('*.f90')) + list((root/'test').glob('*.f90')) +
                 list((root/'example').glob('*.f90')) + list((root/'tools').glob('*.f90')))
errors: list[str] = []

legacy = re.compile(r"double\s+precision|real\s*\*\s*8|kind\s*\(\s*0\.0d0|\b\d+(?:\.\d*)?[dD][+-]?\d+", re.I)
selfcmp = re.compile(r"\b([A-Za-z_]\w*(?:%\w+)*)\s*/=\s*\1\b", re.I)
for p in fortran:
    lines = p.read_text(encoding='utf-8').splitlines()
    for no,line in enumerate(lines,1):
        if len(line) > 132:
            errors.append(f"{p.relative_to(root)}:{no}: line length {len(line)} > 132")
        code = line.split('!')[0]
        if ';' in code:
            errors.append(f"{p.relative_to(root)}:{no}: executable semicolon")
        if legacy.search(code):
            errors.append(f"{p.relative_to(root)}:{no}: legacy real kind or D exponent")
        if selfcmp.search(code):
            errors.append(f"{p.relative_to(root)}:{no}: self-comparison NaN idiom")
        if 'intent(' in line.lower() and '!!' not in line:
            errors.append(f"{p.relative_to(root)}:{no}: dummy declaration lacks trailing FORD comment")

real64_files = [p.relative_to(root).as_posix() for p in fortran if re.search(r'\breal64\b',p.read_text(),re.I)]
if real64_files != ['src/grpreg_kinds.f90']:
    errors.append(f"real64 must appear only in src/grpreg_kinds.f90, found {real64_files}")

hashes: dict[str,list[str]] = {}
for p in (root/'src').glob('*.f90'):
    h = hashlib.sha256(p.read_bytes()).hexdigest()
    hashes.setdefault(h,[]).append(p.name)
for names in hashes.values():
    if len(names) > 1:
        errors.append(f"duplicate Fortran source contents: {names}")

bad_suffix = {'.o','.mod','.smod','.exe','.obj','.pyc','.zip'}
for p in root.rglob('*'):
    if p.is_file() and (p.suffix.lower() in bad_suffix or p.name in {'a.out'}):
        errors.append(f"build/archive artifact in package tree: {p.relative_to(root)}")
    if '__pycache__' in p.parts:
        errors.append(f"cache directory in package tree: {p.relative_to(root)}")

manifest = tomllib.loads((root/'fpm.toml').read_text())
tr = manifest.get('extra',{}).get('translation',{})
maps = tr.get('function',[])
total = tr.get('r_functions_total')
translated = tr.get('r_functions_translated')
untranslated = tr.get('untranslated_r_functions',[])
frac = tr.get('frac_functions_translated')
if total != 19 or translated != 19 or len(maps) != 19 or untranslated != [] or abs(float(frac)-1.0) > 1e-15:
    errors.append(f"coverage manifest inconsistent: total={total}, translated={translated}, mappings={len(maps)}, untranslated={untranslated}, frac={frac}")
readme=(root/'README.md').read_text()
if '19 of 19 (100%)' not in readme:
    errors.append('README coverage does not state 19 of 19 (100%)')

# Forbid copied external dependencies in the maintained source tree.
for name in ['r.f90','r_mod.f90','blas.f90','lapack.f90','arpack.f90']:
    if any(p.name.lower()==name for p in root.rglob('*') if p.is_file() and 'upstream' not in p.parts):
        errors.append(f"forbidden copied dependency source: {name}")

if errors:
    print('\n'.join('ERROR: '+e for e in errors))
    sys.exit(1)
print(f"SOURCE AUDIT PASSED: {len(fortran)} Fortran files; 19/19 coverage mapping consistent.")
