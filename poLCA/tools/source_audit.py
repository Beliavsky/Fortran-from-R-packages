from __future__ import annotations

import hashlib
import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FORTRAN_DIRS = [ROOT / 'src', ROOT / 'test', ROOT / 'example', ROOT / 'tools']
files = sorted(p for d in FORTRAN_DIRS for p in d.glob('*.f90'))
errors: list[str] = []

# Source-level policy checks.
for path in files:
    text = path.read_text()
    rel = path.relative_to(ROOT)
    for lineno, line in enumerate(text.splitlines(), 1):
        if len(line) > 132:
            errors.append(f'{rel}:{lineno}: line length {len(line)} exceeds 132')
        code = line.split('!')[0]
        if ';' in code and not re.search(r'\bcase\s*\([^)]*\)\s*;', code, re.I):
            errors.append(f'{rel}:{lineno}: semicolon-separated statement')
        if re.search(r'\bdouble\s+precision\b|\breal\s*\*\s*8\b', code, re.I):
            errors.append(f'{rel}:{lineno}: legacy real declaration')
        if re.search(r'(?<![A-Za-z0-9_])(?:\d+(?:\.\d*)?|\.\d+)[dD][+-]?\d+', code):
            errors.append(f'{rel}:{lineno}: D-exponent literal')
        m = re.search(r'\b([A-Za-z_]\w*)\s*/=\s*\1\b', code, re.I)
        if m:
            errors.append(f'{rel}:{lineno}: self-comparison NaN test')
        if re.search(r'\b-?f(?:fast-math|finite-math-only)\b|\b-Ofast\b', line, re.I):
            errors.append(f'{rel}:{lineno}: disallowed floating-point compiler option')
        if 'intent(' in code.lower() or re.search(r'\bvalue\b', code, re.I):
            if '::' in code:
                rhs = code.split('::', 1)[1].strip()
                # Dummy declarations must contain exactly one entity before any initialization.
                depth = 0
                commas = 0
                for ch in rhs:
                    if ch == '(':
                        depth += 1
                    elif ch == ')':
                        depth -= 1
                    elif ch == ',' and depth == 0:
                        commas += 1
                if commas:
                    errors.append(f'{rel}:{lineno}: multiple dummy arguments declared on one line')
                if '!!' not in line:
                    errors.append(f'{rel}:{lineno}: dummy declaration lacks trailing FORD !! comment')

    # Collect procedure signatures and verify each dummy declaration explicitly.
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        raw = lines[i]
        code = raw.split('!')[0]
        if re.search(r'\b(subroutine|function)\b', code, re.I) and not re.search(r'^\s*end\s+', code, re.I):
            sig = code.strip()
            j = i
            while ')' not in sig and j + 1 < len(lines):
                j += 1
                sig += ' ' + lines[j].split('!')[0].strip()
            sig = sig.replace('&', ' ')
            mproc = re.search(r'\b(?:subroutine|function)\s+([A-Za-z_]\w*)\s*\(([^)]*)\)', sig, re.I)
            if mproc:
                name = mproc.group(1)
                args = [a.strip() for a in mproc.group(2).split(',') if a.strip()]
                k = j + 1
                block = []
                while k < len(lines):
                    if re.search(rf'^\s*end\s+(?:subroutine|function)\s+{re.escape(name)}\b', lines[k], re.I):
                        break
                    block.append(lines[k])
                    k += 1
                for arg in args:
                    found = False
                    for decl_line in block:
                        before_comment = decl_line.split('!')[0]
                        if '::' not in before_comment:
                            continue
                        rhs = before_comment.split('::', 1)[1]
                        if re.search(rf'\b{re.escape(arg)}\b', rhs, re.I):
                            attrs = before_comment.split('::', 1)[0]
                            if re.search(r'\bintent\s*\([^)]*\)|\bvalue\b', attrs, re.I) and '!!' in decl_line:
                                found = True
                                break
                    if not found:
                        errors.append(
                            f'{rel}: procedure {name}: dummy {arg} lacks explicit INTENT/VALUE and FORD comment'
                        )
                i = max(i, k)
        i += 1

# Duplicate maintained Fortran content.
by_hash: dict[str, list[Path]] = {}
for path in files:
    h = hashlib.sha256(path.read_bytes()).hexdigest()
    by_hash.setdefault(h, []).append(path)
for group in by_hash.values():
    if len(group) > 1:
        errors.append('duplicate Fortran source content: ' + ', '.join(str(p.relative_to(ROOT)) for p in group))

# Coverage consistency.
manifest = tomllib.loads((ROOT / 'fpm.toml').read_text())
tr = manifest['extra']['translation']
entries = tr['function']
if tr['r_functions_total'] != 10:
    errors.append('fpm.toml r_functions_total is not 10')
if tr['r_functions_translated'] != len(entries):
    errors.append('fpm.toml mapping count does not equal r_functions_translated')
frac = tr['r_functions_translated'] / tr['r_functions_total']
if abs(tr['frac_functions_translated'] - frac) > 1e-12:
    errors.append('fpm.toml translation fraction is inconsistent')
if tr['untranslated_r_functions']:
    errors.append('untranslated_r_functions should be empty at 10/10 coverage')
expected = {
    'poLCA', 'poLCA.simdata', 'poLCA.reorder', 'poLCA.table', 'poLCA.entropy',
    'poLCA.predcell', 'poLCA.posterior', 'rmulti', 'coef.poLCA', 'vcov.poLCA'
}
seen = {e['r_name'] for e in entries}
if seen != expected:
    errors.append(f'coverage entry names differ: seen={sorted(seen)}')
readme = (ROOT / 'README.md').read_text()
if '10 of 10 (100%)' not in readme:
    errors.append('README coverage summary does not say 10 of 10 (100%)')

if errors:
    print('SOURCE AUDIT FAILED')
    for e in errors:
        print(' -', e)
    raise SystemExit(1)

print(f'SOURCE AUDIT PASSED: {len(files)} Fortran files; 10/10 coverage mapping consistent.')
