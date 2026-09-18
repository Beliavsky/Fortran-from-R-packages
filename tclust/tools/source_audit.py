#!/usr/bin/env python3
from pathlib import Path
import hashlib
import re
import tomllib

root = Path(__file__).resolve().parents[1]
fortran = []
for sub in ('src', 'test', 'example', 'tools'):
    fortran.extend(sorted((root / sub).glob('*.f90')))
errors = []
seen = {}
legacy = re.compile(r'\bdouble\s+precision\b|\breal\s*\*\s*8\b|\bkind\s*\(\s*0\.0d0\s*\)|(?<![A-Za-z0-9_])\d+(?:\.\d*)?[dD][+-]?\d+', re.I)
selfcmp = re.compile(r'\b([A-Za-z_]\w*)\s*/=\s*\1\b|\b([A-Za-z_]\w*)\s*==\s*\2\b', re.I)

def code_part(line):
    in_quote = None
    out = []
    i = 0
    while i < len(line):
        c = line[i]
        if in_quote:
            out.append(c)
            if c == in_quote:
                if i + 1 < len(line) and line[i+1] == in_quote:
                    out.append(line[i+1]); i += 1
                else:
                    in_quote = None
        else:
            if c in "'\"":
                in_quote = c; out.append(c)
            elif c == '!':
                break
            else:
                out.append(c)
        i += 1
    return ''.join(out)

def has_semicolon_outside_string(text):
    quote = None
    i = 0
    while i < len(text):
        c = text[i]
        if quote:
            if c == quote:
                if i + 1 < len(text) and text[i+1] == quote:
                    i += 1
                else:
                    quote = None
        elif c in "'\"":
            quote = c
        elif c == ';':
            return True
        i += 1
    return False

def top_level_commas(text):
    depth = 0
    quote = None
    count = 0
    for c in text:
        if quote:
            if c == quote:
                quote = None
        elif c in "'\"": quote = c
        elif c == '(': depth += 1
        elif c == ')': depth = max(0, depth-1)
        elif c == ',' and depth == 0: count += 1
    return count

for path in fortran:
    text = path.read_text(encoding='utf-8')
    digest = hashlib.sha256(text.encode()).hexdigest()
    if digest in seen:
        errors.append(f'duplicate Fortran content: {path.relative_to(root)} and {seen[digest]}')
    seen[digest] = path.relative_to(root)
    for lineno, line in enumerate(text.splitlines(), 1):
        if len(line) > 132:
            errors.append(f'{path.relative_to(root)}:{lineno}: line length {len(line)} > 132')
        code = code_part(line)
        if has_semicolon_outside_string(code):
            errors.append(f'{path.relative_to(root)}:{lineno}: semicolon in executable/source code')
        if legacy.search(code):
            errors.append(f'{path.relative_to(root)}:{lineno}: legacy real kind or D exponent')
        if selfcmp.search(code):
            errors.append(f'{path.relative_to(root)}:{lineno}: self-comparison NaN idiom')
        low = code.lower()
        if 'intent(' in low or re.search(r'\bvalue\b', low):
            if '::' in code:
                decl = code.split('::', 1)[1].strip()
                if top_level_commas(decl) > 0:
                    errors.append(f'{path.relative_to(root)}:{lineno}: multiple dummy arguments on one declaration')
                if '!!' not in line:
                    errors.append(f'{path.relative_to(root)}:{lineno}: dummy declaration lacks trailing FORD !! comment')

src_text = '\n'.join(p.read_text(encoding='utf-8') for p in root.joinpath('src').glob('*.f90'))
if re.search(r'\buse\s*,?\s*intrinsic\s*::\s*iso_fortran_env', src_text, re.I):
    occurrences = len(re.findall(r'\biso_fortran_env\b', src_text, re.I))
    if occurrences != 1:
        errors.append(f'iso_fortran_env should define dp once; found {occurrences} references')
if 'real64' not in (root/'src'/'tclust_kinds.f90').read_text().lower():
    errors.append('tclust_kinds.f90 does not define dp from real64')

with open(root/'fpm.toml', 'rb') as f:
    manifest = tomllib.load(f)
tr = manifest['extra']['translation']
entries = tr.get('function', [])
total = tr['r_functions_total']
translated = tr['r_functions_translated']
frac = tr['frac_functions_translated']
untranslated = tr['untranslated_r_functions']
if total != 11 or translated != 11 or len(entries) != 11 or untranslated != []:
    errors.append('fpm.toml translation counts/lists do not equal audited 11/11 coverage')
if abs(frac - translated/total) > 1e-15:
    errors.append('fpm.toml translation fraction inconsistent with counts')
for entry in entries:
    names = entry.get('fortran_names')
    if not isinstance(names, list) or not names:
        errors.append(f"{entry.get('r_name')}: fortran_names must be a nonempty array")

readme = (root/'README.md').read_text(encoding='utf-8')
if '11 of 11 (100%)' not in readme or 'Package status: substantial' not in readme:
    errors.append('README coverage summary does not match manifest')
api = (root/'API_COVERAGE.md').read_text(encoding='utf-8')
if '11 of 11' not in api or '100%' not in api:
    errors.append('API_COVERAGE.md coverage summary does not match manifest')

bad_artifact = re.compile(r'\.(o|obj|mod|smod|exe|a|so|dll|pyc)$', re.I)
for p in root.rglob('*'):
    if p.is_file() and (bad_artifact.search(p.name) or p.name == '__pycache__'):
        errors.append(f'build artifact in package tree: {p.relative_to(root)}')
for dirname in ('build', '.fpm', '__pycache__'):
    if any(p.is_dir() and p.name == dirname for p in root.rglob('*')):
        errors.append(f'cache/build directory in package tree: {dirname}')

if errors:
    print('SOURCE AUDIT FAILED')
    for e in errors:
        print(' -', e)
    raise SystemExit(1)
print(f'SOURCE AUDIT PASSED: {len(fortran)} Fortran files; 11/11 coverage mapping consistent.')
