#!/usr/bin/env python3
from pathlib import Path
import re
import sys

roots = [Path('src'), Path('test'), Path('example')]
files = sorted(p for root in roots for p in root.glob('*.f90'))
errors = []
real64_hits = []
for p in files:
    text = p.read_text(encoding='utf-8')
    low = text.lower()
    if 'real64' in low:
        real64_hits.append(str(p))
    for pattern, name in [
        (r'\bdouble\s+precision\b', 'double precision'),
        (r'\breal\s*\*\s*8\b', 'real*8'),
        (r'kind\s*\(\s*0\.0d0\s*\)', 'kind(0.0d0)'),
        (r'(?<![a-z0-9_])(?:\d+(?:\.\d*)?|\.\d+)[d][+-]?\d+', 'D exponent'),
    ]:
        if re.search(pattern, low):
            errors.append(f'{p}: forbidden {name}')
    for lineno, line in enumerate(text.splitlines(), 1):
        if len(line) > 132:
            errors.append(f'{p}:{lineno}: line length {len(line)} > 132')

if real64_hits != ['src/rpart_kinds.f90']:
    errors.append(f'real64 must occur only in src/rpart_kinds.f90; got {real64_hits}')

if errors:
    print('\n'.join(errors))
    sys.exit(1)
print(f'check_source: PASS ({len(files)} maintained Fortran files)')
