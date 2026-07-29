#!/usr/bin/env python3
from __future__ import annotations
import hashlib
import sys
import tomllib
from pathlib import Path

root = Path(__file__).resolve().parents[1]
errors: list[str] = []

with (root / 'fpm.toml').open('rb') as handle:
    manifest = tomllib.load(handle)
if manifest.get('license') != 'GPL-3.0-or-later':
    errors.append('unexpected license in fpm.toml')

for path in sorted(root.rglob('*.f90')):
    text = path.read_text(encoding='ascii')
    lines = text.splitlines()
    if not any('SPDX-License-Identifier: GPL-3.0-or-later' in line for line in lines[:8]):
        errors.append(f'missing SPDX header: {path.relative_to(root)}')
    if not any('implicit none' in line.lower() for line in lines):
        errors.append(f'missing implicit none: {path.relative_to(root)}')
    for number, line in enumerate(lines, 1):
        if len(line) > 132:
            errors.append(f'line longer than 132 columns: {path.relative_to(root)}:{number}')

for path in root.rglob('*'):
    if path.is_file() and path.suffix.lower() in {'.f90', '.md', '.toml', '.sh', '.bat', '.txt'}:
        try:
            path.read_text(encoding='ascii')
        except UnicodeDecodeError:
            errors.append(f'non-ASCII release text: {path.relative_to(root)}')

for manifest_name, base in [
    ('ORIGINAL_FILES_SHA256.txt', root / 'original' / 'nvmix-0.1-2'),
    ('TRANSLATED_FILES_SHA256.txt', root),
]:
    manifest_path = root / 'provenance' / manifest_name
    if not manifest_path.exists():
        errors.append(f'missing {manifest_name}')
        continue
    for line in manifest_path.read_text(encoding='ascii').splitlines():
        if not line.strip():
            continue
        expected, relative = line.split('  ', 1)
        path = base / relative
        if not path.exists():
            errors.append(f'missing hashed file: {relative}')
            continue
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            errors.append(f'checksum mismatch: {relative}')

if errors:
    print('\n'.join(errors))
    sys.exit(1)
print('archive_audit: PASS')
