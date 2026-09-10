#!/usr/bin/env python3
"""Audit maintained Fortran source rules used by this translation."""
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
FILES = sorted(ROOT.joinpath('src').glob('*.f90')) + sorted(ROOT.joinpath('test').glob('*.f90')) + sorted(ROOT.joinpath('example').glob('*.f90'))
errors = []
proc_re = re.compile(r'(?i)^\s*(?:(?:pure|elemental|recursive)\s+)*(?:integer\s+|real\([^)]*\)\s+|logical\s+|character\([^)]*\)\s+)?(subroutine|function)\s+(\w+)\s*\(')

for path in FILES:
    lines = path.read_text(encoding='utf-8').splitlines()
    for lineno, line in enumerate(lines, 1):
        code = line.split('!')[0]
        if ';' in code:
            errors.append(f'{path}:{lineno}: semicolon in executable/declaration code')
        if re.search(r'(?i)\bdouble\s+precision\b|\breal\s*\*\s*\d+|\b\d+(?:\.\d*)?[dD][+-]?\d+', code):
            errors.append(f'{path}:{lineno}: legacy real-kind syntax')
        if re.search(r'(?i)\b([A-Za-z_]\w*)\s*/=\s*\1\b', code):
            errors.append(f'{path}:{lineno}: self-comparison NaN idiom')

    i = 0
    while i < len(lines):
        match = proc_re.match(lines[i])
        if not match:
            i += 1
            continue
        kind, proc_name = match.group(1).lower(), match.group(2)
        start = i
        signature = lines[i].strip()
        # Continue until the procedure argument list closes, ignoring RESULT after it.
        depth = signature.count('(') - signature.count(')')
        while depth > 0 and i + 1 < len(lines):
            i += 1
            signature += ' ' + lines[i].strip()
            depth += lines[i].count('(') - lines[i].count(')')
        signature = signature.replace('&', ' ')
        arg_match = re.search(r'(?i)(?:subroutine|function)\s+' + re.escape(proc_name) + r'\s*\((.*?)\)', signature)
        args = [] if arg_match is None else [a.strip() for a in arg_match.group(1).split(',') if a.strip()]
        end_re = re.compile(r'(?i)^\s*end\s+' + kind + r'\s+' + re.escape(proc_name) + r'\b')
        end = i + 1
        while end < len(lines) and not end_re.match(lines[end]):
            end += 1
        body = lines[i + 1:end]
        for arg in args:
            matches = []
            for offset, declaration in enumerate(body, i + 2):
                if '::' not in declaration:
                    continue
                left, right = declaration.split('::', 1)
                right_code = right.split('!')[0].strip()
                entity = re.match(r'([A-Za-z_]\w*)', right_code)
                if entity and entity.group(1).lower() == arg.lower():
                    matches.append((offset, declaration, left, right_code))
            if not matches:
                errors.append(f'{path}:{start + 1}: dummy {arg} is not separately declared')
                continue
            offset, declaration, left, right_code = matches[0]
            if not re.search(r'(?i)\bintent\s*\(|\bvalue\b', left):
                errors.append(f'{path}:{offset}: dummy {arg} lacks INTENT or VALUE')
            if '!!' not in declaration:
                errors.append(f'{path}:{offset}: dummy {arg} lacks trailing FORD documentation')
            nesting = 0
            for char in right_code:
                if char == '(':
                    nesting += 1
                elif char == ')':
                    nesting = max(0, nesting - 1)
                elif char == ',' and nesting == 0:
                    errors.append(f'{path}:{offset}: dummy {arg} shares its declaration line')
                    break
        i = end + 1

if errors:
    print('\n'.join(errors))
    sys.exit(1)
print(f'Fortran source audit passed for {len(FILES)} maintained files.')
