#!/usr/bin/env python3
from pathlib import Path
import re
import sys

root = Path(__file__).resolve().parents[1]
files = sorted([*root.joinpath('src').glob('*.f90'), *root.joinpath('test').glob('*.f90'), *root.joinpath('example').glob('*.f90')])
errors = []

header_re = re.compile(r'^\s*(?:pure\s+|recursive\s+|elemental\s+)*(?:[\w()=*,: ]+\s+)?(subroutine|function)\s+(\w+)\s*\(([^)]*)\)', re.I)
end_re = re.compile(r'^\s*end\s+(subroutine|function)\b', re.I)

def split_top_level(text):
    out=[]; cur=[]; depth=0
    for ch in text:
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        if ch == ',' and depth == 0:
            out.append(''.join(cur).strip()); cur=[]
        else:
            cur.append(ch)
    if cur: out.append(''.join(cur).strip())
    return [x for x in out if x]

for path in files:
    lines = path.read_text().splitlines()
    if any(len(line) > 132 for line in lines):
        for i,line in enumerate(lines,1):
            if len(line)>132:
                errors.append(f'{path.name}:{i}: line length {len(line)} > 132')
    if path.name != 'learnbayes_kinds.f90' and re.search(r'\breal64\b', path.read_text(), re.I):
        errors.append(f'{path.name}: real64 must be defined/imported only in learnbayes_kinds.f90')
    text = '\n'.join(line.split('!',1)[0] for line in path.read_text().splitlines())
    for pat,label in [(r'\bdouble\s+precision\b','double precision'),(r'\breal\s*\*\s*8\b','real*8'),
                      (r'\bkind\s*\(\s*0\.0d0\s*\)','kind(0.0d0)'),(r'\d(?:\.\d*)?[dD][+-]?\d+','D exponent')]:
        if re.search(pat,text,re.I): errors.append(f'{path.name}: forbidden {label}')
    i=0
    while i < len(lines):
        m=header_re.match(lines[i])
        if not m or lines[i].lstrip().lower().startswith('end '):
            i+=1; continue
        kind,name,args_text=m.groups()
        args=[a.strip().lower() for a in args_text.split(',') if a.strip()]
        # find end of procedure
        j=i+1
        while j < len(lines) and not end_re.match(lines[j]):
            j+=1
        body=lines[i+1:j]
        for arg in args:
            found=[]
            for offset,line in enumerate(body,i+2):
                code=line.split('!!',1)[0]
                if '::' not in code: continue
                rhs=code.split('::',1)[1]
                names=split_top_level(rhs)
                base_names=[]
                for n in names:
                    nm=re.match(r'\s*([A-Za-z_]\w*)',n)
                    if nm: base_names.append(nm.group(1).lower())
                if arg in base_names:
                    found.append((offset,line,names))
            if len(found)!=1:
                errors.append(f'{path.name}:{i+1}: {name} dummy {arg} has {len(found)} declarations')
                continue
            line_no,line,names=found[0]
            low=line.lower()
            if 'intent(' not in low and not re.search(r'\bvalue\b', low.split('::',1)[0]):
                errors.append(f'{path.name}:{line_no}: dummy {arg} lacks INTENT or VALUE')
            if '!!' not in line:
                errors.append(f'{path.name}:{line_no}: dummy {arg} lacks trailing FORD !! comment')
            else:
                comment=line.split('!!',1)[1].strip()
                if len(comment)<12:
                    errors.append(f'{path.name}:{line_no}: dummy {arg} FORD comment is not meaningful enough')
            if len(names) != 1:
                errors.append(f'{path.name}:{line_no}: dummy {arg} shares a declaration with another entity')
        i=j+1

if errors:
    print('\n'.join(errors))
    sys.exit(1)
print(f'PASS: {len(files)} maintained Fortran files satisfy source rules.')
