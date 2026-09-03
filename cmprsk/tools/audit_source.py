#!/usr/bin/env python3
from pathlib import Path
import re, sys
root=Path(__file__).resolve().parents[1]
files=sorted(root.rglob('*.f90'))
issues=[]
seen_hash={}

def codepart(line): return line.split('!',1)[0]
def split_top(s):
    out=[];cur='';depth=0
    for ch in s:
        if ch=='(': depth+=1
        elif ch==')': depth-=1
        if ch==',' and depth==0:
            out.append(cur.strip());cur=''
        else: cur+=ch
    if cur.strip():out.append(cur.strip())
    return out

for p in files:
    rel=p.relative_to(root)
    text=p.read_text()
    lines=text.splitlines()
    # line length / semicolon / forbidden spellings
    for i,line in enumerate(lines,1):
        if len(line)>132: issues.append(f'{rel}:{i}: line length {len(line)}')
        cp=codepart(line)
        if ';' in cp: issues.append(f'{rel}:{i}: semicolon in code')
        if re.search(r'(?i)\bdouble\s+precision\b|\breal\s*\*\s*8\b|kind\s*\(\s*0\.0d0\s*\)|\b\d+(?:\.\d*)?[dD][+-]?\d+',cp):
            issues.append(f'{rel}:{i}: forbidden real-kind/D exponent form')
        m=re.search(r'\b([A-Za-z_]\w*)\s*/=\s*\1\b',cp,re.I)
        if m: issues.append(f'{rel}:{i}: self-comparison NaN idiom')
    # duplicate exact source content
    import hashlib
    h=hashlib.sha256(text.encode()).hexdigest()
    if h in seen_hash: issues.append(f'{rel}: duplicate of {seen_hash[h]}')
    seen_hash[h]=rel
    # Build logical procedure headers
    logical=[];buf='';start=0
    for i,line in enumerate(lines,1):
        cp=codepart(line).rstrip()
        if not buf:
            buf=cp;start=i
        else:
            buf += ' ' + cp.lstrip('&').strip()
        if cp.endswith('&'):
            buf=buf[:-1].rstrip();continue
        logical.append((start,buf));buf=''
    if buf: logical.append((start,buf))
    proc_dummies=[]
    for ln,s in logical:
        mm=re.search(r'(?i)^\s*(?:pure\s+|elemental\s+|recursive\s+|impure\s+)*(?:[\w()]+\s+)?(?:subroutine|function)\s+(\w+)\s*\((.*?)\)',s)
        if mm:
            args=[a.strip() for a in split_top(mm.group(2)) if a.strip()]
            proc_dummies.extend((a.lower(),mm.group(1),ln) for a in args)
    decl={}
    for i,line in enumerate(lines,1):
        cp=codepart(line)
        if '::' not in cp: continue
        left,right=cp.split('::',1)
        if not re.search(r'(?i)\b(intent\s*\(|value\b)',left): continue
        names=split_top(right)
        if len(names)!=1:
            issues.append(f'{rel}:{i}: shared dummy declaration: {right.strip()}')
        for nm in names:
            base=re.match(r'\s*(\w+)',nm)
            if base: decl[base.group(1).lower()]=(i,left,line)
    for arg,proc,ln in proc_dummies:
        if arg not in decl:
            issues.append(f'{rel}:{ln}: dummy {arg} of {proc} lacks explicit INTENT/VALUE declaration')
        else:
            i,left,line=decl[arg]
            if '!!' not in line:
                issues.append(f'{rel}:{i}: dummy {arg} lacks trailing FORD !! comment')

# Build artifacts / archives inside package
for p in root.rglob('*'):
    if p.is_file() and (p.suffix.lower() in {'.o','.mod','.exe','.zip','.a','.so','.dll','.pyc'} or '__pycache__' in p.parts or 'build' in p.parts):
        issues.append(f'{p.relative_to(root)}: build/archive artifact in package')

if issues:
    print('\n'.join(issues)); sys.exit(1)
print(f'AUDIT PASS: {len(files)} Fortran files, 0 issues')
