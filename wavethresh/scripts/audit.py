#!/usr/bin/env python3
from pathlib import Path
import re, sys, hashlib, tomllib
root=Path(__file__).resolve().parents[1]
errors=[]
fort=list(root.glob('src/*.f90'))+list(root.glob('test/*.f90'))+list(root.glob('app/*.f90'))+list(root.glob('example/*.f90'))
# ASCII, length, forbidden patterns, semicolon executable chains.
for p in fort:
    raw=p.read_bytes()
    try: text=raw.decode('ascii')
    except UnicodeDecodeError: errors.append(f'non-ASCII: {p.relative_to(root)}'); continue
    for ln,line in enumerate(text.splitlines(),1):
        if len(line)>132: errors.append(f'line >132: {p.relative_to(root)}:{ln}:{len(line)}')
        code=line.split('!')[0]
        code_no_strings=re.sub(r'\"(?:[^\"]|\"\")*\"|\'(?:[^\']|\'\')*\'', '', code)
        if ';' in code_no_strings: errors.append(f'semicolon statement: {p.relative_to(root)}:{ln}')
    patterns=[
        (r'(?i)\bdouble\s+precision\b','double precision'),
        (r'(?i)\breal\s*\*\s*8\b','real*8'),
        (r'(?i)kind\s*\(\s*0\.0d0\s*\)','kind(0.0d0)'),
        (r'(?i)(?<![A-Za-z0-9_])[0-9]+(?:\.[0-9]*)?[dD][+-]?[0-9]+','D exponent'),
        (r'\b([A-Za-z][A-Za-z0-9_]*)\s*/=\s*\1\b','self-comparison NaN test')]
    for pat,label in patterns:
        if re.search(pat,text): errors.append(f'{label}: {p.relative_to(root)}')
# duplicate maintained Fortran content
seen={}
for p in fort:
    h=hashlib.sha256(p.read_bytes()).hexdigest()
    if h in seen: errors.append(f'duplicate Fortran source: {seen[h]} and {p.relative_to(root)}')
    seen[h]=p.relative_to(root)
# no dependency source in maintained src
for forbidden in ['r_kinds.f90','waveslim_transform_1d.f90','waveslim_transform_nd.f90','waveslim_types.f90']:
    if (root/'src'/forbidden).exists(): errors.append(f'vendored dependency source: src/{forbidden}')
# no build artifacts
bad_suffix={'.o','.obj','.mod','.smod','.exe','.dll','.so','.a','.pyc','.zip'}
for p in root.rglob('*'):
    if p.is_file() and p.suffix.lower() in bad_suffix: errors.append(f'build/archive artifact in tree: {p.relative_to(root)}')
for p in root.rglob('__pycache__'):
    errors.append(f'cache directory: {p.relative_to(root)}')
# dummy argument declarations: join procedure signature continuation lines, then inspect source lines.
proc_start=re.compile(r'(?i)^\s*(?:(?:pure|elemental|recursive|impure|module)\s+)*(?:function|subroutine)\s+([A-Za-z][A-Za-z0-9_]*)\s*\(')
for p in fort:
    lines=p.read_text().splitlines()
    i=0
    while i < len(lines):
        m=proc_start.match(lines[i])
        if not m:
            i+=1; continue
        pname=m.group(1)
        sig=lines[i].split('!')[0]
        j=i
        while ')' not in sig and j+1<len(lines):
            j+=1; sig += ' '+lines[j].split('!')[0].replace('&',' ')
        mm=re.search(re.escape(pname)+r'\s*\((.*?)\)',sig,re.I)
        if not mm:
            i=j+1; continue
        args=[a.strip() for a in mm.group(1).split(',') if a.strip()]
        # strip possible keyword-ish whitespace; signatures here are simple names.
        args=[re.sub(r'[^A-Za-z0-9_].*$','',a) for a in args]
        # search until contains/end procedure (or 120 lines max)
        search_lines=lines[j+1:min(len(lines),j+140)]
        for arg in args:
            decls=[]
            for off,line in enumerate(search_lines,j+2):
                code=line.split('!')[0]
                if re.match(r'(?i)^\s*(contains|end\s+(function|subroutine))\b',code): break
                if '::' in code:
                    right=code.split('::',1)[1].strip()
                    # first declared entity name; each dummy must be alone
                    nm=re.match(r'([A-Za-z][A-Za-z0-9_]*)',right)
                    if nm and nm.group(1).lower()==arg.lower(): decls.append((off,line,code))
            if len(decls)!=1:
                errors.append(f'dummy {arg} declaration count {len(decls)}: {p.relative_to(root)}:{i+1} {pname}')
                continue
            off,line,code=decls[0]
            left=code.split('::',1)[0]
            if not (re.search(r'(?i)\bintent\s*\(',left) or re.search(r'(?i)\bvalue\b',left)):
                errors.append(f'dummy {arg} lacks INTENT/VALUE: {p.relative_to(root)}:{off}')
            if '!!' not in line:
                errors.append(f'dummy {arg} lacks FORD comment: {p.relative_to(root)}:{off}')
            right=code.split('::',1)[1].strip()
            # detect multiple top-level declared entities after shape/initialization; commas inside parens allowed.
            depth=0; topcomma=False
            for ch in right:
                if ch=='(': depth+=1
                elif ch==')': depth=max(0,depth-1)
                elif ch==',' and depth==0: topcomma=True
            if topcomma:
                errors.append(f'dummy shares declaration line: {p.relative_to(root)}:{off}')
        i=j+1
# manifest/README consistency
try:
    manifest=tomllib.loads((root/'fpm.toml').read_text())
    tr=manifest['extra']['translation']; entries=tr['function']
    mapped=tr['r_functions_translated']; total=tr['r_functions_total']; frac=tr['frac_functions_translated']
    if mapped!=len(entries): errors.append(f'manifest mapped count {mapped} != entries {len(entries)}')
    if abs(frac-mapped/total)>1e-14: errors.append('manifest fraction inconsistent')
    names=[e['r_name'] for e in entries]
    if len(names)!=len(set(names)): errors.append('duplicate r_name mapping entries')
    if len(tr['untranslated_r_functions']) != total-mapped: errors.append('untranslated list length inconsistent')
    readme=(root/'README.md').read_text()
    pct=100*mapped/total
    if f'{mapped} of {total} ({pct:.1f}%)' not in readme: errors.append('README coverage summary inconsistent')
    row_count=sum(1 for line in readme.splitlines() if line.startswith('| `'))
    if row_count!=mapped: errors.append(f'README mapping rows {row_count} != {mapped}')
except Exception as e:
    errors.append(f'manifest/README audit exception: {e}')
if errors:
    print('AUDIT: FAIL')
    for e in errors: print(' -',e)
    sys.exit(1)
print(f'AUDIT: PASS ({len(fort)} Fortran files)')
