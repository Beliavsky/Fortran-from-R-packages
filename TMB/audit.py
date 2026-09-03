from pathlib import Path
import hashlib, re, sys
root=Path(__file__).resolve().parent
files=sorted([*root.glob('src/*.f90'),*root.glob('test/*.f90'),*root.glob('example/*.f90')])
errors=[]
# textual rules
for p in files:
    lines=p.read_text().splitlines()
    for n,line in enumerate(lines,1):
        code=line.split('!')[0]
        if ';' in code:
            errors.append(f'{p.relative_to(root)}:{n}: semicolon in code')
        low=code.lower().replace(' ','')
        if re.search(r'\bdouble\s+precision\b|real\s*\*\s*8|kind\s*\(\s*0\.0d0\s*\)|(?:\d|\.)d[+-]?\d', code, re.I):
            errors.append(f'{p.relative_to(root)}:{n}: disallowed real-kind spelling')
        # detect simple self-comparison NaN tests
        m=re.search(r'\b([a-zA-Z_]\w*)\s*(/=|==)\s*\1\b', code)
        if m:
            errors.append(f'{p.relative_to(root)}:{n}: self comparison {m.group(0)}')
        if len(line)>132:
            errors.append(f'{p.relative_to(root)}:{n}: line length {len(line)} > 132')
# duplicate source by content
seen={}
for p in files:
    h=hashlib.sha256(p.read_bytes()).hexdigest()
    if h in seen:
        errors.append(f'duplicate Fortran source: {seen[h]} and {p.relative_to(root)}')
    seen[h]=p.relative_to(root)
# dummy declaration checks (simple parser sufficient for maintained style)
proc_re=re.compile(r'^\s*(?:pure\s+)?(?:elemental\s+)?(?:recursive\s+)?(?:[\w()]+\s+)?(?:function|subroutine)\s+(\w+)\s*\(([^)]*)\)',re.I)
for p in files:
    lines=p.read_text().splitlines()
    i=0
    while i<len(lines):
        m=proc_re.match(lines[i])
        if not m:
            i+=1; continue
        name=m.group(1)
        args=[a.strip().lower() for a in m.group(2).split(',') if a.strip()]
        # handle continuation header crudely
        if '&' in m.group(2):
            args=[]
        end_re=re.compile(r'^\s*end\s+(?:function|subroutine)\b',re.I)
        block=[]; j=i+1
        while j<len(lines) and not end_re.match(lines[j]):
            block.append((j+1,lines[j])); j+=1
        for arg in args:
            decl=[]
            for ln,line in block:
                if '::' in line and re.search(r'\b'+re.escape(arg)+r'\b',line,re.I):
                    decl.append((ln,line))
            if not decl:
                errors.append(f'{p.relative_to(root)}:{i+1}: dummy {arg} in {name} lacks declaration')
                continue
            ln,line=decl[0]
            pre=line.split('::',1)[0].lower()
            rhs=line.split('::',1)[1].split('!',1)[0]
            # one dummy per declaration: reject commas separating names (shape commas allowed before :: only)
            if re.search(r'(^|,)\s*'+re.escape(arg)+r'\s*(?:,|$)',rhs.strip(),re.I) and ',' in rhs.strip():
                pass
            if 'intent(' not in pre and 'value' not in pre and not pre.strip().startswith('procedure'):
                errors.append(f'{p.relative_to(root)}:{ln}: dummy {arg} lacks INTENT/VALUE')
            if '!!' not in line:
                errors.append(f'{p.relative_to(root)}:{ln}: dummy {arg} lacks trailing FORD comment')
        i=j+1
# dependency-copy / build product checks
for p in root.rglob('*'):
    if not p.is_file(): continue
    rel=p.relative_to(root)
    if any(part in {'build','build_manual','.git','__pycache__'} for part in rel.parts):
        continue
    if p.suffix.lower() in {'.o','.mod','.smod','.exe','.dll','.so','.a','.zip','.pyc'}:
        errors.append(f'build/archive product present: {rel}')
    if p.name.lower() in {'r.f90','r_mod.f90'}:
        errors.append(f'copied shared dependency source: {rel}')
print(f'Checked {len(files)} maintained Fortran files.')
if errors:
    print('\n'.join(errors))
    sys.exit(1)
print('Audit passed: no prohibited semicolons, kind spellings, self-comparison NaN tests, duplicate source, long lines, or packaged build products.')
