from pathlib import Path
import sys

def is_comment(line):
    if not line: return False
    return line[0] in 'cC*!'

def convert(inp,out):
    lines=Path(inp).read_text(errors='ignore').splitlines()
    res=[]
    # whether each code line is continuation
    code_indices=[i for i,l in enumerate(lines) if l.strip() and not is_comment(l)]
    next_code={}
    for a,b in zip(code_indices,code_indices[1:]): next_code[a]=b
    for i,l in enumerate(lines):
        if not l.strip():
            res.append('')
            continue
        if is_comment(l):
            # convert fixed-form comments to free-form comments
            if l[0]=='!': res.append(l)
            else: res.append('!'+l[1:])
            continue
        # pad to col 72 and ignore seq cols
        s=l+' '*(72-len(l)) if len(l)<72 else l
        label=s[:5].strip()
        cont=(len(s)>=6 and s[5] not in ' 0')
        stmt=s[6:72].rstrip()
        # preserve code beyond 72 only when nonblank (rare source exception)
        if len(l)>72 and l[72:].strip():
            stmt += l[72:]
        # identify if next code line continues this statement
        j=next_code.get(i)
        has_next_cont=False
        if j is not None:
            ns=lines[j]+' '*(6-len(lines[j])) if len(lines[j])<6 else lines[j]
            has_next_cont=(ns[5] not in ' 0')
        prefix=''
        if cont:
            prefix='& '
        elif label:
            prefix=label+' '
        text=prefix+stmt.lstrip() if not cont else prefix+stmt.lstrip()
        if has_next_cont:
            text=text.rstrip()+' &'
        res.append(text)
    Path(out).write_text('\n'.join(res)+'\n')

if __name__=='__main__': convert(sys.argv[1],sys.argv[2])
