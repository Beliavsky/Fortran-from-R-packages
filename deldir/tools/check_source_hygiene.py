from pathlib import Path
import re, sys
src = Path(__file__).resolve().parents[1] / "src"
files = sorted(src.glob("*.f90"))
errors=[]
text="\n".join(p.read_text() for p in files)
checks=[
    (r"\bdouble\s+precision\b", "double precision"),
    (r"\breal\s*\*\s*8\b", "real*8"),
    (r"kind\s*\(\s*0\.0d0\s*\)", "kind(0.0d0)"),
    (r"(?i)(?<![A-Za-z0-9_])[0-9.]+d[+-]?[0-9]+", "D exponent literal"),
]
for pat,label in checks:
    if re.search(pat,text,re.I): errors.append(label)
real64_count=sum(p.read_text().count("real64") for p in files)
if real64_count != 2: # import plus dp definition, both in the one kinds file
    errors.append(f"unexpected real64 occurrence count: {real64_count}")
for p in files:
    if p.name != "deldir_kinds.f90" and "real64" in p.read_text():
        errors.append(f"real64 outside kinds module: {p.name}")
if errors:
    print("source hygiene: FAIL")
    for e in errors: print(" -",e)
    sys.exit(1)
print(f"source hygiene: PASS ({len(files)} maintained Fortran files)")
