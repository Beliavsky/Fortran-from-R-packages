from pathlib import Path
import re
import tomllib

root = Path(__file__).resolve().parents[1]
sources = sorted(root.joinpath("src").glob("*.f90"))
errors = []
texts = {}
for path in sources:
    text = path.read_text(encoding="ascii")
    texts[path] = text
    for no, line in enumerate(text.splitlines(), 1):
        code = line.split("!", 1)[0]
        if ";" in code:
            errors.append(f"{path.name}:{no}: semicolon in executable/declaration code")
        if len(line) > 132:
            errors.append(f"{path.name}:{no}: line exceeds 132 columns")
        if re.search(r"\b(?:double\s+precision|real\s*\*\s*8)\b", code, re.I):
            errors.append(f"{path.name}:{no}: legacy real declaration")
        if re.search(r"\d(?:\.\d*)?[dD][+-]?\d+", code):
            errors.append(f"{path.name}:{no}: D exponent literal")
        if re.search(r"\b(\w+)\s*/=\s*\1\b", code):
            errors.append(f"{path.name}:{no}: self-comparison NaN test")
        if "::" in code and "intent(" in code.lower() and "!!" not in line:
            errors.append(f"{path.name}:{no}: dummy declaration lacks trailing FORD comment")

hashes = {}
import hashlib
for path, text in texts.items():
    h = hashlib.sha256(text.encode("ascii")).hexdigest()
    if h in hashes:
        errors.append(f"duplicate Fortran content: {path.name} and {hashes[h].name}")
    hashes[h] = path

manifest = tomllib.loads(root.joinpath("fpm.toml").read_text())
translation = manifest["extra"]["translation"]
mappings = translation["function"]
if translation["r_functions_total"] != 1 or translation["r_functions_translated"] != 1:
    errors.append("translation counts are not 1/1")
if len(mappings) != 1 or mappings[0]["r_name"] != "gamm4":
    errors.append("translation mapping does not contain exactly gamm4")
readme = root.joinpath("README.md").read_text()
if "1 of 1 (100%)" not in readme:
    errors.append("README coverage summary disagrees with manifest")

bad_suffixes = {".o", ".mod", ".smod", ".exe", ".pyc", ".zip"}
for path in root.rglob("*"):
    if path.is_file() and path.suffix.lower() in bad_suffixes:
        errors.append(f"build/archive product in package tree: {path.relative_to(root)}")
    if "__pycache__" in path.parts:
        errors.append(f"cache directory in package tree: {path.relative_to(root)}")

if errors:
    print("SOURCE AUDIT FAILED")
    for error in errors:
        print(error)
    raise SystemExit(1)
print(f"SOURCE AUDIT PASSED: {len(sources)} Fortran files; 1/1 coverage mapping consistent.")
