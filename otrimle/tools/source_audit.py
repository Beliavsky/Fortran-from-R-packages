"""Static policy and coverage audit for the otrimle Fortran translation."""

from __future__ import annotations

import hashlib
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FORTRAN_FILES = sorted(
    list((ROOT / "src").glob("*.f90"))
    + list((ROOT / "test").glob("*.f90"))
    + list((ROOT / "example").glob("*.f90"))
    + list((ROOT / "tools").glob("*.f90"))
)


def main() -> int:
    errors: list[str] = []
    legacy_patterns = [
        (re.compile(r"\bdouble\s+precision\b", re.I), "double precision"),
        (re.compile(r"\breal\s*\*\s*8\b", re.I), "real*8"),
        (re.compile(r"\bkind\s*\(\s*0\.0d0\s*\)", re.I), "kind(0.0d0)"),
        (
            re.compile(r"(?<![A-Za-z0-9_])(?:\d+(?:\.\d*)?|\.\d+)[dD][+-]?\d+"),
            "D exponent",
        ),
    ]
    self_compare = re.compile(
        r"\b([A-Za-z_]\w*)\s*/=\s*\1\b|\b([A-Za-z_]\w*)\s*==\s*\2\b"
    )

    for path in FORTRAN_FILES:
        for line_number, line in enumerate(path.read_text().splitlines(), 1):
            if len(line) > 132:
                errors.append(f"{path.relative_to(ROOT)}:{line_number}: line length {len(line)}")
            code = line.split("!")[0]
            if ";" in code:
                errors.append(f"{path.relative_to(ROOT)}:{line_number}: executable semicolon")
            for pattern, label in legacy_patterns:
                if pattern.search(code):
                    errors.append(f"{path.relative_to(ROOT)}:{line_number}: legacy {label}")
            if self_compare.search(code):
                errors.append(f"{path.relative_to(ROOT)}:{line_number}: self-comparison idiom")

    dp_definitions = 0
    for path in (ROOT / "src").glob("*.f90"):
        for line in path.read_text().splitlines():
            if re.search(r"\binteger\s*,[^:]*\bparameter\b[^:]*::\s*dp\b", line, re.I):
                dp_definitions += 1
    if dp_definitions != 1:
        errors.append(f"dp definitions: {dp_definitions}, expected 1")

    seen: dict[str, Path] = {}
    for path in FORTRAN_FILES:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest in seen:
            errors.append(
                "duplicate Fortran content: "
                f"{path.relative_to(ROOT)} == {seen[digest].relative_to(ROOT)}"
            )
        seen[digest] = path

    procedure_start = re.compile(
        r"^\s*(?:pure\s+|elemental\s+|recursive\s+|module\s+)*\s*"
        r"(?:subroutine|(?:[a-z][\w()=*,:\s]*\s+)?function)\s+"
        r"(\w+)\s*\(([^)]*)\)",
        re.I,
    )
    declaration = re.compile(r"^\s*([^!]*?)::\s*([^!]+?)(?:\s*!!(.*))?$", re.I)

    for path in FORTRAN_FILES:
        lines = path.read_text().splitlines()
        index = 0
        while index < len(lines):
            match = procedure_start.match(lines[index])
            if not match:
                index += 1
                continue
            procedure = match.group(1)
            arguments = [
                arg.strip().lower() for arg in match.group(2).split(",") if arg.strip()
            ]
            declared: set[str] = set()
            scan = index + 1
            while scan < len(lines):
                if re.match(r"\s*end\s+(subroutine|function)\b", lines[scan], re.I):
                    break
                decl_match = declaration.match(lines[scan])
                if decl_match:
                    attributes = decl_match.group(1).lower()
                    variable_part = decl_match.group(2).strip()
                    comment = decl_match.group(3)
                    items: list[str] = []
                    depth = 0
                    token = ""
                    for char in variable_part + ",":
                        if char == "(":
                            depth += 1
                        elif char == ")":
                            depth -= 1
                        if char == "," and depth == 0:
                            items.append(token.strip())
                            token = ""
                        else:
                            token += char
                    dummy_names: list[str] = []
                    for item in items:
                        var_match = re.match(r"([a-z_]\w*)", item, re.I)
                        if var_match and var_match.group(1).lower() in arguments:
                            dummy_names.append(var_match.group(1).lower())
                    for dummy in dummy_names:
                        if len(dummy_names) > 1:
                            errors.append(
                                f"{path.relative_to(ROOT)}:{scan + 1}: "
                                f"dummy {dummy} shares declaration line"
                            )
                        if "intent(" not in attributes and not re.search(
                            r"\bvalue\b", attributes
                        ):
                            errors.append(
                                f"{path.relative_to(ROOT)}:{scan + 1}: "
                                f"dummy {dummy} lacks INTENT/VALUE"
                            )
                        if comment is None or not comment.strip():
                            errors.append(
                                f"{path.relative_to(ROOT)}:{scan + 1}: "
                                f"dummy {dummy} lacks trailing !! comment"
                            )
                        declared.add(dummy)
                scan += 1
            for dummy in arguments:
                if dummy not in declared:
                    errors.append(
                        f"{path.relative_to(ROOT)}:{index + 1}: procedure {procedure} "
                        f"dummy {dummy} not audited/declared"
                    )
            index = scan + 1

    with (ROOT / "fpm.toml").open("rb") as stream:
        manifest = tomllib.load(stream)
    translation = manifest["extra"]["translation"]
    functions = translation.get("function", [])
    total = translation["r_functions_total"]
    translated = translation["r_functions_translated"]
    fraction = translation["frac_functions_translated"]
    untranslated = translation["untranslated_r_functions"]
    if len(functions) != translated:
        errors.append(f"manifest mapping entries {len(functions)} != {translated}")
    if translated + len(untranslated) != total:
        errors.append(
            f"manifest translated + untranslated mismatch: "
            f"{translated} + {len(untranslated)} != {total}"
        )
    if abs(fraction - translated / total) > 1.0e-12:
        errors.append(f"manifest fraction mismatch: {fraction} != {translated / total}")

    readme = (ROOT / "README.md").read_text()
    expected = f"**{translated} of {total} ({100.0 * translated / total:.1f}%)**"
    if expected not in readme:
        errors.append(f"README coverage text does not contain {expected}")

    forbidden_names = {"r.f90", "r_mod.f90"}
    for path in ROOT.rglob("*"):
        if path.is_file() and path.name.lower() in forbidden_names:
            errors.append(f"forbidden copied source: {path.relative_to(ROOT)}")
        if path.is_file() and path.suffix.lower() in {".o", ".mod", ".smod", ".exe", ".obj", ".a"}:
            errors.append(f"build product in package tree: {path.relative_to(ROOT)}")
        if path.is_file() and path.suffix.lower() == ".zip":
            errors.append(f"ZIP nested in package tree: {path.relative_to(ROOT)}")

    if errors:
        print("SOURCE AUDIT FAILED")
        for error in errors:
            print(error)
        return 1

    print(
        f"SOURCE AUDIT PASSED: {len(FORTRAN_FILES)} Fortran files; "
        f"{translated}/{total} coverage mapping consistent."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
