#!/usr/bin/env python3
"""Static policy checks for the maintained Fortran translation."""

from __future__ import annotations

import hashlib
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAINTAINED_DIRS = [ROOT / "src", ROOT / "test", ROOT / "example", ROOT / "validation"]
FORTRAN = sorted(
    path
    for base in MAINTAINED_DIRS
    for path in base.rglob("*")
    if path.is_file() and path.suffix.lower() in {".f90", ".f95", ".f03", ".f08", ".f", ".for"}
)


def code_before_comment(line: str) -> str:
    """Return Fortran code before an unquoted exclamation-mark comment."""
    out: list[str] = []
    quote: str | None = None
    i = 0
    while i < len(line):
        ch = line[i]
        if quote is None:
            if ch in {"'", '"'}:
                quote = ch
                out.append(ch)
            elif ch == "!":
                break
            else:
                out.append(ch)
        else:
            out.append(ch)
            if ch == quote:
                if i + 1 < len(line) and line[i + 1] == quote:
                    out.append(line[i + 1])
                    i += 1
                else:
                    quote = None
        i += 1
    return "".join(out)


def contains_unquoted_semicolon(code: str) -> bool:
    quote: str | None = None
    i = 0
    while i < len(code):
        ch = code[i]
        if quote is None:
            if ch in {"'", '"'}:
                quote = ch
            elif ch == ";":
                return True
        elif ch == quote:
            if i + 1 < len(code) and code[i + 1] == quote:
                i += 1
            else:
                quote = None
        i += 1
    return False


def split_top_level_commas(text: str) -> list[str]:
    parts: list[str] = []
    start = 0
    depth = 0
    quote: str | None = None
    i = 0
    while i < len(text):
        ch = text[i]
        if quote is None:
            if ch in {"'", '"'}:
                quote = ch
            elif ch in "([":
                depth += 1
            elif ch in ")]":
                depth = max(0, depth - 1)
            elif ch == "," and depth == 0:
                parts.append(text[start:i].strip())
                start = i + 1
        elif ch == quote:
            if i + 1 < len(text) and text[i + 1] == quote:
                i += 1
            else:
                quote = None
        i += 1
    parts.append(text[start:].strip())
    return [part for part in parts if part]


def fail(errors: list[str], path: Path, lineno: int | None, message: str) -> None:
    rel = path.relative_to(ROOT)
    loc = f"{rel}:{lineno}" if lineno is not None else str(rel)
    errors.append(f"{loc}: {message}")


def main() -> int:
    errors: list[str] = []

    # Maintained Fortran must all be free-form .f90.
    for path in FORTRAN:
        if path.suffix.lower() != ".f90":
            fail(errors, path, None, "maintained Fortran source is not free-form .f90")

    dp_defs = 0
    hashes: dict[str, Path] = {}
    self_compare = re.compile(r"\b([a-z_]\w*(?:%[a-z_]\w*)?)\b\s*(?:/=|\.ne\.)\s*\b\1\b", re.I)
    d_exp = re.compile(r"(?<![a-z0-9_])(?:\d+(?:\.\d*)?|\.\d+)[dD][+-]?\d+")
    bad_kind = re.compile(r"\bdouble\s+precision\b|\breal\s*\*\s*8\b|kind\s*\(\s*0\.0[dD]0\s*\)", re.I)

    for path in FORTRAN:
        raw = path.read_text(encoding="utf-8")
        digest = hashlib.sha256(raw.encode("utf-8")).hexdigest()
        if digest in hashes:
            fail(errors, path, None, f"duplicate Fortran content also appears in {hashes[digest].relative_to(ROOT)}")
        else:
            hashes[digest] = path

        for lineno, line in enumerate(raw.splitlines(), 1):
            code = code_before_comment(line)
            if contains_unquoted_semicolon(code):
                fail(errors, path, lineno, "semicolon-separated statement is not allowed")
            if bad_kind.search(code) or d_exp.search(code):
                fail(errors, path, lineno, "disallowed real kind or D-exponent literal")
            if self_compare.search(code):
                fail(errors, path, lineno, "self-comparison idiom detected; use ieee_is_nan for NaN tests")
            if re.search(r"\binteger\s*,\s*parameter\s*,\s*public\s*::\s*dp\s*=\s*real64\b", code, re.I):
                dp_defs += 1

            # Each dummy declaration must have INTENT/VALUE, one declarator, and a FORD trailing comment.
            low = code.lower()
            prefix = low.split("::", 1)[0]
            if "::" in code and ("intent(" in prefix or re.search(r"(?:^|,)\s*value(?:\s*,|\s*$)", prefix)):
                decls = split_top_level_commas(code.split("::", 1)[1])
                if len(decls) != 1:
                    fail(errors, path, lineno, "declare each dummy argument on its own line")
                if "!!" not in line:
                    fail(errors, path, lineno, "dummy argument lacks trailing FORD !! documentation")
                else:
                    doc = line.split("!!", 1)[1].strip()
                    if len(doc) < 12 or doc.lower() in {"input value", "output value", "argument"}:
                        fail(errors, path, lineno, "dummy FORD documentation is not meaningful")

    # Verify that every procedure dummy has an INTENT or VALUE declaration.
    procedure_start = re.compile(r"\b(subroutine|function)\s+([a-z_]\w*)\s*\((.*?)\)", re.I)
    for path in FORTRAN:
        lines = path.read_text(encoding="utf-8").splitlines()
        i = 0
        while i < len(lines):
            code = code_before_comment(lines[i]).rstrip()
            if re.search(r"\bend\s+(?:subroutine|function)\b", code, re.I):
                i += 1
                continue
            stmt = code
            end_stmt_line = i
            while stmt.rstrip().endswith("&") and end_stmt_line + 1 < len(lines):
                stmt = stmt.rstrip()[:-1] + " " + code_before_comment(lines[end_stmt_line + 1]).strip()
                end_stmt_line += 1
            match = procedure_start.search(stmt)
            if match:
                kind, name, arg_text = match.groups()
                args = [arg.strip().lower() for arg in split_top_level_commas(arg_text) if arg.strip()]
                declared: set[str] = set()
                k = end_stmt_line + 1
                end_re = re.compile(rf"\bend\s+{kind}\b(?:\s+{re.escape(name)}\b)?", re.I)
                while k < len(lines):
                    body_code = code_before_comment(lines[k])
                    if end_re.search(body_code):
                        break
                    if "::" in body_code:
                        prefix, rhs = body_code.split("::", 1)
                        if "intent(" in prefix.lower() or re.search(
                            r"(?:^|,)\s*value(?:\s*,|\s*$)", prefix.lower()
                        ):
                            decl = re.match(r"\s*([a-z_]\w*)", rhs, re.I)
                            if decl:
                                declared.add(decl.group(1).lower())
                    k += 1
                for arg in args:
                    if arg not in declared:
                        fail(errors, path, i + 1, f"dummy argument {arg} lacks explicit INTENT or VALUE")
                i = end_stmt_line
            i += 1

    if dp_defs != 1:
        errors.append(f"expected exactly one dp = real64 definition, found {dp_defs}")

    # Packaging/build-product policy.
    forbidden_suffixes = {".o", ".obj", ".mod", ".smod", ".exe", ".a", ".so", ".dll", ".pyc", ".zip"}
    for path in ROOT.rglob("*"):
        if path.is_file() and path.suffix.lower() in forbidden_suffixes:
            fail(errors, path, None, "generated/build product must not be packaged")
        if path.is_dir() and path.name in {"build", "__pycache__", ".pytest_cache"}:
            fail(errors, path, None, "build/cache directory must not be packaged")

    # Manifest/coverage consistency.
    manifest = tomllib.loads((ROOT / "fpm.toml").read_text(encoding="utf-8"))
    tr = manifest["extra"]["translation"]
    mappings = tr.get("function", [])
    total = tr["r_functions_total"]
    translated = tr["r_functions_translated"]
    fraction = tr["frac_functions_translated"]
    untranslated = tr["untranslated_r_functions"]
    if translated != len(mappings):
        errors.append("fpm.toml r_functions_translated does not equal the number of mapping entries")
    if total <= 0 or abs(fraction - translated / total) > 1.0e-12:
        errors.append("fpm.toml frac_functions_translated is inconsistent")
    if total != translated + len(untranslated):
        errors.append("fpm.toml translated/untranslated counts do not sum to r_functions_total")
    r_names = [entry["r_name"] for entry in mappings]
    if len(r_names) != len(set(r_names)):
        errors.append("fpm.toml contains duplicate R-function mapping entries")

    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    expected_summary = f"Coverage is **{translated} of {total} ({round(100 * fraction):d}%)**"
    if expected_summary not in readme:
        errors.append("README Translation coverage summary disagrees with fpm.toml")
    for entry in mappings:
        marker = f"| `{entry['r_name']}` |"
        if marker not in readme:
            errors.append(f"README translation table is missing {entry['r_name']}")

    # No external BLAS/LAPACK/system-library linkage or vendored shared dependency source.
    manifest_text = (ROOT / "fpm.toml").read_text(encoding="utf-8").lower()
    if re.search(r"\blink\s*=", manifest_text) or "-llapack" in manifest_text or "-lblas" in manifest_text:
        errors.append("fpm.toml contains disallowed system-library linkage")
    forbidden_names = {"r.f90", "r_mod.f90"}
    for path in ROOT.rglob("*"):
        if path.is_file() and path.name.lower() in forbidden_names:
            fail(errors, path, None, "copied shared compatibility source is disallowed")

    if errors:
        print("Static audit FAILED:")
        for item in errors:
            print(f"  - {item}")
        return 1

    print(f"Static audit passed for {len(FORTRAN)} maintained Fortran source files.")
    print(f"Translation coverage: {translated}/{total} computational R functions.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
