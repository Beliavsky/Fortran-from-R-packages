#!/usr/bin/env python3
"""Release checks for the vars Fortran translation."""

from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FORTRAN_DIRS = (ROOT / "src", ROOT / "test", ROOT / "example")
FORBIDDEN_NAMES = {"r.f90", "r_mod.f90"}
FORBIDDEN_DIRS = {
    "build",
    "__pycache__",
    "rfortran-core",
    "rfortran-linalg",
    "rfortran-arpack",
    "rfortran-compat",
}
FORBIDDEN_SUFFIXES = {".o", ".obj", ".mod", ".smod", ".a", ".so", ".dll", ".exe", ".zip"}
UNSUFFIXED_REAL = re.compile(
    r"(?<![A-Za-z0-9_])(?:(?:\d+\.\d*|\.\d+)(?:[eE][+-]?\d+)?|\d+[eE][+-]?\d+)(?![A-Za-z0-9_])"
)
FORBIDDEN_KIND_PATTERNS = (
    re.compile(r"\bdouble\s+precision\b", re.IGNORECASE),
    re.compile(r"\breal\s*\*\s*8\b", re.IGNORECASE),
    re.compile(r"\bkind\s*\(\s*0\.0d0\s*\)", re.IGNORECASE),
    re.compile(r"(?<![A-Za-z0-9_])(?:\d+(?:\.\d*)?|\.\d+)[dD][+-]?\d+"),
)


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def run(command: list[str]) -> None:
    print("+", " ".join(command))
    completed = subprocess.run(command, cwd=ROOT, check=False)
    if completed.returncode != 0:
        fail(f"command failed with exit status {completed.returncode}: {' '.join(command)}")


def strip_comment(line: str) -> str:
    in_single = False
    in_double = False
    for i, char in enumerate(line):
        if char == "'" and not in_double:
            in_single = not in_single
        elif char == '"' and not in_single:
            in_double = not in_double
        elif char == "!" and not in_single and not in_double:
            return line[:i]
    return line


def has_code_semicolon(line: str) -> bool:
    code = strip_comment(line)
    in_single = False
    in_double = False
    for char in code:
        if char == "'" and not in_double:
            in_single = not in_single
        elif char == '"' and not in_single:
            in_double = not in_double
        elif char == ";" and not in_single and not in_double:
            return True
    return False


def scan_fortran() -> None:
    files = sorted(path for directory in FORTRAN_DIRS for path in directory.rglob("*.f90"))
    if not files:
        fail("no maintained Fortran sources found")

    hashes: dict[str, Path] = {}
    for path in files:
        data = path.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        if digest in hashes:
            fail(f"duplicate Fortran source: {path.relative_to(ROOT)} == {hashes[digest].relative_to(ROOT)}")
        hashes[digest] = path

        text = data.decode("utf-8")
        for line_number, line in enumerate(text.splitlines(), start=1):
            if len(line) > 132:
                fail(f"line longer than 132 characters: {path.relative_to(ROOT)}:{line_number}")
            if has_code_semicolon(line):
                fail(f"semicolon-separated statement: {path.relative_to(ROOT)}:{line_number}")
        for pattern in FORBIDDEN_KIND_PATTERNS:
            if pattern.search(text):
                fail(f"forbidden real-kind form in {path.relative_to(ROOT)}: {pattern.pattern}")

        for line_number, line in enumerate(text.splitlines(), start=1):
            code = strip_comment(line)
            code = re.sub(r"'[^']*'|\"[^\"]*\"", "", code)
            for match in UNSUFFIXED_REAL.finditer(code):
                token = match.group(0)
                tail = code[match.end():]
                if tail.startswith("_dp"):
                    continue
                fail(f"real literal lacks _dp suffix: {path.relative_to(ROOT)}:{line_number}: {token}")
            if re.search(r"(?<![<>=/ ])=(?!=|>)", code) or re.search(r"(?<![<>=/])=(?![= >])", code):
                fail(f"assignment/keyword '=' lacks surrounding spaces: {path.relative_to(ROOT)}:{line_number}")
            if re.search(r"[^ ](?:<=|>=|/=|==)", code) or re.search(r"(?:<=|>=|/=|==)[^ ]", code):
                fail(f"comparison operator lacks surrounding spaces: {path.relative_to(ROOT)}:{line_number}")

        if "real(" in text.lower() and "use r_kinds, only : dp" not in text.lower():
            fail(f"real-valued maintained source does not import dp from r_kinds: {path.relative_to(ROOT)}")

    print(f"Fortran scan: {len(files)} unique files")


def scan_tree() -> None:
    for path in ROOT.rglob("*"):
        rel = path.relative_to(ROOT)
        if any(part in FORBIDDEN_DIRS for part in rel.parts):
            fail(f"forbidden vendored/build directory: {rel}")
        if path.is_file():
            if path.name.lower() in FORBIDDEN_NAMES:
                fail(f"forbidden copied source: {rel}")
            if path.suffix.lower() in FORBIDDEN_SUFFIXES:
                fail(f"build/archive artifact present: {rel}")
    print("Tree scan: no vendored dependencies or build/archive artifacts")


def main() -> None:
    if ROOT.name != "vars":
        fail(f"top-level package directory must be named vars, got {ROOT.name!r}")
    scan_tree()
    scan_fortran()
    if shutil.which("fpm") is None:
        fail("fpm executable not found")

    run(["fpm", "build"])
    run(["fpm", "test"])
    scan_fortran()
    run(["fpm", "clean", "--all"])
    scan_tree()
    print("Release checks passed.")


if __name__ == "__main__":
    main()
