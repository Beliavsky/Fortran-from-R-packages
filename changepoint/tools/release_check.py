#!/usr/bin/env python3
"""Release checks for the changepoint Fortran translation."""

from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAINTAINED = [ROOT / "src", ROOT / "test", ROOT / "example"]


def fail(message: str) -> None:
    raise SystemExit(f"release check failed: {message}")


def run(command: list[str]) -> None:
    print("+", " ".join(command))
    subprocess.run(command, cwd=ROOT, check=True)


def fortran_files() -> list[Path]:
    return sorted(path for directory in MAINTAINED for path in directory.glob("*.f90"))


def scan_sources() -> None:
    files = fortran_files()
    if not files:
        fail("no maintained Fortran sources found")

    hashes: dict[str, Path] = {}
    d_literal = re.compile(r"(?i)(?:\d+(?:\.\d*)?|\.\d+)[d][+-]?\d+")
    forbidden = [
        re.compile(r"(?i)\bdouble\s+precision\b"),
        re.compile(r"(?i)\breal\s*\*\s*8\b"),
        re.compile(r"(?i)kind\s*\(\s*0\.0d0\s*\)"),
        re.compile(r"(?i)\breal64\b"),
        re.compile(r"(?i)\biso_fortran_env\b"),
    ]

    for path in files:
        raw = path.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if digest in hashes:
            fail(f"duplicate Fortran source: {path} and {hashes[digest]}")
        hashes[digest] = path

        text = raw.decode("utf-8")
        lower = text.lower()
        if "use r_kinds" not in lower or "only : dp" not in lower:
            fail(f"{path} does not import the shared dp kind from r_kinds")
        for number, line in enumerate(text.splitlines(), start=1):
            if len(line) > 132:
                fail(f"{path}:{number} exceeds 132 columns")
            if ";" in line:
                fail(f"{path}:{number} contains a semicolon")
            if d_literal.search(line):
                fail(f"{path}:{number} contains a D-exponent real literal")
            if any(pattern.search(line) for pattern in forbidden):
                fail(f"{path}:{number} contains a forbidden real-kind form")

    public_module = (ROOT / "src" / "changepoint.f90").read_text(encoding="utf-8").lower()
    if "public :: dp" not in public_module:
        fail("public changepoint module does not re-export dp")


def scan_tree() -> None:
    forbidden_names = {
        "r.f90",
        "r_mod.f90",
        "r_kinds.f90",
        "r_linalg.f90",
    }
    forbidden_dirs = {
        "rfortran-core",
        "rfortran-linalg",
        "rfortran-arpack",
        "rfortran-compat",
        "blas",
        "lapack",
        "arpack",
        "build",
    }
    artifact_suffixes = {".o", ".obj", ".mod", ".smod", ".exe", ".zip", ".pyc"}

    for path in ROOT.rglob("*"):
        if path.is_dir() and path.name.lower() in forbidden_dirs:
            fail(f"forbidden vendored/build directory: {path.relative_to(ROOT)}")
        if not path.is_file():
            continue
        if path.name.lower() in forbidden_names:
            fail(f"forbidden copied dependency source: {path.relative_to(ROOT)}")
        if path.suffix.lower() in artifact_suffixes:
            fail(f"build/cache/archive artifact in release tree: {path.relative_to(ROOT)}")
        if "__pycache__" in path.parts:
            fail(f"cache artifact in release tree: {path.relative_to(ROOT)}")


def main() -> None:
    if shutil.which("fpm") is None:
        fail("fpm executable not found")

    run(["fpm", "build"])
    run(["fpm", "test"])
    scan_sources()
    scan_tree()
    run(["fpm", "clean", "--all"])
    scan_tree()
    print("release checks passed")


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as exc:
        sys.exit(exc.returncode)
