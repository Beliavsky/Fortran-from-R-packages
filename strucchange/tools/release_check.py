#!/usr/bin/env python3
"""Release checks for the strucchange Fortran translation."""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FORTRAN = sorted(ROOT.rglob("*.f90"))


def run(command: list[str]) -> None:
    print("+", " ".join(command))
    subprocess.run(command, cwd=ROOT, check=True)


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def scan_sources() -> None:
    if not FORTRAN:
        fail("no Fortran source files found")

    hashes: dict[str, Path] = {}
    forbidden_kind = [
        re.compile(r"\bdouble\s+precision\b", re.I),
        re.compile(r"\breal\s*\*\s*8\b", re.I),
        re.compile(r"kind\s*\(\s*0\.0d0\s*\)", re.I),
        re.compile(r"(?<![A-Za-z0-9_])(?:\d+(?:\.\d*)?|\.\d+)[dD][+-]?\d+"),
    ]
    forbidden_vendored_names = {
        "r.f90", "r_mod.f90", "blas.f90", "lapack.f90", "arpack.f90"
    }

    for path in FORTRAN:
        rel = path.relative_to(ROOT)
        text = path.read_text(encoding="utf-8")
        if ";" in text:
            for lineno, line in enumerate(text.splitlines(), 1):
                code = line.split("!", 1)[0]
                if ";" in code:
                    fail(f"semicolon in {rel}:{lineno}: {line.strip()}")
        for pattern in forbidden_kind:
            match = pattern.search(text)
            if match:
                fail(f"forbidden real-kind syntax in {rel}: {match.group(0)}")
        if "real(" in text.lower() and path.parent.name in {"src", "test", "example"}:
            if "r_kinds, only : dp" not in text and path.name != "strucchange.f90":
                fail(f"{rel} uses real declarations without importing dp from r_kinds")
        for lineno, line in enumerate(text.splitlines(), 1):
            if len(line) > 132:
                fail(f"line longer than 132 characters in {rel}:{lineno}")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest in hashes:
            fail(f"duplicate Fortran source content: {hashes[digest]} and {rel}")
        hashes[digest] = rel
        if path.name.lower() in forbidden_vendored_names:
            fail(f"forbidden vendored source: {rel}")

    dependency_markers = ["rfortran-core", "rfortran-linalg", "fortran-lapack"]
    for directory in ROOT.iterdir():
        if directory.is_dir() and directory.name.lower() in dependency_markers:
            fail(f"dependency source directory copied into package: {directory.name}")

    artifact_suffixes = {".o", ".obj", ".mod", ".smod", ".exe", ".zip", ".pyc"}
    for path in ROOT.rglob("*"):
        if path.is_file() and path.suffix.lower() in artifact_suffixes:
            fail(f"build/cache/archive artifact present: {path.relative_to(ROOT)}")
        if path.is_dir() and path.name in {"build", "__pycache__", ".pytest_cache"}:
            fail(f"build/cache directory present: {path.relative_to(ROOT)}")


def main() -> None:
    scan_sources()
    fpm = shutil.which("fpm")
    if fpm is None:
        fail("fpm is not available on PATH")
    run([fpm, "build"])
    run([fpm, "test"])
    scan_sources()
    run([fpm, "clean", "--all"])
    scan_sources()
    print("Release checks passed.")


if __name__ == "__main__":
    main()
