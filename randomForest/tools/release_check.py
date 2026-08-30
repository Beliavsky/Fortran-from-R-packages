#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAINTAINED = [ROOT / "src", ROOT / "test", ROOT / "example"]
FORBIDDEN_SOURCE_NAMES = {
    "r.f90", "r_mod.f90", "r_kinds.f90", "r_quantiles.f90", "r_robust.f90",
    "r_linalg.f90", "blas.f90", "lapack.f90", "arpack.f90",
}
ARTIFACT_SUFFIXES = {".o", ".obj", ".mod", ".smod", ".a", ".so", ".dll", ".dylib", ".exe", ".zip"}


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, cwd=ROOT, check=True)


def fail(message: str) -> None:
    print("release check failed:", message, file=sys.stderr)
    raise SystemExit(1)


def fortran_files() -> list[Path]:
    files: list[Path] = []
    for directory in MAINTAINED:
        files.extend(sorted(directory.glob("*.f90")))
    return files


def scan_sources() -> None:
    files = fortran_files()
    if not files:
        fail("no maintained Fortran sources found")

    hashes: dict[str, Path] = {}
    max_line = 0
    for path in files:
        text = path.read_text(encoding="utf-8")
        digest = hashlib.sha256(text.encode()).hexdigest()
        if digest in hashes:
            fail(f"duplicate Fortran sources: {hashes[digest]} and {path}")
        hashes[digest] = path
        for lineno, line in enumerate(text.splitlines(), start=1):
            max_line = max(max_line, len(line))
            code = line.split("!", 1)[0]
            if ";" in code:
                fail(f"semicolon-separated statement in {path}:{lineno}")
        lower = text.lower()
        if re.search(r"\bdouble\s+precision\b|\breal\s*\*\s*8\b|kind\s*\(\s*0\.0d0\s*\)", lower):
            fail(f"forbidden real-kind form in {path}")
        if re.search(r"(?i)(?:\d(?:\.\d*)?|\.\d+)d[+-]?\d+", text):
            fail(f"D exponent literal in {path}; use _dp literals")
        if re.search(r"\breal\s*(?:,|::)", lower):
            fail(f"default-real declaration in {path}")
        if "real64" in lower or "iso_fortran_env" in lower:
            fail(f"maintained source defines/imports another real kind in {path}")

    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(ROOT)
        if path.name.lower() in FORBIDDEN_SOURCE_NAMES:
            fail(f"copied shared/dependency source found: {rel}")
        if path.suffix.lower() in ARTIFACT_SUFFIXES:
            fail(f"build/archive artifact found: {rel}")
        if any(part.lower() in {"build", "__pycache__", ".pytest_cache", ".cache"} for part in rel.parts):
            fail(f"cache/build directory found: {rel}")

    print(f"maintained Fortran files: {len(files)}")
    print(f"unique source hashes: {len(hashes)}")
    print(f"maximum maintained Fortran line length: {max_line}")
    if max_line > 132:
        fail("maintained source contains a line longer than 132 characters")


def main() -> None:
    if shutil.which("fpm") is None:
        fail("fpm not found")
    run(["fpm", "build"])
    run(["fpm", "test"])
    scan_sources()
    run(["fpm", "clean", "--all"])


if __name__ == "__main__":
    main()
