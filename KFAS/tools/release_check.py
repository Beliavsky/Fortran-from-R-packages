#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Run the release checks required before packaging KFAS."""
from __future__ import annotations

import hashlib
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXCLUDED_DIRECTORIES = {"build", ".git", "__pycache__", ".cache"}


def maintained_files(pattern: str) -> list[Path]:
    return [
        path
        for path in ROOT.rglob(pattern)
        if not any(part in EXCLUDED_DIRECTORIES for part in path.relative_to(ROOT).parts)
    ]


def run(*args: str) -> None:
    print("+", " ".join(args))
    subprocess.run(args, cwd=ROOT, check=True)


def fail(message: str) -> None:
    raise SystemExit(message)


def check_sources() -> None:
    fixed_form = []
    for suffix in ("*.f", "*.for", "*.f77"):
        fixed_form.extend(maintained_files(suffix))
    if fixed_form:
        fail("fixed-form Fortran source found: " + repr(fixed_form))

    fortran = sorted(maintained_files("*.f90"))
    if not fortran:
        fail("no Fortran source found")

    hashes: dict[str, list[Path]] = defaultdict(list)
    for path in fortran:
        data = path.read_bytes()
        hashes[hashlib.sha256(data).hexdigest()].append(path)
        for lineno, line in enumerate(data.decode("utf-8").splitlines(), 1):
            if ";" in line:
                fail(f"semicolon found in {path.relative_to(ROOT)}:{lineno}")
            if len(line) > 132:
                fail(f"line longer than 132 columns in {path.relative_to(ROOT)}:{lineno}")

    duplicates = [paths for paths in hashes.values() if len(paths) > 1]
    if duplicates:
        fail("duplicate Fortran source content: " + repr(duplicates))

    disallowed_names = {"r.f90", "r_mod.f90"}
    for path in maintained_files("*"):
        if path.is_file() and path.name.lower() in disallowed_names:
            fail(f"disallowed copied helper: {path.relative_to(ROOT)}")
        if path.is_file() and path.suffix.lower() in {".f90", ".f", ".for", ".f77"}:
            low = path.name.lower()
            if any(name in low for name in ("blas", "lapack", "arpack")):
                fail(f"disallowed copied numerical dependency source: {path.relative_to(ROOT)}")
        if path.is_dir() and (path.name.lower().startswith("rfortran-") or path.name.lower() in {"vendor", "third_party"}):
            fail(f"vendored dependency directory: {path.relative_to(ROOT)}")


def check_clean_tree() -> None:
    bad_suffixes = {".o", ".obj", ".mod", ".smod", ".a", ".so", ".dll", ".dylib", ".exe", ".zip"}
    bad_dirs = {"build", "__pycache__", ".cache"}
    for path in ROOT.rglob("*"):
        if path.is_file() and path.suffix.lower() in bad_suffixes:
            fail(f"build/package artifact remains: {path.relative_to(ROOT)}")
        if path.is_dir() and path.name in bad_dirs:
            fail(f"build/cache directory remains: {path.relative_to(ROOT)}")


def main() -> int:
    if shutil.which("fpm") is None:
        fail("fpm is required for the release check but was not found on PATH")
    run("fpm", "build")
    run("fpm", "test")
    check_sources()
    run("fpm", "clean", "--all")
    check_clean_tree()
    print("KFAS release checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
