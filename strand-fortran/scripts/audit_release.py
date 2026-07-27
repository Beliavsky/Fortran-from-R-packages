#!/usr/bin/env python3
"""Audit the strand-fortran release tree."""
from __future__ import annotations

import hashlib
from pathlib import Path
import tomllib

ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check_manifest(path: Path) -> None:
    if not path.exists():
        return
    for line in path.read_text(encoding="ascii").splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        target = ROOT / relative
        if not target.is_file():
            raise SystemExit(f"missing checksum target: {relative}")
        actual = sha256(target)
        if actual != expected:
            raise SystemExit(f"checksum mismatch: {relative}")


with (ROOT / "fpm.toml").open("rb") as stream:
    manifest = tomllib.load(stream)
if manifest.get("license") != "GPL-3.0-only":
    raise SystemExit("unexpected manifest license")

for path in [*ROOT.glob("src/*.f90"), *ROOT.glob("test/*.f90"),
             *ROOT.glob("app/*.f90"), *ROOT.glob("example/*.f90")]:
    text = path.read_text(encoding="ascii")
    lines = text.splitlines()
    if "SPDX-License-Identifier: GPL-3.0-only" not in "\n".join(lines[:5]):
        raise SystemExit(f"missing SPDX identifier: {path.relative_to(ROOT)}")
    if not any("implicit none" in line.lower() for line in lines):
        raise SystemExit(f"missing implicit none: {path.relative_to(ROOT)}")
    for number, line in enumerate(lines, 1):
        if len(line) > 132:
            raise SystemExit(
                f"line exceeds 132 columns: {path.relative_to(ROOT)}:{number}"
            )

for pattern in ("*.o", "*.mod", "*.smod", "*.exe"):
    for path in ROOT.rglob(pattern):
        if "build" not in path.parts:
            raise SystemExit(f"compiled artifact in release tree: {path.relative_to(ROOT)}")

check_manifest(ROOT / "provenance" / "SOURCE_ARCHIVE.sha256")
check_manifest(ROOT / "provenance" / "ORIGINAL_FILES.sha256")
check_manifest(ROOT / "provenance" / "TRANSLATED_FILES.sha256")
print("archive_audit: PASS")
