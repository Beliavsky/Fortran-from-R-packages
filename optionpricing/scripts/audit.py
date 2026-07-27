#!/usr/bin/env python3
"""Audit the release tree for manifest, source, and provenance requirements."""
from __future__ import annotations
import hashlib
from pathlib import Path
import tomllib

ROOT = Path(__file__).resolve().parents[1]

with (ROOT / "fpm.toml").open("rb") as fh:
    manifest = tomllib.load(fh)
assert manifest["name"] == "optionpricing"
assert manifest["license"] == "GPL-2.0-only OR GPL-3.0-only"
print("manifest_audit: PASS")

text_suffixes = {".f90", ".md", ".toml", ".py", ".sh", ".txt", ""}
for path in ROOT.rglob("*"):
    if not path.is_file() or "original" in path.parts:
        continue
    if path.name.endswith(".zip") or path.name.endswith(".sha256"):
        continue
    if path.suffix.lower() in text_suffixes or path.name in {"NOTICE"}:
        data = path.read_bytes()
        data.decode("ascii")

for path in list((ROOT / "src").glob("*.f90")) + \
            list((ROOT / "test").glob("*.f90")) + \
            list((ROOT / "app").glob("*.f90")) + \
            list((ROOT / "example").glob("*.f90")):
    lines = path.read_text(encoding="ascii").splitlines()
    assert lines and "SPDX-License-Identifier:" in lines[0]
    assert any("implicit none" in line.lower() for line in lines)
    assert max(map(len, lines), default=0) <= 132
print("source_audit: PASS")

for forbidden in ["build", ".git", "__pycache__"]:
    assert not (ROOT / forbidden).exists()
for path in ROOT.rglob("*"):
    assert path.suffix.lower() not in {".o", ".mod", ".a", ".exe"}
print("release_tree_audit: PASS")

def verify_manifest(manifest_path: Path) -> None:
    for line in manifest_path.read_text(encoding="ascii").splitlines():
        if not line.strip():
            continue
        digest, rel = line.split("  ", 1)
        data = (ROOT / rel).read_bytes()
        actual = hashlib.sha256(data).hexdigest()
        assert actual == digest, rel

verify_manifest(ROOT / "provenance" / "SOURCE_ARCHIVE.sha256")
verify_manifest(ROOT / "provenance" / "ORIGINAL_FILES.sha256")
verify_manifest(ROOT / "provenance" / "TRANSLATED_FILES.sha256")
print("checksum_audit: PASS")
