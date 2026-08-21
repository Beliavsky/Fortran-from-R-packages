#!/usr/bin/env python3
from __future__ import annotations
import hashlib
import pathlib
import sys
import tomllib

root = pathlib.Path(__file__).resolve().parents[1]
with (root / "fpm.toml").open("rb") as handle:
    manifest = tomllib.load(handle)
assert manifest["license"] == "GPL-2.0-only"

for path in sorted(root.glob("src/*.f90")) + sorted(root.glob("test/*.f90")) + sorted(root.glob("app/*.f90")) + sorted(root.glob("example/*.f90")):
    text = path.read_text(encoding="ascii")
    assert "SPDX-License-Identifier: GPL-2.0-only" in text, path
    assert "implicit none" in text.lower(), path
    for number, line in enumerate(text.splitlines(), 1):
        assert len(line) <= 132, f"{path}:{number}: line length {len(line)}"

for name in ["README.md", "COVERAGE.md", "PORTING_NOTES.md", "VALIDATION.md", "NOTICE"]:
    (root / name).read_text(encoding="ascii")

print("source_audit: PASS")
