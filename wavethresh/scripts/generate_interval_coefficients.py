#!/usr/bin/env python3
"""Regenerate the interval-wavelet coefficient tables from upstream C."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "upstream" / "src" / "functions.c"
TARGET = ROOT / "src" / "wavethresh_interval_coefficients.f90"
NAMES = ("Interior", "Left", "Right", "LeftPre", "RightPre")


def values_for(text: str, name: str) -> list[str]:
    """Return normalized real literals from one upstream C table."""
    match = re.search(rf"const double {name}\[\] = \{{(.*?)\}};", text, re.S)
    if match is None:
        raise RuntimeError(f"could not find {name}")
    values = re.findall(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?", match.group(1))
    return [value.replace("E", "e").replace("e", "e") + "_dp" for value in values]


def declaration(name: str, values: list[str]) -> list[str]:
    """Format one coefficient table as a readable Fortran parameter array."""
    identifier = re.sub(r"(?<!^)(?=[A-Z])", "_", name).lower()
    lines = [f"   real(dp), parameter :: {identifier}({len(values)}) = [ &"]
    width = 4
    for start in range(0, len(values), width):
        group = ", ".join(values[start : start + width])
        suffix = ", &" if start + width < len(values) else " ]"
        lines.append(f"      {group}{suffix}")
    return lines


def main() -> None:
    """Write the generated Fortran module deterministically."""
    text = SOURCE.read_text(encoding="utf-8", errors="replace")
    output = [
        "! SPDX-License-Identifier: GPL-2.0-or-later",
        "! Generated from the interval-wavelet tables in upstream/src/functions.c.",
        "module wavethresh_interval_coefficients",
        "   use r_kinds, only : dp",
        "   implicit none",
        "   private",
        "",
        "   public :: interior, left, right, left_pre, right_pre",
        "",
    ]
    for name in NAMES:
        output.extend(declaration(name, values_for(text, name)))
        output.append("")
    output.append("end module wavethresh_interval_coefficients")
    TARGET.write_text("\n".join(output) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
