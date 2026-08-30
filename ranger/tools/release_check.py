#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import pathlib
import re
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]


def run(*args: str) -> None:
    subprocess.run(args, cwd=ROOT, check=True)


def fail(message: str) -> None:
    raise SystemExit(message)


def code_without_comments_or_strings(line: str) -> str:
    out: list[str] = []
    quote: str | None = None
    i = 0
    while i < len(line):
        ch = line[i]
        if quote is not None:
            if ch == quote:
                if i + 1 < len(line) and line[i + 1] == quote:
                    i += 2
                    continue
                quote = None
            i += 1
            continue
        if ch in ("'", '"'):
            quote = ch
            i += 1
            continue
        if ch == "!":
            break
        out.append(ch)
        i += 1
    return "".join(out)


def main() -> None:
    if shutil.which("fpm") is None:
        fail("fpm is required for the release gate")
    run("fpm", "build")
    run("fpm", "test")

    maintained = sorted(list((ROOT / "src").glob("*.f90")) + list((ROOT / "test").glob("*.f90")) +
                        list((ROOT / "example").glob("*.f90")))
    if not maintained:
        fail("no maintained Fortran files found")

    hashes: dict[str, pathlib.Path] = {}
    forbidden = [
        re.compile(r"\bdouble\s+precision\b", re.I),
        re.compile(r"\breal\s*\*\s*8\b", re.I),
        re.compile(r"kind\s*\(\s*0\.0d0\s*\)", re.I),
        re.compile(r"(?:\d|\.)[dD][+-]?\d"),
        re.compile(r"\breal64\b", re.I),
        re.compile(r"\biso_fortran_env\b", re.I),
    ]
    for path in maintained:
        data = path.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        if digest in hashes:
            fail(f"duplicate maintained source: {path} and {hashes[digest]}")
        hashes[digest] = path
        text = data.decode("utf-8")
        if "use r_kinds, only : dp" not in text and "use r_kinds, only : dp," not in text:
            fail(f"maintained Fortran file does not import shared dp: {path}")
        real_literal = re.compile(
            r"(?<![A-Za-z0-9_])(?:\d+\.\d*|\.\d+)(?:[eE][+-]?\d+)?(?!_dp)(?![A-Za-z0-9_])"
            r"|(?<![A-Za-z0-9_])\d+[eE][+-]?\d+(?!_dp)(?![A-Za-z0-9_])"
        )
        for lineno, line in enumerate(text.splitlines(), 1):
            if len(line) > 132:
                fail(f"line longer than 132 characters: {path}:{lineno}")
            code = code_without_comments_or_strings(line)
            if ";" in code:
                fail(f"semicolon found in maintained Fortran code: {path}:{lineno}")
            for pattern in forbidden:
                if pattern.search(code):
                    fail(f"forbidden real-kind form in {path}:{lineno}: {pattern.pattern}")
            match = real_literal.search(code)
            if match:
                fail(f"real literal without _dp suffix: {path}:{lineno}: {match.group(0)}")

    forbidden_names = {
        "r.f90", "r_mod.f90", "r_kinds.f90", "r_linalg.f90",
        "blas.f90", "lapack.f90", "arpack.f90"
    }
    for path in (ROOT / "src").glob("*.f90"):
        if not path.name.lower().startswith("ranger"):
            fail(f"unexpected non-ranger source in src/: {path.name}")

    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(ROOT)
        lower = path.name.lower()
        if lower in forbidden_names:
            fail(f"vendored shared source found: {rel}")
        if lower.endswith((".o", ".obj", ".mod", ".smod", ".exe", ".a", ".so", ".dll", ".zip")):
            fail(f"build/archive artifact found: {rel}")
        if "build" in rel.parts or "__pycache__" in rel.parts or ".pytest_cache" in rel.parts:
            fail(f"cache/build directory found: {rel}")

    run("fpm", "clean", "--all")
    print(f"release checks passed for {len(maintained)} maintained Fortran files")


if __name__ == "__main__":
    main()
