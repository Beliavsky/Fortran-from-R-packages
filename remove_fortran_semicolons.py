#!/usr/bin/env python3
"""Put semicolon-separated free-form Fortran statements on separate lines.

The scanner recognizes single- and double-quoted character literals, doubled
quote escapes, character literals continued across physical lines, and inline
comments. Semicolons in character literals or comments are never changed.

The default mode reports what would change. Pass ``--write`` to update files or
``--check`` to return a nonzero status when changes would be required.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
import sys
from typing import Iterable


FREE_FORM_SUFFIXES = frozenset({".f90", ".f95", ".f03", ".f08"})
SKIPPED_DIRECTORIES = frozenset(
    {".git", ".fpm", "build", "dependencies", "orig", "original", "upstream"}
)
CASE_RE = re.compile(r"^\s*case(?:\s+default\b|\s*\()", re.IGNORECASE)


class FortranScanError(ValueError):
    """Raised when a character literal remains open at end of file."""


@dataclass(frozen=True)
class TransformResult:
    text: str
    separators: int


def split_eol(line: str) -> tuple[str, str]:
    if line.endswith("\r\n"):
        return line[:-2], "\r\n"
    if line.endswith("\n") or line.endswith("\r"):
        return line[:-1], line[-1]
    return line, ""


def scan_line(line: str, quote: str | None) -> tuple[list[int], int | None, str | None]:
    """Find statement separators and the comment start on one physical line."""
    separators: list[int] = []
    comment_start: int | None = None
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
        if ch in {"'", '"'}:
            quote = ch
        elif ch == "!":
            comment_start = i
            break
        elif ch == ";":
            separators.append(i)
        i += 1
    return separators, comment_start, quote


def split_statement_line(line: str, separators: list[int], comment_start: int | None) -> list[str]:
    """Split one physical line at already-validated separator positions."""
    code_end = len(line) if comment_start is None else comment_start
    comment = "" if comment_start is None else line[comment_start:]
    boundaries = [position for position in separators if position < code_end]
    starts = [0, *(position + 1 for position in boundaries)]
    ends = [*boundaries, code_end]
    indent = re.match(r"[ \t]*", line).group(0)
    output: list[str] = []

    for index, (start, end) in enumerate(zip(starts, ends)):
        segment = line[start:end]
        if index == 0:
            statement = segment.rstrip()
        else:
            statement = indent + segment.strip()
        if statement.strip():
            output.append(statement)

    if comment:
        if output:
            output[-1] = output[-1].rstrip() + " " + comment.lstrip()
        else:
            output.append(indent + comment.lstrip())
    if not output:
        output.append(indent.rstrip())
    return output


def transform_text(text: str, keep_simple_case: bool = False) -> TransformResult:
    """Return transformed free-form Fortran source and the split count."""
    lines = text.splitlines(keepends=True)
    if not lines and text:
        lines = [text]
    default_eol = "\r\n" if "\r\n" in text else "\n"
    output: list[str] = []
    quote: str | None = None
    split_count = 0

    for raw in lines:
        body, eol = split_eol(raw)
        if quote is None and body.lstrip().startswith("#"):
            output.append(raw)
            continue
        separators, comment_start, quote = scan_line(body, quote)
        if keep_simple_case and len(separators) == 1 and CASE_RE.match(body):
            output.append(raw)
            continue
        if not separators:
            output.append(raw)
            continue

        replacements = split_statement_line(body, separators, comment_start)
        split_count += len(separators)
        for index, replacement in enumerate(replacements):
            if eol:
                replacement_eol = eol
            elif index < len(replacements) - 1:
                replacement_eol = default_eol
            else:
                replacement_eol = ""
            output.append(replacement + replacement_eol)

    if quote is not None:
        raise FortranScanError(f"unterminated {quote} character literal")
    return TransformResult("".join(output), split_count)


def is_skipped(path: Path, root: Path) -> bool:
    try:
        relative = path.relative_to(root)
    except ValueError:
        relative = path
    return any(part.lower() in SKIPPED_DIRECTORIES for part in relative.parts[:-1])


def iter_source_files(paths: Iterable[Path], include_upstream: bool) -> list[Path]:
    files: set[Path] = set()
    for supplied in paths:
        path = supplied.resolve()
        if path.is_file():
            if path.suffix.lower() in FREE_FORM_SUFFIXES:
                files.add(path)
            continue
        if not path.is_dir():
            raise FileNotFoundError(f"path does not exist: {supplied}")
        for candidate in path.rglob("*"):
            if not candidate.is_file() or candidate.suffix.lower() not in FREE_FORM_SUFFIXES:
                continue
            if not include_upstream and is_skipped(candidate, path):
                continue
            files.add(candidate.resolve())
    return sorted(files, key=lambda item: str(item).lower())


def decode_source(data: bytes) -> tuple[str, bool]:
    has_bom = data.startswith(b"\xef\xbb\xbf")
    payload = data[3:] if has_bom else data
    return payload.decode("utf-8", errors="surrogateescape"), has_bom


def encode_source(text: str, has_bom: bool) -> bytes:
    payload = text.encode("utf-8", errors="surrogateescape")
    return (b"\xef\xbb\xbf" if has_bom else b"") + payload


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", type=Path, default=[Path(".")])
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true", help="update affected files in place")
    mode.add_argument("--check", action="store_true", help="exit 1 if any file would change")
    parser.add_argument(
        "--keep-simple-case",
        action="store_true",
        help="retain lines of the form 'case (...); one_statement'",
    )
    parser.add_argument(
        "--include-upstream",
        action="store_true",
        help="also scan orig, original, upstream, and dependencies directories",
    )
    parser.add_argument("--verbose", action="store_true", help="list affected files")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        files = iter_source_files(args.paths, args.include_upstream)
    except FileNotFoundError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    updates: list[tuple[Path, bytes, int]] = []
    errors: list[str] = []
    for path in files:
        original = path.read_bytes()
        text, has_bom = decode_source(original)
        try:
            result = transform_text(text, keep_simple_case=args.keep_simple_case)
        except FortranScanError as error:
            errors.append(f"{path}: {error}")
            continue
        transformed = encode_source(result.text, has_bom)
        if transformed != original:
            updates.append((path, transformed, result.separators))

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        print("No files were changed because scanning failed.", file=sys.stderr)
        return 2

    if args.verbose:
        for path, _data, separators in updates:
            print(f"{path}: {separators}")
    if args.write:
        for path, transformed, _separators in updates:
            path.write_bytes(transformed)

    action = "Updated" if args.write else "Would update"
    separator_count = sum(item[2] for item in updates)
    print(f"Scanned {len(files)} free-form Fortran file(s).")
    print(f"{action} {len(updates)} file(s) and split {separator_count} separator(s).")
    if updates and not args.write and not args.check:
        print("Run again with --write to apply the changes.")
    return 1 if args.check and updates else 0


if __name__ == "__main__":
    raise SystemExit(main())
