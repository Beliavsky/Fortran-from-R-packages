#!/usr/bin/env python3
"""Find byte-identical maintained Fortran source files in the repository."""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
from collections import defaultdict
from pathlib import Path


FORTRAN_SUFFIXES = {".f", ".for", ".f77", ".f90", ".f95", ".f03", ".f08"}
EXCLUDED_PARTS = {".git", "build", "orig", "original", "reference"}


def excluded_path(path: Path) -> bool:
    for part in path.parts:
        name = part.casefold()
        if name in EXCLUDED_PARTS or name.startswith("original_"):
            return True
        if name.endswith("-reference") or name.endswith("_reference"):
            return True
    return False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--minimum-bytes", type=int, default=1)
    parser.add_argument("--top", type=int, default=0, help="show only the largest N groups")
    return parser.parse_args()


def maintained_fortran_files(root: Path, minimum_bytes: int) -> list[Path]:
    tracked = subprocess.run(
        ["git", "-C", str(root), "ls-files", "-z"],
        capture_output=True,
        check=False,
    )
    if tracked.returncode == 0:
        paths = (root / name.decode() for name in tracked.stdout.split(b"\0") if name)
        return [
            path
            for path in paths
            if path.suffix.casefold() in FORTRAN_SUFFIXES
            and not excluded_path(path.relative_to(root))
            and path.is_file()
            and path.stat().st_size >= minimum_bytes
        ]
    files = []
    for directory, subdirectories, names in os.walk(root):
        subdirectories[:] = [
            name for name in subdirectories if not excluded_path(Path(name))
        ]
        for name in names:
            path = Path(directory) / name
            if path.suffix.casefold() in FORTRAN_SUFFIXES and path.stat().st_size >= minimum_bytes:
                files.append(path)
    return files


def tracked_content_ids(root: Path, minimum_bytes: int) -> list[tuple[Path, str]]:
    """Return paths and Git-compatible content IDs, honoring unstaged edits."""
    index = subprocess.run(
        ["git", "-C", str(root), "ls-files", "-s", "-z"],
        capture_output=True,
        check=False,
    )
    if index.returncode != 0:
        return []
    changed_result = subprocess.run(
        ["git", "-C", str(root), "diff", "--name-only", "-z"],
        capture_output=True,
        check=True,
    )
    changed = {os.fsdecode(name) for name in changed_result.stdout.split(b"\0") if name}
    records = []
    for entry in index.stdout.split(b"\0"):
        if not entry:
            continue
        metadata, encoded_name = entry.split(b"\t", 1)
        digest = metadata.split()[1].decode()
        name = os.fsdecode(encoded_name)
        path = root / name
        if not path.is_file() or path.suffix.casefold() not in FORTRAN_SUFFIXES:
            continue
        if excluded_path(Path(name)):
            continue
        size = path.stat().st_size
        if size < minimum_bytes:
            continue
        if name in changed:
            content = path.read_bytes()
            header = f"blob {len(content)}\0".encode()
            digest = hashlib.sha1(header + content).hexdigest()
        records.append((path, digest))
    return records


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    groups: dict[tuple[int, str], list[Path]] = defaultdict(list)
    records = tracked_content_ids(root, args.minimum_bytes)
    if not records:
        records = [
            (path, hashlib.sha256(path.read_bytes()).hexdigest())
            for path in maintained_fortran_files(root, args.minimum_bytes)
        ]
    for path, digest in records:
        groups[(path.stat().st_size, digest)].append(path.relative_to(root))
    duplicates = [(size, paths) for (size, _), paths in groups.items() if len(paths) > 1]
    duplicates.sort(key=lambda item: item[0] * (len(item[1]) - 1), reverse=True)
    if args.top > 0:
        duplicates = duplicates[: args.top]
    for size, paths in duplicates:
        print(f"{len(paths)} copies, {size:,} bytes each, {size*(len(paths)-1):,} duplicate bytes")
        for path in sorted(paths, key=lambda item: str(item).casefold()):
            print(f"  {path}")
    duplicate_bytes = sum(size * (len(paths)-1) for size, paths in duplicates)
    print(f"\n{len(duplicates)} duplicate group(s); {duplicate_bytes:,} duplicate working-tree bytes shown.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
