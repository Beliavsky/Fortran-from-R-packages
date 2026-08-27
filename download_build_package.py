#!/usr/bin/env python3
"""Sparse-clone one translated package, then build and run it with FPM."""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path
from typing import Iterator


REPOSITORY_URL = "https://github.com/Beliavsky/Fortran-from-R-packages.git"
VALID_PACKAGE_NAME = re.compile(r"[A-Za-z0-9._-]+\Z")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "package",
        help="top-level package directory to download, for example rugarch",
    )
    return parser.parse_args()


def require_program(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(f"{name} was not found on PATH")


def run(command: list[str], cwd: Path | None = None) -> None:
    location = f" in {cwd}" if cwd is not None else ""
    print(f"Running: {subprocess.list2cmdline(command)}{location}", flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def local_dependency_paths(value: object) -> Iterator[str]:
    """Yield local FPM dependency paths from nested manifest tables."""
    if not isinstance(value, dict):
        return
    for key, child in value.items():
        if key in {"dependencies", "dev-dependencies"} and isinstance(child, dict):
            for specification in child.values():
                if isinstance(specification, dict):
                    path = specification.get("path")
                    if isinstance(path, str):
                        yield path
        yield from local_dependency_paths(child)


def checkout_with_dependencies(clone: Path, package: str) -> Path:
    """Check out a package and recursively include sibling path dependencies."""
    selected = {package}
    run(["git", "-C", str(clone), "sparse-checkout", "set", package])
    package_directory = clone / package
    if not (package_directory / "fpm.toml").is_file():
        raise RuntimeError(
            f'package "{package}" was not found in the repository'
        )
    pending = [package_directory]
    visited: set[Path] = set()
    clone_root = clone.resolve()

    while pending:
        project = pending.pop()
        project = project.resolve()
        if project in visited:
            continue
        visited.add(project)
        manifest = project / "fpm.toml"
        if not manifest.is_file():
            raise RuntimeError(f'expected local dependency manifest "{manifest}"')
        try:
            with manifest.open("rb") as stream:
                data = tomllib.load(stream)
        except (OSError, tomllib.TOMLDecodeError) as exc:
            raise RuntimeError(f'could not read manifest "{manifest}": {exc}') from exc

        for dependency_path in local_dependency_paths(data):
            dependency = (project / dependency_path).resolve()
            try:
                relative = dependency.relative_to(clone_root)
            except ValueError as exc:
                raise RuntimeError(
                    f'local dependency path escapes the clone: "{dependency_path}"'
                ) from exc
            if not relative.parts:
                raise RuntimeError(f'invalid local dependency path: "{dependency_path}"')
            top_directory = relative.parts[0]
            if top_directory not in selected:
                selected.add(top_directory)
                run(
                    [
                        "git",
                        "-C",
                        str(clone),
                        "sparse-checkout",
                        "set",
                        *sorted(selected),
                    ]
                )
            pending.append(dependency)

    if len(selected) > 1:
        dependencies = ", ".join(sorted(selected - {package}))
        print(f"Included local dependencies: {dependencies}")
    return package_directory


def main() -> int:
    args = parse_arguments()
    package = args.package
    if VALID_PACKAGE_NAME.fullmatch(package) is None:
        print(
            "ERROR: the package name may contain only letters, digits, '.', '_', and '-'",
            file=sys.stderr,
        )
        return 2

    destination = Path.cwd() / f"{package}-download"
    if destination.exists():
        print(f'ERROR: destination already exists: "{destination}"', file=sys.stderr)
        return 2

    try:
        require_program("git")
        require_program("fpm")
        run(
            [
                "git",
                "clone",
                "--depth",
                "1",
                "--filter=blob:none",
                "--sparse",
                REPOSITORY_URL,
                str(destination),
            ]
        )
        package_directory = checkout_with_dependencies(destination, package)
        manifest = package_directory / "fpm.toml"
        if not manifest.is_file():
            raise RuntimeError(
                f'package "{package}" was not found; expected manifest "{manifest}"'
            )

        run(["fpm", "build"], cwd=package_directory)
        run(["fpm", "run"], cwd=package_directory)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except subprocess.CalledProcessError as exc:
        print(
            f"ERROR: command failed with exit status {exc.returncode}",
            file=sys.stderr,
        )
        return exc.returncode or 1

    print()
    print(f'Package "{package}" was downloaded, built, and run successfully.')
    print(f'Location: "{package_directory}"')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
