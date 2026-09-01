#!/usr/bin/env python3
"""Refactor FPM packages one at a time and reject regressions.

The package order is dependency-first.  For each package, the driver first
establishes a passing FPM baseline, snapshots the maintained free-form Fortran
sources, applies mechanical transformations, and repeats the same FPM checks.
If a post-refactoring check fails, every changed source is restored byte for
byte and a log is retained.

Downloaded upstream sources, vendored dependencies, and fixed-form Fortran are
intentionally excluded.  Semantic documentation is outside the scope of this
mechanical workflow.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10
    tomllib = None


FREE_FORM_SUFFIXES = {".f90", ".f95", ".f03", ".f08"}
EXCLUDED_PARTS = {
    ".git",
    ".fpm",
    "build",
    "dependencies",
    "orig",
    "original",
    "upstream",
    "vendor",
}
DEFAULT_XPURE = Path(r"C:\python\fortran\xpure.py")
DEFAULT_XPARAM = Path(r"C:\python\fortran\xparam.py")
NUMERIC_STATEMENT_LABEL_RE = re.compile(r"^\s*\d{1,5}\s+\S")
MAX_FORMATTABLE_CONTINUED_STATEMENT = 1000


@dataclass(frozen=True)
class Package:
    path: Path
    dependencies: tuple[str, ...]
    local_dependencies: tuple[str, ...]

    @property
    def name(self) -> str:
        return self.path.name


@dataclass
class CommandResult:
    command: list[str]
    returncode: int | None
    output: str
    seconds: float


def load_manifest(path: Path) -> dict:
    """Read an FPM manifest, returning an empty mapping on parse failure."""
    if tomllib is None:
        raise RuntimeError("Python 3.11 or newer is required for TOML parsing")
    try:
        with path.open("rb") as stream:
            return tomllib.load(stream)
    except (OSError, tomllib.TOMLDecodeError):
        return {}


def discover_packages(root: Path) -> dict[str, Package]:
    """Discover top-level FPM packages and their top-level path dependencies."""
    package_paths = sorted(
        (p.resolve() for p in root.iterdir() if p.is_dir() and (p / "fpm.toml").is_file()),
        key=lambda p: p.name.casefold(),
    )
    by_resolved_path = {path: path.name for path in package_paths}
    packages: dict[str, Package] = {}
    for path in package_paths:
        manifest = load_manifest(path / "fpm.toml")
        dependencies = manifest.get("dependencies", {})
        if not isinstance(dependencies, dict):
            dependencies = {}
        local: set[str] = set()
        for specification in dependencies.values():
            if not isinstance(specification, dict) or "path" not in specification:
                continue
            target = (path / str(specification["path"])).resolve()
            name = by_resolved_path.get(target)
            if name is not None and name != path.name:
                local.add(name)
        dependency_names = tuple(sorted((str(name) for name in dependencies), key=str.casefold))
        packages[path.name.casefold()] = Package(
            path,
            dependency_names,
            tuple(sorted(local, key=str.casefold)),
        )
    return packages


def dependency_order(packages: dict[str, Package]) -> list[Package]:
    """Return a topological order, prioritizing foundations with many users."""
    prerequisites = {
        key: {name.casefold() for name in package.local_dependencies if name.casefold() in packages}
        for key, package in packages.items()
    }
    dependents: dict[str, set[str]] = {key: set() for key in packages}
    for key, dependencies in prerequisites.items():
        for dependency in dependencies:
            dependents[dependency].add(key)

    def priority(key: str) -> tuple[int, int, int, str]:
        dependency_count = len(packages[key].dependencies)
        return (dependency_count != 0, dependency_count, -len(dependents[key]), packages[key].name.casefold())

    ready = sorted((key for key, dependencies in prerequisites.items() if not dependencies), key=priority)
    ordered: list[Package] = []
    emitted: set[str] = set()
    while ready:
        next_ready: set[str] = set()
        for key in ready:
            if key in emitted:
                continue
            emitted.add(key)
            ordered.append(packages[key])
            for dependent in dependents[key]:
                prerequisites[dependent].discard(key)
                if not prerequisites[dependent]:
                    next_ready.add(dependent)
        ready = sorted(next_ready - emitted, key=priority)

    # Cycles are not expected, but list any cycle members deterministically.
    ordered.extend(packages[key] for key in sorted(set(packages) - emitted, key=priority))
    return ordered


def maintained_sources(package: Path) -> list[Path]:
    """Return maintained free-form Fortran sources, excluding imported trees."""
    sources = []
    for path in package.rglob("*"):
        if not path.is_file() or path.suffix.casefold() not in FREE_FORM_SUFFIXES:
            continue
        relative_parts = {part.casefold() for part in path.relative_to(package).parts[:-1]}
        if relative_parts & EXCLUDED_PARTS:
            continue
        sources.append(path)
    return sorted(sources, key=lambda path: str(path.relative_to(package)).casefold())


def package_has_tests(package: Path) -> bool:
    """Return whether FPM should have tests to run for this package."""
    manifest = load_manifest(package / "fpm.toml")
    if manifest.get("test"):
        return True
    if manifest.get("build", {}).get("auto-tests", True) is False:
        return False
    test_dir = package / "test"
    return test_dir.is_dir() and any(path in maintained_sources(package) for path in test_dir.rglob("*"))


def stop_process_tree(process: subprocess.Popen[str]) -> None:
    """Terminate a timed-out command and its children."""
    if os.name == "nt":
        subprocess.run(
            ["taskkill", "/PID", str(process.pid), "/T", "/F"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    else:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        process.kill()


def run_command(command: list[str], cwd: Path, timeout: float) -> CommandResult:
    """Run a command with combined captured output and a timeout."""
    started = time.perf_counter()
    creationflags = subprocess.CREATE_NEW_PROCESS_GROUP if os.name == "nt" else 0
    process = subprocess.Popen(
        command,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        creationflags=creationflags,
        start_new_session=os.name != "nt",
    )
    try:
        output, _ = process.communicate(timeout=timeout)
        return CommandResult(command, process.returncode, output, time.perf_counter() - started)
    except subprocess.TimeoutExpired:
        stop_process_tree(process)
        output, _ = process.communicate()
        return CommandResult(command, None, output, time.perf_counter() - started)


def fpm_checks(package: Path, timeout: float) -> list[CommandResult]:
    """Build examples and library, then run tests when the package has tests."""
    results = [run_command(["fpm", "build"], package, timeout)]
    if results[-1].returncode == 0 and package_has_tests(package):
        results.append(run_command(["fpm", "test"], package, timeout))
    return results


def checks_passed(results: list[CommandResult]) -> bool:
    return bool(results) and all(result.returncode == 0 for result in results)


def has_numeric_statement_labels(source: Path) -> bool:
    """Return whether a source uses legacy numeric statement labels.

    Fprettify 0.3.7 can silently remove a label from statements such as
    ``60 CALL TRSBOX(...)`` while preserving branches to that label. Skip the
    formatter until those control-flow labels have been modernized safely.
    """
    try:
        lines = source.read_text(encoding="utf-8", errors="ignore").splitlines()
    except OSError:
        return False
    return any(
        NUMERIC_STATEMENT_LABEL_RE.match(line) and not line.lstrip().startswith("!")
        for line in lines
    )


def has_oversized_continued_statement(source: Path) -> bool:
    """Return whether a continued statement is too large for safe fprettify use.

    Fprettify 0.3.7 can collapse very large continued array constructors onto
    one physical line instead of reflowing them. Such generated data tables do
    not benefit from formatting and may then exceed the compiler line limit.
    """
    try:
        lines = source.read_text(encoding="utf-8", errors="ignore").splitlines()
    except OSError:
        return False
    continued_length = 0
    in_continuation = False
    for line in lines:
        code = line.split("!", 1)[0].rstrip()
        continues = code.endswith("&")
        if in_continuation or continues:
            continued_length += len(code.strip().lstrip("&").rstrip("&"))
            if continued_length > MAX_FORMATTABLE_CONTINUED_STATEMENT:
                return True
        if not continues:
            continued_length = 0
        in_continuation = continues
    return False


def format_sources(sources: list[Path], package: Path, timeout: float) -> list[CommandResult]:
    """Format each file separately to avoid Windows command-line limits."""
    results = []
    for source in sources:
        if has_numeric_statement_labels(source) or has_oversized_continued_statement(source):
            continue
        result = run_command(
            [
                "fprettify",
                "--indent",
                "3",
                "--line-length",
                "120",
                "--whitespace",
                "3",
                "--enable-decl",
                "--strict-indent",
                str(source),
            ],
            package,
            timeout,
        )
        results.append(result)
        if result.returncode != 0:
            break
    return results


def xpure_command(xpure: Path, sources: list[Path], elemental: bool) -> list[str]:
    command = [
        sys.executable,
        str(xpure),
        "--fix",
        "--all-candidates",
        "--iterate",
        "--no-backup",
        "--strict-unknown-calls",
    ]
    if elemental:
        command.append("--suggest-elemental")
    command.extend(str(path) for path in sources)
    return command


def xparam_command(xparam: Path, sources: list[Path]) -> list[str]:
    """Build the conservative initialized-constant promotion command."""
    return [
        sys.executable,
        str(xparam),
        "--fix",
        "--initialized-only",
        "--no-backup",
        *[str(path) for path in sources],
    ]


def write_log(log_path: Path, heading: str, results: list[CommandResult]) -> None:
    """Write enough command output to diagnose a baseline failure or regression."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    sections = [heading]
    for result in results:
        status = "TIMEOUT" if result.returncode is None else str(result.returncode)
        sections.extend(
            [
                "",
                f"> {subprocess.list2cmdline(result.command)}",
                f"exit={status}; seconds={result.seconds:.1f}",
                result.output.rstrip(),
            ]
        )
    log_path.write_text("\n".join(sections).rstrip() + "\n", encoding="utf-8")


def restore_sources(snapshot: Path, package: Path, sources: list[Path]) -> None:
    """Restore the exact source bytes captured before transformation."""
    for source in sources:
        saved = snapshot / source.relative_to(package)
        source.write_bytes(saved.read_bytes())


def changed_sources(snapshot: Path, package: Path, sources: list[Path]) -> list[Path]:
    return [source for source in sources if source.read_bytes() != (snapshot / source.relative_to(package)).read_bytes()]


def transform_package(
    package: Package,
    root: Path,
    xpure: Path,
    xparam: Path,
    timeout: float,
    clean: bool,
) -> bool:
    """Transform one package transactionally; return True on success."""
    sources = maintained_sources(package.path)
    log_path = root / "refactor_logs" / f"{package.name}.log"
    if not sources:
        print(f"SKIP {package.name}: no maintained free-form sources")
        return True

    print(f"BASELINE {package.name}: {len(sources)} source file(s)", flush=True)
    baseline = fpm_checks(package.path, timeout)
    if not checks_passed(baseline):
        write_log(log_path, f"Baseline failure: {package.name}", baseline)
        print(f"BASELINE FAIL {package.name}; unchanged; log: {log_path}")
        return False

    with tempfile.TemporaryDirectory(prefix=f"r_fpm_{package.name}_") as temporary:
        snapshot = Path(temporary)
        for source in sources:
            target = snapshot / source.relative_to(package.path)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(source.read_bytes())

        transformations: list[CommandResult] = []
        transformations.append(
            run_command(
                [
                    sys.executable,
                    str(root / "remove_fortran_semicolons.py"),
                    "--keep-simple-case",
                    "--write",
                    *[str(path) for path in sources],
                ],
                package.path,
                timeout,
            )
        )
        if checks_passed(transformations):
            transformations.append(
                run_command(
                    [
                        sys.executable,
                        str(root / "modernize_fortran_syntax.py"),
                        "--write",
                        *[str(path) for path in sources],
                    ],
                    package.path,
                    timeout,
                )
            )
        if checks_passed(transformations):
            transformations.append(
                run_command(xparam_command(xparam, sources), package.path, timeout)
            )
        if checks_passed(transformations):
            transformations.extend(format_sources(sources, package.path, timeout))
        if checks_passed(transformations):
            transformations.append(run_command(xpure_command(xpure, sources, False), package.path, timeout))
        if checks_passed(transformations):
            transformations.append(run_command(xpure_command(xpure, sources, True), package.path, timeout))
        if checks_passed(transformations):
            transformations.extend(format_sources(sources, package.path, timeout))

        if not checks_passed(transformations):
            restore_sources(snapshot, package.path, sources)
            write_log(log_path, f"Transformation-tool failure; sources restored: {package.name}", transformations)
            print(f"TOOL FAIL {package.name}; restored; log: {log_path}")
            return False

        changed = changed_sources(snapshot, package.path, sources)
        if not changed:
            print(f"UNCHANGED {package.name}")
            if log_path.exists():
                log_path.unlink()
            return True

        verification = fpm_checks(package.path, timeout)
        if not checks_passed(verification):
            restore_sources(snapshot, package.path, sources)
            write_log(log_path, f"Post-refactoring regression; sources restored: {package.name}", verification)
            print(f"REGRESSION {package.name}; restored; log: {log_path}")
            return False

        if log_path.exists():
            log_path.unlink()
        print(f"PASS {package.name}: retained changes in {len(changed)} file(s)")

    if clean:
        cleanup = run_command(["fpm", "clean", "--all"], package.path, timeout)
        if cleanup.returncode != 0:
            write_log(log_path, f"Refactoring passed, but cleanup failed: {package.name}", [cleanup])
            print(f"WARNING: cleanup failed for {package.name}; log: {log_path}")
    return True


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packages", nargs="*", help="top-level package directory names")
    parser.add_argument("--all", action="store_true", help="process every package in dependency-first order")
    parser.add_argument("--list", action="store_true", help="show dependency-first order without changing files")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--xpure", type=Path, default=DEFAULT_XPURE)
    parser.add_argument("--xparam", type=Path, default=DEFAULT_XPARAM)
    parser.add_argument("--timeout", type=float, default=900.0, help="seconds per command (default: 900)")
    parser.add_argument("--keep-build", action="store_true", help="do not run 'fpm clean --all' after success")
    parser.add_argument("--continue-on-failure", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.root.resolve()
    if args.timeout <= 0:
        print("ERROR: --timeout must be positive", file=sys.stderr)
        return 2
    if not (root / "remove_fortran_semicolons.py").is_file():
        print("ERROR: remove_fortran_semicolons.py is missing", file=sys.stderr)
        return 2
    if not (root / "modernize_fortran_syntax.py").is_file():
        print("ERROR: modernize_fortran_syntax.py is missing", file=sys.stderr)
        return 2
    if not args.xpure.is_file():
        print(f"ERROR: xpure.py was not found: {args.xpure}", file=sys.stderr)
        return 2
    if not args.xparam.is_file():
        print(f"ERROR: xparam.py was not found: {args.xparam}", file=sys.stderr)
        return 2
    for executable in ("fpm", "fprettify"):
        if shutil.which(executable) is None:
            print(f"ERROR: {executable} is not on PATH", file=sys.stderr)
            return 2

    packages = discover_packages(root)
    ordered = dependency_order(packages)
    if args.list:
        for number, package in enumerate(ordered, start=1):
            dependencies = ", ".join(package.dependencies) or "-"
            print(f"{number:>3}  {package.name:<32} dependencies: {dependencies}")
        print(f"\n{len(ordered)} package(s).")
        return 0

    if args.all and args.packages:
        print("ERROR: specify package names or --all, not both", file=sys.stderr)
        return 2
    if not args.all and not args.packages:
        print("ERROR: specify one or more package names, or use --all/--list", file=sys.stderr)
        return 2

    if args.all:
        selected = ordered
    else:
        requested = {name.casefold() for name in args.packages}
        missing = sorted(requested - set(packages))
        if missing:
            print(f"ERROR: unknown package(s): {', '.join(missing)}", file=sys.stderr)
            return 2
        selected = [package for package in ordered if package.name.casefold() in requested]

    failures = 0
    for number, package in enumerate(selected, start=1):
        print(f"\n[{number}/{len(selected)}] {package.name}")
        if not transform_package(
            package,
            root,
            args.xpure.resolve(),
            args.xparam.resolve(),
            args.timeout,
            not args.keep_build,
        ):
            failures += 1
            if not args.continue_on_failure:
                break
    print(f"\nCompleted {len(selected) if not failures or args.continue_on_failure else number} package(s); failures: {failures}.")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
