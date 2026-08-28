#!/usr/bin/env python3
"""Build or test each top-level FPM package and summarize the results."""

from __future__ import annotations

import argparse
import os
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 compatibility
    tomllib = None


FORTRAN_SUFFIXES = {".f", ".for", ".f77", ".f90", ".f95", ".f03", ".f08"}


@dataclass
class Result:
    package: str
    status: str
    seconds: float
    log_path: Path | None = None
    cleanup_failed: bool = False


def package_test_prerequisites(package: Path) -> tuple[list[Path], str]:
    """Return missing package test prerequisites and an optional preparation hint."""
    if tomllib is None:
        return [], ""
    try:
        with (package / "fpm.toml").open("rb") as stream:
            config = tomllib.load(stream)
    except (OSError, tomllib.TOMLDecodeError):
        return [], ""
    test_config = config.get("extra", {}).get("r-fpm-test", {})
    required_files = test_config.get("required-files", [])
    missing = [package / path for path in required_files if not (package / path).is_file()]
    return missing, test_config.get("preparation", "")


def discover_packages(root: Path) -> list[Path]:
    """Return immediate child directories containing an FPM manifest."""
    return sorted(
        (path for path in root.iterdir() if path.is_dir() and (path / "fpm.toml").is_file()),
        key=lambda path: path.name.casefold(),
    )


def package_has_tests(package: Path) -> bool:
    """Use the manifest and conventional test directory to detect FPM tests."""
    manifest = package / "fpm.toml"
    if tomllib is not None:
        try:
            with manifest.open("rb") as stream:
                config = tomllib.load(stream)
            if config.get("test"):
                return True
            if config.get("build", {}).get("auto-tests", True) is False:
                return False
        except (OSError, tomllib.TOMLDecodeError):
            # Let FPM report a malformed or unreadable manifest.
            pass

    test_dir = package / "test"
    if not test_dir.is_dir():
        return False
    return any(
        path.is_file() and path.suffix.casefold() in FORTRAN_SUFFIXES
        for path in test_dir.rglob("*")
    )


def stop_process_tree(process: subprocess.Popen[str]) -> None:
    """Stop a timed-out FPM process and any compiler/test children."""
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


def run_fpm(package: Path, command: list[str], timeout: float) -> tuple[int | None, str, float]:
    """Run FPM, returning (exit code, combined output, elapsed seconds)."""
    creationflags = subprocess.CREATE_NEW_PROCESS_GROUP if os.name == "nt" else 0
    started = time.perf_counter()
    process = subprocess.Popen(
        command,
        cwd=package,
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
        return process.returncode, output, time.perf_counter() - started
    except subprocess.TimeoutExpired:
        stop_process_tree(process)
        output, _ = process.communicate()
        return None, output, time.perf_counter() - started
    except KeyboardInterrupt:
        stop_process_tree(process)
        process.communicate()
        raise


def safe_log_name(package_name: str) -> str:
    return "".join(character if character.isalnum() or character in "._-" else "_" for character in package_name)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "packages",
        nargs="*",
        help="package directory names to test; omit to test every top-level FPM package",
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="repository root (default: directory containing this script)",
    )
    parser.add_argument("--timeout", type=float, default=300.0, help="seconds allowed per package (default: 300)")
    parser.add_argument("--profile", choices=("debug", "release"), default="debug")
    parser.add_argument("--build-only", action="store_true", help="run fpm build even when tests are available")
    parser.add_argument("--stop-on-failure", action="store_true")
    parser.add_argument(
        "--clean",
        choices=("all", "skip", "none"),
        default="all",
        help="cleanup after each package: all removes dependencies, skip preserves them, none preserves the build (default: all)",
    )
    parser.add_argument("--list", action="store_true", help="list selected packages without running FPM")
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=Path("package_test_logs"),
        help="failure-log directory, relative to the root unless absolute",
    )
    return parser.parse_args(argv)


def select_packages(root: Path, requested: list[str]) -> tuple[list[Path], list[str]]:
    discovered = discover_packages(root)
    if not requested:
        return discovered, []
    by_name = {path.name.casefold(): path for path in discovered}
    selected = []
    missing = []
    for name in requested:
        package = by_name.get(name.casefold())
        if package is None:
            missing.append(name)
        elif package not in selected:
            selected.append(package)
    return selected, missing


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.root.resolve()
    if args.timeout <= 0:
        print("ERROR: --timeout must be positive.", file=sys.stderr)
        return 2
    if not root.is_dir():
        print(f"ERROR: repository root does not exist: {root}", file=sys.stderr)
        return 2

    packages, missing = select_packages(root, args.packages)
    if missing:
        print(f"ERROR: no top-level FPM package(s): {', '.join(missing)}", file=sys.stderr)
        return 2
    if not packages:
        print("ERROR: no top-level FPM packages found.", file=sys.stderr)
        return 2
    if args.list:
        for package in packages:
            print(package.name)
        print(f"\n{len(packages)} package(s).")
        return 0
    if shutil.which("fpm") is None:
        print("ERROR: fpm is not on PATH.", file=sys.stderr)
        return 2

    log_dir = args.log_dir if args.log_dir.is_absolute() else root / args.log_dir
    results: list[Result] = []
    total = len(packages)
    print(f"Testing {total} top-level FPM package(s) with a {args.timeout:g}-second timeout.")
    if args.clean == "none":
        print("Build products will be preserved.")
    else:
        print(f"Each package will be cleaned with 'fpm clean --{args.clean}' after it runs.")

    for number, package in enumerate(packages, start=1):
        has_tests = package_has_tests(package) and not args.build_only
        action = "test" if has_tests else "build"
        command = ["fpm", action, "--profile", args.profile]
        print(f"[{number:>3}/{total}] {action.upper():5} {package.name} ... ", end="", flush=True)
        missing_prerequisites, preparation = package_test_prerequisites(package)
        if missing_prerequisites:
            detail = f"; {preparation}" if preparation else ""
            print(f"SKIP (missing {missing_prerequisites[0].relative_to(package)}{detail})")
            stale_log = log_dir / f"{safe_log_name(package.name)}.log"
            if stale_log.exists():
                stale_log.unlink()
            results.append(Result(package.name, "SKIP", 0.0))
            continue
        try:
            returncode, output, seconds = run_fpm(package, command, args.timeout)
        except KeyboardInterrupt:
            print("INTERRUPTED")
            if args.clean != "none":
                print(f"Cleaning {package.name} before stopping ...", flush=True)
                run_fpm(package, ["fpm", "clean", f"--{args.clean}"], args.timeout)
            return 130

        if returncode is None:
            status = "TIMEOUT"
        elif returncode == 0:
            status = "PASS" if has_tests else "BUILD"
        else:
            status = "FAIL"

        cleanup_failed = False
        cleanup_output = ""
        cleanup_seconds = 0.0
        if args.clean != "none":
            cleanup_command = ["fpm", "clean", f"--{args.clean}"]
            cleanup_code, cleanup_output, cleanup_seconds = run_fpm(package, cleanup_command, args.timeout)
            cleanup_failed = cleanup_code != 0

        failure_log_path = log_dir / f"{safe_log_name(package.name)}.log"
        log_path = None
        if status in {"FAIL", "TIMEOUT"} or cleanup_failed:
            log_dir.mkdir(parents=True, exist_ok=True)
            log_path = failure_log_path
            command_text = subprocess.list2cmdline(command)
            cleanup_text = ""
            if args.clean != "none":
                cleanup_text = (
                    f"\n\nCleanup command: fpm clean --{args.clean}"
                    f"\nCleanup status: {'FAILED' if cleanup_failed else 'passed'}"
                    f"\n\n{cleanup_output}"
                )
            log_path.write_text(
                f"Package: {package.name}\nCommand: {command_text}\nStatus: {status}\n\n{output}{cleanup_text}",
                encoding="utf-8",
            )
        elif failure_log_path.exists():
            failure_log_path.unlink()
        results.append(Result(package.name, status, seconds + cleanup_seconds, log_path, cleanup_failed))
        suffix = f"; log: {log_path}" if log_path is not None else ""
        cleanup_suffix = "; CLEANUP FAILED" if cleanup_failed else ""
        print(f"{status} ({seconds:.1f}s){cleanup_suffix}{suffix}")
        if args.stop_on_failure and (status in {"FAIL", "TIMEOUT"} or cleanup_failed):
            break

    statuses = ("PASS", "BUILD", "SKIP", "FAIL", "TIMEOUT")
    counts = {status: sum(result.status == status for result in results) for status in statuses}
    cleanup_failures = sum(result.cleanup_failed for result in results)
    elapsed = sum(result.seconds for result in results)
    print("\nSummary")
    print(f"  Packages run:       {len(results)} of {total}")
    print(f"  Tests passed:       {counts['PASS']}")
    print(f"  Build-only passed:  {counts['BUILD']}")
    print(f"  Prerequisite skips: {counts['SKIP']}")
    print(f"  Failed:             {counts['FAIL']}")
    print(f"  Timed out:          {counts['TIMEOUT']}")
    print(f"  Cleanup failures:   {cleanup_failures}")
    print(f"  Total package time: {elapsed:.1f}s")
    if counts["FAIL"] or counts["TIMEOUT"] or cleanup_failures:
        print(f"  Failure logs:       {log_dir}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
