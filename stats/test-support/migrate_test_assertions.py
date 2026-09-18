"""Replace duplicated local assertion helpers in stats tests with the shared module."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUPPORTED = {
    "assert_approx",
    "assert_close",
    "assert_close_tolerance",
    "assert_equal",
    "assert_integer",
    "assert_integer_matrix",
    "assert_integer_vector",
    "assert_matrix_close",
    "assert_true",
    "assert_vector_approx",
    "assert_vector_close",
    "assert_vector_close_tolerance",
}
REAL_DEFAULT_ASSERTIONS = {"assert_close", "assert_matrix_close", "assert_vector_close"}


def matching_parenthesis(text: str, opening: int) -> int:
    """Return the closing parenthesis while ignoring quoted text."""
    depth = 0
    quote = ""
    index = opening
    while index < len(text):
        character = text[index]
        if quote:
            if character == quote:
                if index + 1 < len(text) and text[index + 1] == quote:
                    index += 1
                else:
                    quote = ""
        elif character in "'\"":
            quote = character
        elif character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
            if depth == 0:
                return index
        index += 1
    raise ValueError("unclosed assertion call")


def argument_count(arguments: str) -> int:
    """Count top-level Fortran actual arguments."""
    if not arguments.strip():
        return 0
    parentheses = brackets = 0
    quote = ""
    count = 1
    index = 0
    while index < len(arguments):
        character = arguments[index]
        if quote:
            if character == quote:
                if index + 1 < len(arguments) and arguments[index + 1] == quote:
                    index += 1
                else:
                    quote = ""
        elif character in "'\"":
            quote = character
        elif character == "(":
            parentheses += 1
        elif character == ")":
            parentheses -= 1
        elif character == "[":
            brackets += 1
        elif character == "]":
            brackets -= 1
        elif character == "," and parentheses == 0 and brackets == 0:
            count += 1
        index += 1
    return count


def preserve_tolerance(text: str) -> str:
    """Pass a test's former host-associated tolerance to three-argument calls."""
    declared = {
        name
        for name in ("tolerance", "parameter_tolerance", "covariance_tolerance")
        if re.search(
            rf"\breal\s*\([^)]*\)\s*,\s*parameter\s*::\s*{name}\b",
            text,
            re.I,
        )
    }
    if not declared:
        return text
    pattern = re.compile(r"\bcall\s+(assert_close|assert_matrix_close|assert_vector_close)\s*\(", re.I)
    additions: list[tuple[int, str]] = []
    for match in pattern.finditer(text):
        closing = matching_parenthesis(text, match.end() - 1)
        arguments = text[match.end():closing]
        if argument_count(arguments) == 3 and "allowed_error" not in arguments.lower():
            procedure = match.group(1).lower()
            if procedure == "assert_matrix_close" and "covariance_tolerance" in declared:
                additions.append((closing, "covariance_tolerance"))
            elif "parameter_tolerance" in declared:
                additions.append((closing, "parameter_tolerance"))
            elif "tolerance" in declared:
                additions.append((closing, "tolerance"))
    for position, tolerance_name in reversed(additions):
        text = text[:position] + f", allowed_error={tolerance_name}" + text[position:]
    return text


def format_use(names: set[str], newline: str) -> str:
    """Build a compact continued ONLY import."""
    ordered = sorted(names)
    lines = ["   use stats_test_assertions, only: "]
    for index, name in enumerate(ordered):
        suffix = ", " if index + 1 < len(ordered) else ""
        if len(lines[-1]) + len(name) + len(suffix) > 96:
            lines[-1] += "&"
            lines.append("                                      ")
        lines[-1] += name + suffix
    return newline.join(lines) + newline


def wrap_tolerance_keywords(text: str, newline: str) -> str:
    """Wrap generated tolerance keywords when a call would exceed 100 columns."""
    output: list[str] = []
    for line in text.splitlines(keepends=True):
        content = line.rstrip("\r\n")
        ending = line[len(content):]
        if len(content) > 100 and ", allowed_error=" in content:
            content = content.replace(", allowed_error=", ", &" + newline + "      allowed_error=", 1)
        output.append(content + ending)
    return "".join(output)


def migrate(path: Path) -> bool:
    """Migrate one test source, returning whether it changed."""
    raw = path.read_bytes()
    newline = "\r\n" if b"\r\n" in raw else "\n"
    text = raw.decode("utf-8")
    use_pattern = re.compile(
        r"^   use stats_test_assertions, only: "
        r"(?:(?:[^\r\n]*&\r?\n[ \t]+)*)[^\r\n]*\r?\n",
        re.I | re.M,
    )
    use_match = use_pattern.search(text)
    imported = set(re.findall(r"\bassert_\w+\b", use_match.group(0), re.I)) if use_match else set()
    if use_match:
        text = text[:use_match.start()] + text[use_match.end():]
    removed: set[str] = set(imported)
    newly_removed = False
    for name in SUPPORTED:
        pattern = re.compile(
            rf"^[ \t]*subroutine {name}\b.*?^[ \t]*end subroutine {name}[ \t]*(?:\r?\n)?",
            re.I | re.M | re.S,
        )
        text, count = pattern.subn("", text)
        if count:
            removed.add(name)
            newly_removed = True
    if not newly_removed and not imported:
        return False
    text = preserve_tolerance(text)
    text = wrap_tolerance_keywords(text, newline)
    program_match = re.search(r"^program\s+\w+[^\r\n]*\r?\n", text, re.I | re.M)
    if not program_match:
        raise ValueError(f"program statement not found in {path}")
    text = text[:program_match.end()] + format_use(removed, newline) + text[program_match.end():]
    text = re.sub(
        r"\r?\ncontains\s*\r?\n\s*(end program\b)",
        lambda match: newline + match.group(1),
        text,
        flags=re.I,
    )
    updated = text.encode("utf-8")
    if updated == raw:
        return False
    path.write_bytes(updated)
    return True


def main() -> None:
    """Migrate every auto-discovered stats test."""
    changed = [path for path in sorted((ROOT / "test").glob("test_*.f90")) if migrate(path)]
    print(f"Migrated {len(changed)} test file(s).")
    for path in changed:
        print(path.relative_to(ROOT))


if __name__ == "__main__":
    main()
