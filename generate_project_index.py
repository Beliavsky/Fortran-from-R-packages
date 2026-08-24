#!/usr/bin/env python3
"""Generate subject-tag and NIST GAMS Markdown indexes for FPM projects."""

from __future__ import annotations

import argparse
import re
import sys
import tomllib
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote


REPOSITORY_ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = "PROJECT_INDEX.md"
DEFAULT_GAMS_OUTPUT = "GAMS_INDEX.md"
GAMS_CLASSES = {
    "A": "Arithmetic, error analysis",
    "B": "Number theory",
    "C": "Elementary and special functions",
    "D": "Linear algebra",
    "E": "Interpolation",
    "F": "Solution of nonlinear equations",
    "G": "Optimization",
    "H": "Differentiation, integration",
    "I": "Differential and integral equations",
    "J": "Integral transforms",
    "K": "Approximation",
    "L": "Statistics, probability",
    "M": "Simulation, stochastic modeling",
    "N": "Data handling",
    "O": "Symbolic computation",
    "P": "Computational geometry",
    "Q": "Graphics",
    "R": "Service routines",
    "S": "Software development tools",
    "Z": "Other",
}


def normalize_tag(tag: str) -> str:
    """Return a stable, case-insensitive spelling for a category or keyword."""
    tag = tag.strip().casefold()
    tag = re.sub(r"[\s_]+", "-", tag)
    return re.sub(r"-+", "-", tag).strip("-")


def markdown_text(text: str) -> str:
    """Escape text used inside a Markdown link label or heading."""
    return text.replace("\\", "\\\\").replace("`", "\\`").replace("|", "\\|")


def read_tags(manifest: Path, field: str, errors: list[str]) -> list[str]:
    try:
        with manifest.open("rb") as stream:
            data = tomllib.load(stream)
    except (OSError, tomllib.TOMLDecodeError) as exc:
        errors.append(f"{manifest}: {exc}")
        return []

    value = data.get(field)
    if value is None:
        return []
    if not isinstance(value, list) or any(not isinstance(item, str) for item in value):
        errors.append(f'{manifest}: "{field}" must be an array of strings')
        return []

    tags = {normalize_tag(item) for item in value}
    tags.discard("")
    return sorted(tags)


def collect_projects(
    root: Path,
) -> tuple[
    dict[str, set[str]],
    dict[str, set[str]],
    dict[str, set[str]],
    list[str],
    list[str],
    list[str],
    list[str],
]:
    categories: dict[str, set[str]] = defaultdict(set)
    keywords: dict[str, set[str]] = defaultdict(set)
    gams: dict[str, set[str]] = defaultdict(set)
    missing_categories: list[str] = []
    missing_keywords: list[str] = []
    missing_gams: list[str] = []
    errors: list[str] = []

    project_dirs = sorted(
        (
            path
            for path in root.iterdir()
            if path.is_dir() and (path / "fpm.toml").is_file()
        ),
        key=lambda path: (path.name.casefold(), path.name),
    )

    for project_dir in project_dirs:
        manifest = project_dir / "fpm.toml"
        project_categories = read_tags(manifest, "categories", errors)
        project_keywords = read_tags(manifest, "keywords", errors)
        project_gams = [code.upper() for code in read_tags(manifest, "gams", errors)]

        if project_categories:
            for category in project_categories:
                categories[category].add(project_dir.name)
        else:
            missing_categories.append(project_dir.name)

        if project_keywords:
            for keyword in project_keywords:
                keywords[keyword].add(project_dir.name)
        else:
            missing_keywords.append(project_dir.name)

        if project_gams:
            for code in project_gams:
                if code not in GAMS_CLASSES:
                    errors.append(f'{manifest}: unknown GAMS class "{code}"')
                else:
                    gams[code].add(project_dir.name)
        else:
            missing_gams.append(project_dir.name)

    return (
        categories,
        keywords,
        gams,
        missing_categories,
        missing_keywords,
        missing_gams,
        errors,
    )


def project_link(project: str) -> str:
    return f"[`{markdown_text(project)}`]({quote(project, safe='')}/)"


def render_grouped_projects(title: str, groups: dict[str, set[str]]) -> list[str]:
    lines = [f"## {title}", ""]
    if not groups:
        lines.extend(["No tags were found.", ""])
        return lines

    for tag in sorted(groups):
        projects = sorted(groups[tag], key=lambda name: (name.casefold(), name))
        lines.extend(
            [
                f"### {markdown_text(tag)}",
                "",
                ", ".join(project_link(project) for project in projects),
                "",
            ]
        )
    return lines


def render_gams_projects(groups: dict[str, set[str]]) -> list[str]:
    lines: list[str] = []
    for code, description in GAMS_CLASSES.items():
        projects = sorted(groups.get(code, set()), key=lambda name: (name.casefold(), name))
        if projects:
            lines.extend(
                [
                    f"### {code} — {description}",
                    "",
                    ", ".join(project_link(project) for project in projects),
                    "",
                ]
            )
    return lines


def count_projects(groups: dict[str, set[str]], missing: list[str]) -> int:
    projects = set(missing)
    for group_projects in groups.values():
        projects.update(group_projects)
    return len(projects)


def render_project_index(
    categories: dict[str, set[str]],
    keywords: dict[str, set[str]],
    missing_categories: list[str],
    missing_keywords: list[str],
) -> str:
    all_projects = set(missing_categories) | set(missing_keywords)
    for projects in categories.values():
        all_projects.update(projects)
    for projects in keywords.values():
        all_projects.update(projects)

    lines = [
        "# Project index",
        "",
        "This file is generated by `generate_project_index.py` from the",
        "`categories` and `keywords` fields in each package's `fpm.toml`.",
        "Do not edit it manually.",
        "",
        f"Projects scanned: {len(all_projects)}.",
        "",
    ]
    lines.extend(render_grouped_projects("Projects by category", categories))
    lines.extend(render_grouped_projects("Projects by keyword", keywords))

    if missing_categories or missing_keywords:
        lines.extend(["## Untagged projects", ""])
        if missing_categories:
            projects = ", ".join(project_link(name) for name in missing_categories)
            lines.extend([f"Missing `categories`: {projects}", ""])
        if missing_keywords:
            projects = ", ".join(project_link(name) for name in missing_keywords)
            lines.extend([f"Missing `keywords`: {projects}", ""])

    return "\n".join(lines).rstrip() + "\n"


def render_gams_index(gams: dict[str, set[str]], missing_gams: list[str]) -> str:
    lines = [
        "# NIST GAMS project index",
        "",
        "This file is generated by `generate_project_index.py` from the `gams` field",
        "in each package's `fpm.toml`. Do not edit it manually.",
        "",
        "The classifications describe broad mathematical-software capabilities and",
        "supplement the repository's subject-based categories and keywords. See the",
        "[NIST GAMS Classification Scheme](https://gams.nist.gov/cgi-bin/serve.cgi).",
        "",
        f"Projects scanned: {count_projects(gams, missing_gams)}.",
        "",
        "## Projects by classification",
        "",
    ]
    lines.extend(render_gams_projects(gams))

    if missing_gams:
        projects = ", ".join(project_link(name) for name in missing_gams)
        lines.extend(["## Unclassified projects", "", f"Missing `gams`: {projects}", ""])

    return "\n".join(lines).rstrip() + "\n"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate subject-tag and NIST GAMS indexes for FPM projects."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(DEFAULT_OUTPUT),
        help=f"output path relative to the repository root (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--gams-output",
        type=Path,
        default=Path(DEFAULT_GAMS_OUTPUT),
        help=(
            "GAMS output path relative to the repository root "
            f"(default: {DEFAULT_GAMS_OUTPUT})"
        ),
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="exit unsuccessfully if either output file is absent or out of date",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    output = args.output if args.output.is_absolute() else REPOSITORY_ROOT / args.output
    gams_output = (
        args.gams_output
        if args.gams_output.is_absolute()
        else REPOSITORY_ROOT / args.gams_output
    )

    (
        categories,
        keywords,
        gams,
        missing_categories,
        missing_keywords,
        missing_gams,
        errors,
    ) = collect_projects(REPOSITORY_ROOT)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    outputs = (
        (
            output,
            render_project_index(
                categories,
                keywords,
                missing_categories,
                missing_keywords,
            ),
        ),
        (gams_output, render_gams_index(gams, missing_gams)),
    )

    if args.check:
        failed = False
        for path, generated in outputs:
            try:
                current = path.read_text(encoding="utf-8")
            except FileNotFoundError:
                print(
                    f"ERROR: {path.name} does not exist; run this script to generate it.",
                    file=sys.stderr,
                )
                failed = True
                continue
            if current != generated:
                print(
                    f"ERROR: {path.name} is out of date; run this script to regenerate it.",
                    file=sys.stderr,
                )
                failed = True
            else:
                print(f"{path.name} is up to date.")
        return int(failed)

    for path, generated in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="\n") as stream:
            stream.write(generated)
        print(f"Wrote {path.name}.")

    if missing_categories:
        print(f"WARNING: {len(missing_categories)} project(s) have no categories.", file=sys.stderr)
    if missing_keywords:
        print(f"WARNING: {len(missing_keywords)} project(s) have no keywords.", file=sys.stderr)
    if missing_gams:
        print(f"WARNING: {len(missing_gams)} project(s) have no GAMS classifications.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
