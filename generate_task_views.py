#!/usr/bin/env python3
"""Generate filtered task views for translated top-level FPM packages."""

from __future__ import annotations

import argparse
import re
import urllib.error
import urllib.request
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path


PACKAGE_PATTERN = re.compile(r"pkg\([\"']([^\"']+)[\"']")
README_ROW_PATTERN = re.compile(
    r"^\| \[`([^`]+)`\]\([^)]*\) \| (.*?) \|", re.MULTILINE
)
HEADING_PATTERN = re.compile(r"^(#{1,4})\s+(.+?)\s*$")


@dataclass(frozen=True)
class ViewSpec:
    name: str
    filename: str
    skip_sections: frozenset[str] = frozenset({"Introduction", "Links"})

    @property
    def source_url(self) -> str:
        return (
            f"https://raw.githubusercontent.com/cran-task-views/{self.name}/"
            f"main/{self.name}.md"
        )

    @property
    def repository_url(self) -> str:
        return f"https://github.com/cran-task-views/{self.name}/blob/main/{self.name}.md"

    @property
    def cran_url(self) -> str:
        return f"https://CRAN.R-project.org/view={self.name}"


VIEW_SPECS = (
    ViewSpec("Finance", "Finance.md"),
    ViewSpec("Optimization", "Optimization.md"),
    ViewSpec("NumericalMathematics", "NumericalMathematics.md"),
    ViewSpec(
        "Distributions",
        "Distributions.md",
        frozenset(
            {
                "Introduction",
                "Table of contents",
                "Bibliography",
                "General books",
                "Books dedicated to a distribution family",
                "Books with applications",
                "Links",
            }
        ),
    ),
)


def clean_heading(text: str) -> str:
    text = re.sub(r"\{#[^}]+\}", "", text).strip()
    match = re.fullmatch(r"\[([^]]+)\](?:\([^)]*\))?", text)
    if match:
        text = match.group(1)
    return text.strip()


def metadata(source: str) -> dict[str, str]:
    result = {}
    for key in ("name", "topic", "maintainer", "version"):
        match = re.search(rf"^{key}:\s*(.+?)\s*$", source, re.MULTILINE)
        if not match:
            raise ValueError(f"source task view has no {key!r} metadata")
        result[key] = match.group(1)
    return result


def translated_packages(root: Path) -> dict[str, str]:
    return {
        path.name.casefold(): path.name
        for path in root.iterdir()
        if path.is_dir() and (path / "fpm.toml").is_file()
    }


def readme_descriptions(root: Path) -> dict[str, str]:
    text = (root / "README.md").read_text(encoding="utf-8")
    return {
        name.casefold(): description.strip()
        for name, description in README_ROW_PATTERN.findall(text)
    }


def group_packages(
    source: str, local: dict[str, str], skip_sections: frozenset[str]
) -> tuple[OrderedDict[str, list[str]], set[str]]:
    groups: OrderedDict[str, list[str]] = OrderedDict()
    current_section = "Introduction"
    all_matches: set[str] = set()
    assigned: set[str] = set()

    for line in source.splitlines():
        heading = HEADING_PATTERN.match(line)
        if heading:
            current_section = clean_heading(heading.group(2))
        for upstream_name in PACKAGE_PATTERN.findall(line):
            local_name = local.get(upstream_name.casefold())
            if local_name is None:
                continue
            all_matches.add(local_name)
            if current_section in skip_sections or local_name in assigned:
                continue
            groups.setdefault(current_section, []).append(local_name)
            assigned.add(local_name)

    missing = all_matches - assigned
    return groups, missing


def download(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "r-fpm-task-view-generator"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return response.read().decode("utf-8")


def render_view(
    spec: ViewSpec,
    source: str,
    local: dict[str, str],
    descriptions: dict[str, str],
) -> tuple[str, int]:
    info = metadata(source)
    groups, missing = group_packages(source, local, spec.skip_sections)
    if missing:
        raise ValueError(f"{spec.name}: matched but unassigned packages: {sorted(missing)}")

    package_names = [name for names in groups.values() for name in names]
    missing_descriptions = [name for name in package_names if name.casefold() not in descriptions]
    if missing_descriptions:
        raise ValueError(
            f"{spec.name}: packages missing from the README table: {missing_descriptions}"
        )

    lines = [
        f"# {info['topic']}",
        "",
        f"This is an independent, filtered adaptation of the",
        f"[CRAN Task View: {info['topic']}]({spec.cran_url}), maintained by",
        f"{info['maintainer']}, version {info['version']}. The",
        f"[source task view]({spec.repository_url}) provides the broader annotated",
        "guide to R packages.",
        "",
        f"This page includes {len(package_names)} translated packages. Package membership and broad",
        "topic organization follow the source task view; the concise descriptions are",
        "the high-level summaries maintained by this project. They describe translated",
        "computational scope and intentionally omit plotting, interactive interfaces, R",
        "classes, and internal Fortran organization. See each package README for precise",
        "API coverage and validation status.",
        "",
        "<!-- Generated by generate_task_views.py. Update package summaries in README.md. -->",
    ]
    for section, names in groups.items():
        lines.extend(("", f"## {section}", ""))
        for name in names:
            description = descriptions[name.casefold()]
            lines.append(f"- [`{name}`](../{name}/) - {description}")
    lines.extend(
        (
            "",
            "Packages are listed under their first applicable source-view section to keep",
            "the filtered view compact; many packages support methods relevant to additional",
            "sections.",
            "",
        )
    )
    return "\n".join(lines), len(package_names)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "views",
        nargs="*",
        help="view names to generate; omit to generate all configured views",
    )
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parent, help="repository root"
    )
    parser.add_argument("--check", action="store_true", help="report stale files without writing")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    selected = {name.casefold() for name in args.views}
    specs = [spec for spec in VIEW_SPECS if not selected or spec.name.casefold() in selected]
    unknown = selected - {spec.name.casefold() for spec in VIEW_SPECS}
    if unknown:
        print(f"ERROR: unknown task view(s): {', '.join(sorted(unknown))}")
        return 2

    local = translated_packages(root)
    descriptions = readme_descriptions(root)
    output_dir = root / "TASK_VIEWS"
    stale = 0
    for spec in specs:
        try:
            source = download(spec.source_url)
            rendered, count = render_view(spec, source, local, descriptions)
        except (OSError, ValueError, urllib.error.URLError) as error:
            print(f"ERROR: {spec.name}: {error}")
            return 1
        output = output_dir / spec.filename
        old = output.read_text(encoding="utf-8") if output.exists() else None
        if old == rendered:
            print(f"CURRENT {spec.name}: {count} package(s)")
        elif args.check:
            print(f"STALE   {spec.name}: {count} package(s)")
            stale += 1
        else:
            output.write_text(rendered, encoding="utf-8", newline="\n")
            print(f"UPDATED {spec.name}: {count} package(s)")
    return 1 if stale else 0


if __name__ == "__main__":
    raise SystemExit(main())
