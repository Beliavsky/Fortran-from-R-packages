#!/usr/bin/env python3
"""Apply conservative, semantics-preserving Fortran syntax cleanups.

Currently implemented transformations:

* remove a procedure-level ``implicit none`` inherited from its enclosing
  module or procedure;
* remove names from block DO constructs when no EXIT or CYCLE targets them.
* make the implicit ``SAVE`` semantics of initialized local variables explicit.

Interface bodies, standalone procedures, legacy labelled DO loops, imported
source trees, fixed-form source, comments, and character literals are left
alone.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path


FORTRAN_SUFFIXES = {".f90", ".f95", ".f03", ".f08"}
EXCLUDED_PARTS = {".git", ".fpm", "build", "dependencies", "orig", "original", "upstream", "vendor"}
IMPLICIT_NONE_RE = re.compile(r"^\s*implicit\s+none(?:\s*\([^)]*\))?\s*$", re.IGNORECASE)
MODULE_START_RE = re.compile(r"^\s*module\s+(?!procedure\b|subroutine\b|function\b)([a-z]\w*)\b", re.IGNORECASE)
PROGRAM_START_RE = re.compile(r"^\s*program\s+([a-z]\w*)\b", re.IGNORECASE)
PROCEDURE_START_RE = re.compile(
    r"^\s*(?:(?:pure|impure|elemental|recursive|module|non_recursive)\s+)*"
    r"(?:(?:integer|real|logical|complex)(?:\s*\([^)]*\)|\s*\*\s*\d+)?\s+"
    r"|(?:double\s+precision|character(?:\s*\([^)]*\))?"
    r"|type\s*\([^)]*\)|class\s*\([^)]*\))\s+)?"
    r"(subroutine|function)\s+([a-z]\w*)\b",
    re.IGNORECASE,
)
END_SCOPE_RE = re.compile(
    r"^\s*end(?:\s+(module|program|subroutine|function)(?:\s+[a-z]\w*)?)?\s*$",
    re.IGNORECASE,
)
INTERFACE_START_RE = re.compile(r"^\s*(?:abstract\s+)?interface\b", re.IGNORECASE)
INTERFACE_END_RE = re.compile(r"^\s*end\s*interface\b", re.IGNORECASE)
DO_START_RE = re.compile(r"^(?P<indent>\s*)(?:(?P<name>[a-z]\w*)\s*:\s*)?do\b(?P<tail>.*)$", re.IGNORECASE)
DO_END_RE = re.compile(r"^(?P<indent>\s*)end\s*do(?:\s+(?P<name>[a-z]\w*))?\s*$", re.IGNORECASE)
TYPE_SPEC_RE = re.compile(
    r"^(?:integer|real|logical|complex)(?:\s*\([^)]*\)|\s*\*\s*\w+)?\b"
    r"|^double\s+precision\b"
    r"|^character(?:\s*\([^)]*\)|\s*\*\s*\w+)?\b"
    r"|^(?:type|class)\s*\([^)]*\)",
    re.IGNORECASE,
)
DERIVED_TYPE_START_RE = re.compile(
    r"^\s*type\s*(?:,\s*[^:]*)?::\s*[a-z]\w*\b|^\s*type\s+[a-z]\w*\s*$",
    re.IGNORECASE,
)
DERIVED_TYPE_END_RE = re.compile(r"^\s*end\s*type\b", re.IGNORECASE)


@dataclass
class Scope:
    kind: str
    effective_implicit_none: bool
    inherited_implicit_none: bool


@dataclass
class DoConstruct:
    start: int
    end: int
    name: str | None
    end_name: str | None


def split_code_comment(line: str) -> tuple[str, str]:
    """Split a physical line at a comment marker outside a character literal."""
    in_single = False
    in_double = False
    i = 0
    while i < len(line):
        character = line[i]
        if character == "'" and not in_double:
            if in_single and i + 1 < len(line) and line[i + 1] == "'":
                i += 2
                continue
            in_single = not in_single
        elif character == '"' and not in_single:
            if in_double and i + 1 < len(line) and line[i + 1] == '"':
                i += 2
                continue
            in_double = not in_double
        elif character == "!" and not in_single and not in_double:
            return line[:i], line[i:]
        i += 1
    return line, ""


def eol_of(line: str) -> str:
    if line.endswith("\r\n"):
        return "\r\n"
    if line.endswith("\n"):
        return "\n"
    return ""


def mask_character_literals(code: str) -> str:
    """Replace character-literal contents so keywords inside them are ignored."""
    characters = list(code)
    quote: str | None = None
    i = 0
    while i < len(characters):
        character = characters[i]
        if quote is None and character in {"'", '"'}:
            quote = character
            characters[i] = " "
        elif quote is not None:
            characters[i] = " "
            if character == quote:
                if i + 1 < len(characters) and characters[i + 1] == quote:
                    characters[i + 1] = " "
                    i += 1
                else:
                    quote = None
        i += 1
    return "".join(characters)


def source_paths(paths: list[Path]) -> list[Path]:
    """Expand files and directories to maintained free-form sources."""
    selected: set[Path] = set()
    for path in paths:
        if path.is_file():
            candidates = [path]
            base = path.parent
        elif path.is_dir():
            candidates = path.rglob("*")
            base = path
        else:
            continue
        for candidate in candidates:
            if not candidate.is_file() or candidate.suffix.casefold() not in FORTRAN_SUFFIXES:
                continue
            relative_parts = {part.casefold() for part in candidate.relative_to(base).parts[:-1]}
            if relative_parts & EXCLUDED_PARTS:
                continue
            selected.add(candidate.resolve())
    return sorted(selected, key=lambda path: str(path).casefold())


def remove_redundant_implicit_none(lines: list[str]) -> tuple[list[str], int]:
    """Remove procedure IMPLICIT NONE inherited through host association."""
    updated = list(lines)
    scopes: list[Scope] = []
    interface_depth = 0
    removed = 0

    for index, line in enumerate(lines):
        raw = line.rstrip("\r\n")
        code, comment = split_code_comment(raw)
        statement = code.strip()
        if not statement:
            continue

        if INTERFACE_END_RE.match(statement):
            interface_depth = max(0, interface_depth - 1)
            continue
        if INTERFACE_START_RE.match(statement):
            interface_depth += 1
            continue
        if interface_depth:
            continue

        module_match = MODULE_START_RE.match(statement)
        program_match = PROGRAM_START_RE.match(statement)
        procedure_match = PROCEDURE_START_RE.match(statement)
        if module_match:
            scopes.append(Scope("module", False, False))
            continue
        if program_match:
            scopes.append(Scope("program", False, False))
            continue
        if procedure_match:
            inherited = bool(scopes and scopes[-1].effective_implicit_none)
            scopes.append(Scope("procedure", inherited, inherited))
            continue

        if IMPLICIT_NONE_RE.match(statement) and scopes:
            scope = scopes[-1]
            if scope.kind == "procedure" and scope.inherited_implicit_none:
                indent = code[: len(code) - len(code.lstrip())]
                ending = eol_of(line)
                updated[index] = f"{indent}{comment.lstrip()}{ending}" if comment else ""
                removed += 1
            else:
                scope.effective_implicit_none = True
            continue

        end_match = END_SCOPE_RE.match(statement)
        if end_match and scopes:
            end_kind = (end_match.group(1) or "").casefold()
            if not end_kind or end_kind == scopes[-1].kind or (
                scopes[-1].kind == "procedure" and end_kind in {"subroutine", "function"}
            ):
                scopes.pop()

    return updated, removed


def find_do_constructs(lines: list[str]) -> list[DoConstruct]:
    """Pair block DO starts and ends, including unnamed nesting."""
    stack: list[tuple[int, str | None]] = []
    constructs: list[DoConstruct] = []
    for index, line in enumerate(lines):
        code, _ = split_code_comment(line.rstrip("\r\n"))
        end_match = DO_END_RE.match(code)
        if end_match:
            if stack:
                start, name = stack.pop()
                constructs.append(DoConstruct(start, index, name, end_match.group("name")))
            continue
        start_match = DO_START_RE.match(code)
        if not start_match:
            continue
        if re.match(r"^\s*\d+\s+", start_match.group("tail")):
            continue
        stack.append((index, start_match.group("name")))
    return constructs


def remove_unused_do_names(lines: list[str]) -> tuple[list[str], int]:
    """Remove block-DO names that are not targeted by EXIT or CYCLE."""
    updated = list(lines)
    removed = 0
    code_lines = [mask_character_literals(split_code_comment(line.rstrip("\r\n"))[0]) for line in lines]
    for construct in find_do_constructs(lines):
        if not construct.name:
            continue
        name = construct.name
        if construct.end_name and construct.end_name.casefold() != name.casefold():
            continue
        target_re = re.compile(rf"\b(?:exit|cycle)\s+{re.escape(name)}\b", re.IGNORECASE)
        if any(target_re.search(code_lines[index]) for index in range(construct.start + 1, construct.end)):
            continue

        start_line = updated[construct.start]
        start_raw = start_line.rstrip("\r\n")
        start_code, start_comment = split_code_comment(start_raw)
        start_match = DO_START_RE.match(start_code)
        if start_match is None:
            continue
        updated[construct.start] = (
            f"{start_match.group('indent')}do{start_match.group('tail')}{start_comment}{eol_of(start_line)}"
        )

        end_line = updated[construct.end]
        end_raw = end_line.rstrip("\r\n")
        end_code, end_comment = split_code_comment(end_raw)
        end_match = DO_END_RE.match(end_code)
        if end_match is not None:
            updated[construct.end] = f"{end_match.group('indent')}end do{end_comment}{eol_of(end_line)}"
        removed += 1
    return updated, removed


def declaration_has_initializer(entities: str) -> bool:
    """Return whether an entity list contains explicit initialization."""
    return any("=" in entity for entity in split_top_level_commas(entities))


def join_initialized_continued_declarations(
    lines: list[str], max_line_length: int = 120
) -> list[str]:
    """Join simple continued declarations so entity lifetimes can be split safely."""
    updated: list[str] = []
    index = 0
    while index < len(lines):
        first = lines[index]
        first_raw = first.rstrip("\r\n")
        first_code, first_comment = split_code_comment(first_raw)
        statement = first_code.strip()
        if (
            first_comment
            or not first_code.rstrip().endswith("&")
            or "::" not in statement
            or TYPE_SPEC_RE.match(statement) is None
            or "'" in first_code
            or '"' in first_code
        ):
            updated.append(first)
            index += 1
            continue

        pieces = [first_code.rstrip()[:-1].rstrip()]
        final_index = index
        safe = True
        complete = False
        while final_index + 1 < len(lines):
            final_index += 1
            next_raw = lines[final_index].rstrip("\r\n")
            next_code, next_comment = split_code_comment(next_raw)
            if next_comment or "'" in next_code or '"' in next_code:
                safe = False
                break
            continuation = next_code.strip()
            if continuation.startswith("&"):
                continuation = continuation[1:].lstrip()
            continues = continuation.rstrip().endswith("&")
            if continues:
                continuation = continuation.rstrip()[:-1].rstrip()
            pieces.append(continuation)
            if not continues:
                complete = True
                break

        combined = " ".join(pieces)
        if safe and complete:
            _attributes, entities = combined.strip().split("::", 1)
            has_array_constructor = "[" in entities or "(/" in entities
            reasonably_sized = len(combined) <= 4 * max_line_length
            if (
                declaration_has_initializer(entities)
                and not has_array_constructor
                and reasonably_sized
            ):
                updated.append(f"{combined}{eol_of(lines[final_index])}")
                index = final_index + 1
                continue
        updated.extend(lines[index : final_index + 1])
        index = final_index + 1
    return updated


def add_explicit_save_to_initialized_locals(lines: list[str]) -> tuple[list[str], int]:
    """Add SAVE where local initialization already implies SAVE semantics."""
    updated: list[str] = []
    scopes: list[str] = []
    interface_depth = 0
    derived_type_depth = 0
    changed = 0

    for index, line in enumerate(lines):
        raw = line.rstrip("\r\n")
        code, comment = split_code_comment(raw)
        statement = code.strip()
        if not statement:
            updated.append(line)
            continue

        if INTERFACE_END_RE.match(statement):
            interface_depth = max(0, interface_depth - 1)
            updated.append(line)
            continue
        if INTERFACE_START_RE.match(statement):
            interface_depth += 1
            updated.append(line)
            continue
        if interface_depth:
            updated.append(line)
            continue

        if DERIVED_TYPE_END_RE.match(statement):
            derived_type_depth = max(0, derived_type_depth - 1)
            updated.append(line)
            continue
        if DERIVED_TYPE_START_RE.match(statement):
            derived_type_depth += 1
            updated.append(line)
            continue

        module_match = MODULE_START_RE.match(statement)
        program_match = PROGRAM_START_RE.match(statement)
        procedure_match = PROCEDURE_START_RE.match(statement)
        if module_match:
            scopes.append("module")
            updated.append(line)
            continue
        if program_match:
            scopes.append("program")
            updated.append(line)
            continue
        if procedure_match:
            scopes.append("procedure")
            updated.append(line)
            continue

        end_match = END_SCOPE_RE.match(statement)
        if end_match and scopes:
            end_kind = (end_match.group(1) or "").casefold()
            if not end_kind or end_kind == scopes[-1] or (
                scopes[-1] == "procedure" and end_kind in {"subroutine", "function"}
            ):
                scopes.pop()
            updated.append(line)
            continue

        if not scopes or scopes[-1] not in {"procedure", "program"} or derived_type_depth:
            updated.append(line)
            continue
        type_match = TYPE_SPEC_RE.match(statement)
        if type_match is None or "::" not in statement or "&" in code:
            updated.append(line)
            continue
        attributes, entities = statement.split("::", 1)
        if re.search(r"\b(?:parameter|save)\b", attributes, re.IGNORECASE):
            updated.append(line)
            continue
        if not declaration_has_initializer(entities):
            updated.append(line)
            continue

        indent = code[: len(code) - len(code.lstrip())]
        entity_groups: list[tuple[bool, list[str]]] = []
        for entity in split_top_level_commas(entities):
            initialized = "=" in entity
            if entity_groups and entity_groups[-1][0] == initialized:
                entity_groups[-1][1].append(entity)
            else:
                entity_groups.append((initialized, [entity]))
        ending = eol_of(line)
        for group_index, (initialized, group) in enumerate(entity_groups):
            save_attribute = ", save" if initialized else ""
            trailing_comment = comment if group_index == len(entity_groups) - 1 else ""
            updated.append(
                f"{indent}{attributes.rstrip()}{save_attribute} :: "
                f"{', '.join(group)}{trailing_comment}{ending}"
            )
            if initialized:
                changed += 1

    return updated, changed


def split_top_level_commas(text: str) -> list[str]:
    """Split an entity list at commas outside delimiters and literals."""
    parts: list[str] = []
    start = 0
    delimiters: list[str] = []
    quote: str | None = None
    pairs = {")": "(", "]": "[", "}": "{"}
    for index, character in enumerate(text):
        if quote is not None:
            if character == quote:
                if index + 1 < len(text) and text[index + 1] == quote:
                    continue
                quote = None
            continue
        if character in {"'", '"'}:
            quote = character
        elif character in "([{":
            delimiters.append(character)
        elif character in ")]}" and delimiters and delimiters[-1] == pairs[character]:
            delimiters.pop()
        elif character == "," and not delimiters:
            parts.append(text[start:index].strip())
            start = index + 1
    parts.append(text[start:].strip())
    return [part for part in parts if part]


def split_long_declarations(lines: list[str], max_length: int) -> tuple[list[str], int]:
    """Split unambiguous overlong entity lists into declaration statements."""
    updated: list[str] = []
    split_count = 0
    # Declarations may gain module/procedure indentation in the formatter.
    declaration_limit = max(40, max_length - 12)
    for line in lines:
        raw = line.rstrip("\r\n")
        if len(raw) <= declaration_limit:
            updated.append(line)
            continue
        code, comment = split_code_comment(raw)
        if comment or code.rstrip().endswith("&"):
            updated.append(line)
            continue
        indent = code[: len(code) - len(code.lstrip())]
        body = code.strip()
        type_match = TYPE_SPEC_RE.match(body)
        if type_match is None:
            updated.append(line)
            continue

        if "::" in body:
            declaration_prefix, entities_text = body.split("::", 1)
            prefix = declaration_prefix.rstrip() + " :: "
        else:
            type_spec = type_match.group(0).rstrip()
            remainder = body[type_match.end():]
            if not remainder or not remainder[0].isspace():
                updated.append(line)
                continue
            prefix = type_spec + " :: "
            entities_text = remainder.strip()

        entities = split_top_level_commas(entities_text)
        if len(entities) < 2:
            updated.append(line)
            continue
        groups: list[list[str]] = []
        current: list[str] = []
        for entity in entities:
            candidate = ", ".join([*current, entity])
            if current and len(indent) + len(prefix) + len(candidate) > declaration_limit:
                groups.append(current)
                current = [entity]
            else:
                current.append(entity)
        if current:
            groups.append(current)
        if len(groups) < 2:
            updated.append(line)
            continue
        ending = eol_of(line)
        updated.extend(f"{indent}{prefix}{', '.join(group)}{ending}" for group in groups)
        split_count += 1
    return updated, split_count


def transform_text(text: str, max_line_length: int = 120) -> tuple[str, int, int, int, int]:
    """Apply all transformations and return text plus change counts."""
    lines = text.splitlines(keepends=True)
    lines, implicit_count = remove_redundant_implicit_none(lines)
    lines, do_count = remove_unused_do_names(lines)
    lines = join_initialized_continued_declarations(lines, max_line_length)
    lines, save_count = add_explicit_save_to_initialized_locals(lines)
    lines, declaration_count = split_long_declarations(lines, max_line_length)
    return "".join(lines), implicit_count, do_count, save_count, declaration_count


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true", help="rewrite changed files")
    mode.add_argument("--check", action="store_true", help="return 1 when changes are available")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--max-line-length", type=int, default=120)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    files = source_paths(args.paths)
    if not files:
        print("No maintained free-form Fortran sources found.")
        return 2

    changed_files = 0
    implicit_total = 0
    do_total = 0
    save_total = 0
    declaration_total = 0
    for path in files:
        original = path.read_text(encoding="utf-8")
        transformed, implicit_count, do_count, save_count, declaration_count = transform_text(
            original, args.max_line_length
        )
        if transformed == original:
            continue
        changed_files += 1
        implicit_total += implicit_count
        do_total += do_count
        save_total += save_count
        declaration_total += declaration_count
        if args.verbose:
            print(
                f"{path}: implicit-none={implicit_count}, do-names={do_count}, "
                f"explicit-save={save_count}, long-declarations={declaration_count}"
            )
        if args.write:
            path.write_text(transformed, encoding="utf-8", newline="")

    action = "Updated" if args.write else "Would update"
    print(
        f"{action} {changed_files} file(s): removed {implicit_total} redundant IMPLICIT NONE statement(s), "
        f"removed {do_total} unused DO construct name(s), added {save_total} explicit SAVE attribute(s), "
        f"and split {declaration_total} long declaration(s)."
    )
    if args.check and changed_files:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
