from __future__ import annotations

import re
import sys
from pathlib import Path

UNIT_RE = re.compile(
    r"^\s*(?:(?:real\s*\([^)]*\)|complex\s*\([^)]*\)|double\s+precision|double\s+complex|integer|real|logical|complex|character(?:\s*\([^)]*\)|\s*\*\s*\d+)?)\s+)?"
    r"(?:recursive\s+|pure\s+|elemental\s+|impure\s+|module\s+)*"
    r"(subroutine|function|program|block\s+data)\b",
    re.IGNORECASE,
)


def lowercase_code(s: str) -> str:
    """Lowercase outside character literals and inline comments."""
    out: list[str] = []
    quote: str | None = None
    i = 0
    while i < len(s):
        ch = s[i]
        if quote is None:
            if ch in ("'", '"'):
                quote = ch
                out.append(ch)
            elif ch == '!':
                out.append(s[i:])
                break
            else:
                out.append(ch.lower())
        else:
            out.append(ch)
            if ch == quote:
                # Doubled quote escapes itself in Fortran strings.
                if i + 1 < len(s) and s[i + 1] == quote:
                    out.append(s[i + 1])
                    i += 1
                else:
                    quote = None
        i += 1
    return ''.join(out)


def modernize_types(s: str) -> str:
    # Only applied to code outside comments; strings are very unlikely to contain
    # declarations, but protect them by operating segment-wise.
    parts: list[str] = []
    quote: str | None = None
    buf: list[str] = []
    i = 0
    while i < len(s):
        ch = s[i]
        if quote is None and ch in ("'", '"'):
            text = ''.join(buf)
            text = re.sub(r'\bdouble\s+complex\b', 'complex(dp)', text, flags=re.I)
            text = re.sub(r'\bdouble\s+precision\b', 'real(dp)', text, flags=re.I)
            parts.append(text)
            buf = [ch]
            quote = ch
        elif quote is not None:
            buf.append(ch)
            if ch == quote:
                if i + 1 < len(s) and s[i + 1] == quote:
                    buf.append(s[i + 1]); i += 1
                else:
                    parts.append(''.join(buf)); buf = []; quote = None
        elif ch == '!':
            text = ''.join(buf)
            text = re.sub(r'\bdouble\s+complex\b', 'complex(dp)', text, flags=re.I)
            text = re.sub(r'\bdouble\s+precision\b', 'real(dp)', text, flags=re.I)
            parts.append(text)
            parts.append(s[i:])
            return ''.join(parts)
        else:
            buf.append(ch)
        i += 1
    text = ''.join(buf)
    if quote is None:
        text = re.sub(r'\bdouble\s+complex\b', 'complex(dp)', text, flags=re.I)
        text = re.sub(r'\bdouble\s+precision\b', 'real(dp)', text, flags=re.I)
    parts.append(text)
    return ''.join(parts)


def fixed_to_free(text: str) -> list[str]:
    out: list[str] = []
    last_code_index: int | None = None

    for raw in text.splitlines():
        line = raw.rstrip('\r\n')
        if not line.strip():
            out.append('')
            continue

        first = line[0]
        if first in 'cC*!':
            out.append('!' + line[1:])
            continue

        # Source contains no active tab-format lines. Expand defensively.
        line = line.expandtabs(8)
        if len(line) < 6:
            line = line.ljust(6)
        label = line[:5].strip()
        cont = line[5:6]
        stmt = line[6:].rstrip()

        is_cont = bool(cont.strip() and cont != '0')
        if is_cont:
            if last_code_index is None:
                raise ValueError(f'continuation without previous statement: {raw!r}')
            # Append a trailing continuation marker before any inline comment.
            prev = out[last_code_index]
            bang = prev.find('!')
            if bang >= 0:
                code, comment = prev[:bang].rstrip(), prev[bang:]
                out[last_code_index] = code + '& ' + comment
            else:
                out[last_code_index] = prev.rstrip() + '&'
            converted = '&' + stmt.lstrip()
        else:
            converted = (label + ' ' if label else '') + stmt.lstrip()

        converted = modernize_types(lowercase_code(converted))
        converted = re.sub(r'\.\s*(eq|ne|lt|le|gt|ge|and|or|not|eqv|neqv)\s*\.', r'.\1.', converted, flags=re.I)
        if re.fullmatch(r'\s*e\s+n\s+d\s*', converted, flags=re.I):
            converted = 'end'
        out.append(converted.rstrip())
        last_code_index = len(out) - 1

    # Insert the kind import after the complete (possibly continued) program-unit
    # statement. Comments and blank lines may occur between continuation lines.
    result: list[str] = []
    pending_unit = False
    for line in out:
        result.append(line)
        code = line.split('!', 1)[0].rstrip()
        if not pending_unit:
            if UNIT_RE.match(code) and not re.match(r'^\s*end\s+', code, re.I):
                if code.endswith('&'):
                    pending_unit = True
                else:
                    result.append('  use timsac_kinds, only: dp')
                    result.append('  implicit none')
        else:
            if code and not code.lstrip().startswith('!'):
                if not code.endswith('&'):
                    result.append('  use timsac_kinds, only: dp')
                    result.append('  implicit none')
                    pending_unit = False
    if pending_unit:
        raise ValueError('unterminated program-unit declaration')
    return result


def main() -> None:
    src = Path(sys.argv[1])
    dst = Path(sys.argv[2])
    dst.parent.mkdir(parents=True, exist_ok=True)
    lines = fixed_to_free(src.read_text(errors='replace'))
    dst.write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
