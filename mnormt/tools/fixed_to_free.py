"""Convert the upstream fixed-form probability kernels to free-form Fortran.

This is a provenance/build-maintenance helper, not a runtime dependency.
"""
from pathlib import Path
import re
import sys


def normalize(code: str) -> str:
    code = re.sub(r"\.\s*([A-Za-z]+)\s*\.", lambda m: "." + m.group(1) + ".", code)
    while re.search(r"(?<=[0-9.])\s+(?=[0-9])", code):
        code = re.sub(r"(?<=[0-9.])\s+(?=[0-9])", "", code)
    return code


def convert(source: Path, target: Path) -> None:
    output = []
    current = []
    label = ""

    def flush() -> None:
        nonlocal current, label
        if not current:
            return
        for i, segment in enumerate(current):
            last = i == len(current) - 1
            prefix = (label + " ") if i == 0 and label else ("  &" if i else "")
            output.append(prefix + normalize(segment).rstrip() + ("" if last else " &"))
        current = []
        label = ""

    for raw in source.read_text(errors="replace").splitlines():
        line = raw.expandtabs(8)
        if not line.strip():
            flush()
            output.append("")
            continue
        if line[0] in "cC*!":
            flush()
            output.append("!" + line[1:].rstrip())
            continue
        fixed = line.ljust(72)
        new_label = fixed[:5].strip()
        continuation = fixed[5:6]
        code = fixed[6:72].rstrip()
        if continuation not in (" ", "0"):
            if not current:
                label = new_label
            current.append(code)
        else:
            flush()
            label = new_label
            current = [code]
    flush()
    target.write_text("\n".join(output) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: fixed_to_free.py input.f output.f90")
    convert(Path(sys.argv[1]), Path(sys.argv[2]))
