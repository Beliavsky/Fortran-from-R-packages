"""Embed the checked-in real R output as portable Fortran test constants."""
from pathlib import Path

root = Path(__file__).resolve().parent.parent
lines = [
    "! Generated from actual CRAN CreditRisk 0.1.7 output by reference/generate_reference.R.",
    "module creditrisk_reference",
    "  use r_kinds, only: dp",
    "  implicit none",
    "  private",
]
for path in sorted((root / "reference").glob("*.txt")):
    if path.name.startswith("R-"):
        continue
    rows = [[float(value) for value in row.split()] for row in path.read_text().splitlines()]
    if not rows:
        continue
    name = "ref_" + path.stem
    nr, nc = len(rows), len(rows[0])
    assert all(len(row) == nc for row in rows)
    values = [rows[i][j] for j in range(nc) for i in range(nr)]
    lines.append(f"  real(dp), parameter, public :: {name}({nr},{nc}) = reshape([ &")
    for i, value in enumerate(values):
        ending = ", &" if i < len(values) - 1 else f" ], [{nr},{nc}])"
        lines.append(f"    {value:.17e}_dp" + ending)
lines.append("end module creditrisk_reference")
(root / "test/reference_values.f90").write_text("\n".join(lines) + "\n")
