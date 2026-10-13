"""Convert checked-in, genuine R-generated fixtures into compile-time Fortran constants."""
from collections import defaultdict
from pathlib import Path
base = Path(__file__).resolve().parents[1]
rows = defaultdict(list)
for line in (base / 'reference/expected.tsv').read_text().splitlines():
    name, index, value = line.split('\t')
    rows[name].append(value.strip())
text = ['! SPDX-License-Identifier: GPL-2.0-or-later',
        '! Generated from reference/expected.tsv by reference/embed_reference.py.',
        'module esg_reference_values', '   use r_kinds, only: dp', '   implicit none']
for name, values in rows.items():
    text.append(f'   real(dp), parameter :: ref_{name}({len(values)}) = [ &')
    for index, value in enumerate(values):
        text.append('      ' + value + '_dp' + (', &' if index < len(values)-1 else ']'))
text += ['end module', '']
(base / 'test/esg_reference_values.f90').write_text('\n'.join(text))
