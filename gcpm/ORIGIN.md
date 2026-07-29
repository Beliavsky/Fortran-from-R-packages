# Origin and licensing

## Original package

- Package: GCPM
- Version: 1.2.2
- Title: Generalized Credit Portfolio Model
- Date: 2016-12-29
- Author and maintainer listed in `DESCRIPTION`: Kevin Jakob
- Startup copyright notice: Copyright (C) 2015 Kevin Jakob and Dr. Matthias
  Fischer
- Original license field: GPL-2

A copy of the original package `DESCRIPTION` and `MD5` metadata is retained in
`original/`. The original example portfolio files are retained in `data/`.

## Translation license

The source package's startup notice grants redistribution and modification
under GNU General Public License version 2. The translated source files carry:

```text
SPDX-License-Identifier: GPL-2.0-only
```

The full GNU GPL version 2 text is in `LICENSE`.

## Translation scope

This project translates the numerical methods and data flow into modern
Fortran. It does not reproduce R's S4 object system, graphics, progress bars,
R serialization, or R-specific file-export interface. The original C++ loss
simulation and the computational R methods were used as the source material
for the corresponding Fortran routines.
