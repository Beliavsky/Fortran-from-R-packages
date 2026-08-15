# Licensing

This repository is a source translation of the computational code in the R package `distr` 2.9.7.

## distr-derived Fortran code

The upstream package declares `License: LGPL-3`. The translated modules

- `src/distr_kinds.f90`
- `src/distr_special.f90`
- `src/distr_fft.f90`
- `src/distr_rng.f90`
- `src/distr_core.f90`
- `src/distr_matrix.f90`
- `src/distr.f90`

are distributed under **LGPL-3.0-only**. The full LGPL v3 text is in `COPYING.LESSER`; the GPL v3 text referenced by the LGPL is in `COPYING.GPL-3`.

Upstream authors/copyright holders listed by `distr` 2.9.7 include Matthias Kohl and Peter Ruckdeschel, with contributions from Florian Camphausen and Thomas Stabla. See `upstream/DESCRIPTION` and `upstream/CITATION` for the authoritative upstream metadata and citation requests.

## Kolmogorov-Smirnov code

`src/distr_ks.f90` is a translation of `distr/src/ks.c`, which states that it was taken from R Core `stats/src/ks.c` and is copyright the R Core Team. The original file permits redistribution/modification under **GPL-2.0-or-later**. That module therefore remains separately GPL-2.0-or-later.

The original C file is retained verbatim as `upstream/ks.c`. The GPL v2 and GPL v3 license texts are in `COPYING.GPL-2` and `COPYING.GPL-3`.

The umbrella module `distr` intentionally does **not** import `distr_ks`; applications that need the KS routines should explicitly `use distr_ks`. This keeps the licensing boundary visible to downstream users.

## QQ confidence-band code

`src/distr_qq.f90` translates the non-graphical calculations from `qqbounds.R` and `internals-qqplot.R`. Because the module imports `distr_ks`, it is distributed under **GPL-3.0-or-later**. It is also kept outside the LGPL umbrella module and must be imported explicitly with `use distr_qq`. The relevant upstream R sources are retained under `upstream/`.

No license change or relicensing of upstream material is intended by this translation.
