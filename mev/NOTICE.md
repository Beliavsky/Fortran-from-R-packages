# NOTICE

This is a modern Fortran translation of the computational code in the R package
`mev` version 2.2 (Modelling of Extreme Values).

Upstream declares `License: GPL-3`. This translated work preserves that license
as GPL-3.0-only. Original package metadata and source files used for translation
are retained under `orig/` for provenance.

Version 0.3.0 retains two user-supplied compatible Fortran dependency
translations:

- `nleqslv-fortran` 0.1.1, GPL-2.0-or-later;
- `expint-fortran` 0.1.0, GPL-3.0-or-later.

Their license and attribution files are retained under their corresponding
`vendor/` directories. The separately supplied GPL-2.0-only Rsolnp and mvtnorm
translations are not included or linked; see `PORTING_NOTES.md`.

Plotting, S3 presentation, R formula/model-frame code and packaged datasets are
not part of the Fortran API.
