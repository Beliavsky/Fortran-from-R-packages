# Attribution and license notice

This project is a modern Fortran/FPM translation of the computational code in
the R package **tweedie 3.1.0**, authored and maintained by Peter K. Dunn.
The upstream package declares `License: GPL (>= 2)`, so source derived from it
is distributed as **GPL-2.0-or-later**.

The Fourier-inversion implementation in `src/00tweedie_params.f90` through
`src/twcomputation_loop.f90` comes directly from the upstream package's native
Fortran source. Those numerical files are retained unchanged; only the R
printing/interface shims around them were replaced with standalone Fortran
modules.

The density algorithms should be cited as appropriate using the upstream
references:

- Peter K. Dunn and Gordon K. Smyth (2005), "Series evaluation of Tweedie
  exponential dispersion models", *Statistics and Computing* 15(4), 267-280.
- Peter K. Dunn and Gordon K. Smyth (2008), "Evaluation of Tweedie exponential
  dispersion models using Fourier inversion", *Statistics and Computing*
  18(1), 73-86.

The original `DESCRIPTION`, `NAMESPACE`, and `CITATION` are retained under
`upstream/` for provenance.

`src/r_mod.F90` is the separately supplied R-compatibility helper module. It
remains under its stated **MIT License**. The build copy differs from the
supplied `upstream/r_mod-original.f90` only by free-form line wrapping needed
to keep source lines within the standard 132-column limit; a normalized
lexical comparison is identical.
