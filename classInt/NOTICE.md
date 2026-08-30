# Notices and provenance

## classInt

This project is a modern Fortran translation of the computational functionality
of the R package `classInt` 0.4-11 by Roger Bivand and contributors.

Upstream project:

- https://github.com/r-spatial/classInt
- https://r-spatial.github.io/classInt/

The upstream package is licensed GPL (>= 2).  A complete source snapshot used
for this translation is retained in `upstream/`.

The Fisher exact classification routine in this translation is a modernization
of the algorithm embodied by upstream `src/fish1.f`; the original file is
retained unchanged in the upstream snapshot for provenance and comparison.

## KernSmooth direct plug-in bandwidth algorithm

The implementation of style `dpih` is based on the computational algorithm in
`KernSmooth::dpih`, including linear binning, recursive binned kernel functional
estimation, and the power-of-two FFT convolution used by `bkfe`.  The
KernSmooth source states that its code may be used and distributed without
restriction.  The Fortran transform is a native radix-2 implementation rather
than a call into R's FFT runtime.

KernSmooth project source:

- https://cran.r-project.org/package=KernSmooth

## R pretty-break behavior

The native `pretty_breaks` routine follows the default numerical break-selection
logic of R's `pretty` implementation for ordinary finite ranges.  R is
GPL-2.0-or-later; this package's GPL-2.0-or-later distribution is compatible
with that provenance.

## e1071-fortran and proxy-fortran

The `bclust` classification style uses the previously translated
`e1071-fortran` package, which is vendored under `dependencies/e1071-fortran`.
Its own notices, licenses, upstream source, and the vendored `proxy-fortran`
dependency are retained in their respective directories.

## Translation

The modern Fortran translation was generated and reviewed as a source-level
port.  It is not an official release of the upstream R package.
