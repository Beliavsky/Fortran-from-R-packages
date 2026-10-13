# Copyright, attribution, and provenance

This is a modern Fortran translation and adaptation of the computational code
of R package **randtests 1.0.2**, dated 2024-04-22, published on CRAN 2024-04-23.

Upstream authors: **Frederico Caeiro** and **Ayana Mateus**.
Upstream maintainer: Frederico Caeiro <fac@fct.unl.pt>.
Upstream license declaration: **GPL (>= 2)**. The translated/adapted code is
provided under **GPL-2.0-or-later**. Full GPL version 2 and version 3 texts are
included in LICENSE and LICENSE-GPL-3. No narrower copyright ownership claim is
made for the upstream material. The source files contain no additional explicit
copyright lines; all upstream notices have been retained unmodified.

Source:
https://github.com/cran/randtests/tree/7244d86764445e657634c9ae4d59ce942a5fcbc8

CRAN project: https://cran.r-project.org/package=randtests

Source revision: `7244d86764445e657634c9ae4d59ce942a5fcbc8`.
Destination repository baseline used for dependency verification:
`cbd437d6d92f5481706fd8ea5a7829aaff7b6e24`.
Translation date: 2026-10-11.

The unmodified DESCRIPTION, NAMESPACE, Changes, R sources, and man pages are
retained under upstream/. The sweetpotato dataset is not translated or included.
The upstream MD5 index is omitted because this distribution is a selected source
snapshot rather than a CRAN installation tarball.

## Earlier attribution

Upstream R/misc.R and man/permut.Rd state that `permut` is based on the
`permutations` function from the **gtools** R package. This provenance is
preserved verbatim in those files and acknowledged in the translated source.
No gtools package source is copied, vendored, or required as a dependency.

## Shared dependency

`rfortran-core` is reused via `path = "../rfortran-core"`. It supplies `dp`, normal
probabilities, regularized beta, log binomial coefficients, stable log sums,
ordering, and medians. It is MIT-licensed with author attribution to Vivek Rao
and contributors in its own distribution. No dependency sources, R runtime,
BLAS/LAPACK, or translated third-party dependency packages are embedded here.

## Modifications

The Fortran implementation uses typed results, explicit status codes, pure and
elemental interfaces, finite-value validation, bounded exact enumeration, and
Fortran's random number generator. Deliberate corrections and compatibility
limits are detailed in README.md, API_COVERAGE.md, and the manifest mappings.
These include exact Bartels statistics/tails, log runs CDFs, quantile tail flags,
small-permutation callbacks, and p-values capped at one. Tests distinguish R
parity cases from mathematically corrected cases.
