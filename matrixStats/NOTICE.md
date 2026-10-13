# Provenance, modifications and licensing

This package is a modified modern free-form Fortran translation of **matrixStats
1.5.0**, by Henrik Bengtsson and the contributors listed in upstream/DESCRIPTION.
It is not an official upstream release and makes no endorsement claim.

## Pinned sources

- matrixStats CRAN mirror: https://github.com/cran/matrixStats
- Revision: d229f6365e7b5ad0351c4ce98948304e8d23fcae
- Upstream project: https://github.com/HenrikBengtsson/matrixStats
- Package version: 1.5.0, CRAN publication 2025-01-07
- Target repository: https://github.com/Beliavsky/Fortran-from-R-packages
- Tested sibling core revision: cbd437d6d92f5481706fd8ea5a7829aaff7b6e24
- R oracle version: 4.5.0, using the above matrixStats source revision

The exact computational R wrapper files are retained under upstream/R, together
with upstream/DESCRIPTION and NAMESPACE. They are provenance references, not
runtime dependencies. UPSTREAM-NOTICES.txt retains original C implementation
copyright and attribution comment blocks used during translation. No original
C code, r.f90/r_mod.f90, BLAS/LAPACK/ARPACK source, or dependency package is embedded
in the maintained Fortran source. No source from rfortran-core is copied here.

## Authors and preserved notices

matrixStats copyright holder and primary author: Henrik Bengtsson.
Contributors: Constantin Ahlmann-Eltze, Hector Corrada Bravo, Robert Gentleman,
Jan Gleixner, Peter Hickey, Ola Hossjer, Harris Jaffee, Dongcan Jiang,
Peter Langfelder, Brian Montgomery, Angelina Panagopoulou, Hugh Parsonage,
and Jakob Peder Pettersen. Per-routine historical years and attribution are
preserved verbatim in UPSTREAM-NOTICES.txt and the original R sources.

The private R-compatible quicksort in src/matrixstats_order.f90 is adapted from:
https://github.com/wch/r-source/blob/b4322f3207e8d99d9dc6e6fa6b80ad23c4dc2f5f/src/main/qsort-body.c
Source SHA-256: 5fa898a4bf0c94a261256d4427fd7380336ea939145d6a0c847efe5dd8651b65.
Copyright (C) 2002-2012 The R Core Team. GPL-2.0-or-later. The source algorithm
is R.C. Singleton, Algorithm 347, Communications of the ACM 12(3), 1969.
Its exact unstable equal-key ordering materially affects matrixStats weighted
median/MAD interpolation. Shared r_order is stable and therefore incompatible
for this one operation; all compatible sorting operations reuse rfortran-core.
The original R copyright/license notice is retained beside the private helper.

matrixStats-derived work remains Artistic-2.0. The combined distribution also
contains this GPL-2.0-or-later adaptation, so both licenses apply. Full texts:
LICENSE-Artistic-2.0 and LICENSE-GPL-2.0. The sibling rfortran-core dependency
has its own MIT license, which is not replaced by this package's license.

## Translation changes

2026-10-11: introduced independent pure Fortran modules, allocatable numerical
interfaces, explicit double-precision kind imported from the shared core,
R-style numeric selections, deterministic numerical tests and examples.
R runtime allocation entry points are excluded; R attributes, names, type
coercion, callbacks, random rank tie behavior and exception/warning objects
are not reproduced. NA/NaN use one IEEE NaN representation. Some compensated
and rescaled dp arithmetic improves extreme finite overflow behavior; this is
not a bitwise-identical R implementation. Full mapping limits are documented
in README.md, API_COVERAGE.md and the FPM translation metadata.

## Explicit corrections of upstream defects

The following small, intentional differences are tested and disclosed rather
than silently reproducing upstream implementation defects:

1. With an odd number of infinite weights, upstream weightedMedian can read the
   wrong partially sorted element. For `x=[2,0,2], w=[Inf,Inf,Inf]`, matrixStats
   1.5.0 returns 0; this translation returns the mathematically correct median 2.
   The same repair applies to row/column weighted medians. Finite-weight tie
   ordering still follows R's exact unstable quicksort behavior.
2. Scalar weightedMad with no explicit weights ignores upstream `idxs` entirely;
   weightedVar/weightedSd with no explicit weights retain the wrong original
   weight-vector length after subsetting. Here all three consistently select
   observations first and then generate equal weights. For
   `x=[4,1,7,1,9,2], idxs=[1,2,6]`, this gives MAD=1.4826 and variance=7/3,
   versus upstream 2.9652 and 1.8666666667. Matrix wrappers are unaffected because
   their selection occurs before invoking vector summaries.
3. Compensated/scaled dp arithmetic recovers some finite results that upstream's
   specific accumulation order overflows. Mean residual corrections are skipped
   if they are nonfinite. These changes do not collapse actual NaNs or infinities.

