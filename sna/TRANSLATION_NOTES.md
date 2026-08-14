# Translation notes

The source distribution in `upstream/` identifies the original package as
`sna` 2.8 (2024-09-07), authored by Carter T. Butts and licensed GPL-2 or
later.  `COPYING` is copied from the upstream source archive.

The original package combines R orchestration with C kernels.  The Fortran
port consolidates the computational code into typed modules rather than
mirroring the R/C boundary.  Adjacency matrices are `real(dp)`, a stack is
`real(dp) :: graphs(m,n,n)`, and missing numeric values are represented by
quiet IEEE NaNs.

Plotting and graphical layout code was intentionally skipped per the
translation request.  R S3 methods and conversion helpers were also omitted
when they only serve R object infrastructure.

The specialized 64-state biased-net triadic pseudolikelihood was translated
from the upstream C kernel and checked against that implementation.  L-NAM
standard errors use a self-contained finite-difference Hessian.

See the README for the routine mapping and explicit API adaptations.
