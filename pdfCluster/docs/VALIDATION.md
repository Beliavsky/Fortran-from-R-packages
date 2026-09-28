# Validation

This file records release validation for the Fortran translation.

The deterministic test suite contains:

- `test/unit/main.f90`: bandwidth rules, KDE, valley measure, one-dimensional
  graph construction, connected components, and adjusted Rand index including
  the degenerate NaN case.
- `test/clustering/main.f90`: one-dimensional two-cluster fitting and staged
  classification, density-based silhouette, pair-graph refitting, and the
  Delaunay path through the sibling `geometry` package.
- `example/basic.f90`: a minimal end-to-end clustering example.

Release validation also checks coverage metadata, public procedure existence,
free-form/ASCII source, 132-column line limits, dummy-argument `INTENT`/`VALUE`
and FORD comments, the package-wide real kind, real literal suffixes, prohibited
semicolon-separated statements, unsafe NaN self-comparisons, forbidden
BLAS/LAPACK links, duplicate sources, copied dependency source, and build
artifacts.

The exact compiler/FPM command outcomes for this release are appended during
final packaging.

## Release results

Compiler: GNU Fortran 14.2.0.

Two direct compiler validations were run against the current sibling `geometry`
source tree:

1. Strict checking:
   `-std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace -O0`
2. Optimized:
   `-std=f2018 -Wall -Wextra -Werror -O2`

In both configurations:

- all maintained `geometry` dependency modules required by its public facade
  compiled successfully;
- all `pdfCluster` source modules compiled successfully;
- `test/unit/main.f90` printed `All pdfCluster unit tests passed.`;
- `test/clustering/main.f90` printed `All pdfCluster clustering tests passed.`;
- `example/basic.f90` completed with two detected cores and groups 1/2 on the
  deterministic separated fixture.

The requested literal FPM commands were attempted from the package root:

```text
fpm build
fpm test
fpm run --example basic
fpm clean --all
```

Each returned exit code 127 with `fpm: command not found` because this execution
environment does not contain an FPM executable. These outcomes are an
environment limitation and are not reported as successful FPM validation.
