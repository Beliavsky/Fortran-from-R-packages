# Provenance

## Upstream input

The translation was produced from the user-supplied archive
`mitools-master.zip` containing CRAN `mitools` version 2.4.

SHA-256 of the supplied archive:

```text
775c808fe98a3d1e4b64ed57bfc8b91ad49642809f3ed7dde511010123720b73
```

Relevant upstream files and SHA-256 values:

```text
a8c6a98997494929d954ceb7527385279fd1bf2899da5f1cf49016eb24905bb8  DESCRIPTION
82947ddd0044f08e963c4ace017038a8e00548432ba848fb0ba42a12f7e36fd6  NAMESPACE
0fd157776c26a656540f760597829766423f1cc94b706f47c2e217a46a9e8089  R/MI.R
d687bf3b481d2460f1546c53d24f84286f3f1cc4bd02805bfc4fe73ff24e7ea7  R/PV.R
c9ca17ea4efdb2ad95d5b09e10a366aa818551a2edfe8f913b6a16644eb74982  R/update.R
```

Exact reference copies are retained below `reference/upstream/`.

## Translation mapping

- `R/MI.R` -> `src/mitools_combine.f90` for `MIcombine` and numerical summary
  calculations; `src/mitools_imputation.f90` for numeric imputation-list
  storage and binding.
- `R/PV.R` -> `src/mitools_pv.f90` for deterministic plausible-value
  selection/materialization.
- `R/update.R` -> no computational Fortran translation because its routines
  plan SQL queries and evaluate stored R expressions.

The public facade is `src/mitools.f90` and shared result/container types are in
`src/mitools_types.f90`.

## Dependency review

Before implementation, the root of
`Beliavsky/Fortran-from-R-packages` was checked. No existing top-level `mitools`
translation was present. The repository's `rfortran-core` package documents
compatible APIs for the common real kind, covariance, and central Student-t
quantiles. This package therefore uses the sibling path dependency
`../rfortran-core` and does not duplicate those implementations.

No BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, translated R-package dependency,
or other dependency source is included in this directory.
