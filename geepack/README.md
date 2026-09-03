# geepack

Modern free-form Fortran translation of the computational code in R package
`geepack` 1.3.13.

The library implements generalized estimating equations for clustered and
longitudinal responses, including separate mean, scale, and association
estimating equations; ordinary and ordinal GEE; multiple working-correlation
structures; sandwich and jackknife covariance estimates; QIC/CIC; and relative
risk regression by the COPY method.

## Build

Place this directory at the repository root beside `rfortran-core` and
`rfortran-linalg`, then run:

```text
fpm build
fpm test
fpm run --example geepack_example
```

`fpm.toml` uses sibling path dependencies:

```toml
rfortran-core = { path = "../rfortran-core" }
rfortran-linalg = { path = "../rfortran-linalg" }
```

No system BLAS or LAPACK library is required by `geepack`; linear algebra is
routed through the repository's `rfortran-linalg` package.

## Public API

Use the facade module:

```fortran
use geepack
```

Important entry points are:

- `fit_geese` — ordinary GEE for mean/scale/association models
- `fit_ordgee` — clustered ordinal GEE with local odds-ratio association
- `gen_zcor`, `gen_zodds`, `fixed_to_zcor` — association-design helpers
- `compute_qic`, `quasi_likelihood` — QIC/CIC calculations
- `compare_coefficients` — influence-based coefficient comparison
- `make_relative_risk_copy`, `fit_relative_risk` — COPY relative-risk fit
- `coefficient_wald_summary`, `wald_contrast` — numerical inference helpers

`gee_spec` controls links, variance functions, scale handling, working
correlation, convergence tolerance, and AJS/J1S/FIJ jackknife calculations.

## Numerical codes

Working correlations:

- `COR_INDEPENDENCE`
- `COR_EXCHANGEABLE`
- `COR_AR1`
- `COR_UNSTRUCTURED`
- `COR_USERDEFINED`
- `COR_FIXED`

Mean/scale/association links include identity, logit, probit, cloglog, log,
reciprocal, Fisher-z, and both Lin–Wei–Ying links. The upstream link called
`inverse` is named `LINK_RECIPROCAL` here to avoid colliding with the Fortran
procedure `link_inverse`.

## Input convention

The numerical layer deliberately does not implement R formula or data-frame
semantics. Supply design matrices directly. Rows must be grouped into contiguous
clusters, with `cluster_sizes` giving the row count for every cluster. `waves` selects wave-specific links/variances; optional real-valued `cor_param` supplies the upstream `corp` coordinates used by AR(1) distances.

See `API_COVERAGE.md` for a detailed R-to-Fortran map, `PROVENANCE.md` for
source mapping, and `NOTICE.md`/`LICENSE` for attribution and licensing.
