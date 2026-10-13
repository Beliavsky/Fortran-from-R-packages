# Upstream provenance and dependency review

- Upstream package: ineq 0.2-13, dated 2014-07-21.
- Official archive: https://cran.r-project.org/src/contrib/ineq_0.2-13.tar.gz
- Archive SHA-256: `e0876403f59a3dfc2ea7ffc0d965416e1ecfdecf154e5856e5f54800b3efda25`.
- Original license expression: `GPL-2 | GPL-3`; equivalent SPDX choice used here:
  `GPL-2.0-only OR GPL-3.0-only`.
- Upstream authorship: Achim Zeileis (author, maintainer), Christian Kleiber
  (contributor), as recorded verbatim in `upstream/DESCRIPTION` and `upstream/THANKS`.
- Translation target: https://github.com/Beliavsky/Fortran-from-R-packages .
- Repository snapshot inspected for compatible shared and top-level APIs:
  `cbd437d6d92f5481706fd8ea5a7829aaff7b6e24`.
- Tested sibling dependency: `rfortran-core` from that snapshot, unchanged.

## Source retention

`upstream/` is a reference-only subset of the official archive: DESCRIPTION,
NAMESPACE, MD5, THANKS, all six R files, and all manual pages. These files are
byte-for-byte copies. `upstream/MD5` remains the original full-package manifest
and includes omitted data/tests. It must not be interpreted as a manifest of this
translation. `validation/upstream-sha256.json` records each retained file.

R `ineq.R` and `conc.R` map to `src/ineq_indices.f90` (inequality/concentration).
R `pov.R` maps to the same module (poverty).
R `Lc.R` and `major.R` map to `src/ineq_lorenz.f90`.
`src/ineq_common.f90` provides translation-specific validation/scaling support;
`src/ineq_api.f90` re-exports the public numerical interfaces.
Plotting code in `Pen.R` and the plotting/lines methods in `Lc.R` are not translated.

## Dependency compatibility review

The target repository already provides compatible public APIs in rfortran-core:

| Required operation | Reused API |
|---|---|
| Named real kind | `r_kinds: dp` |
| Stable numeric sorting and permutation | `r_ordering: r_sort_values_in_place, r_order` |
| Normal distribution and inverse | `r_distributions: r_pnorm, r_qnorm` |
| Regularized incomplete beta | `r_special: r_regularized_beta` |
| Stable elementary transformations | `r_transforms: r_log1p, r_expm1` |
| Log-mean-exponential | `r_stability: r_log_mean_exp` |

Theoretical Lorenz formulas use the real beta-CDF and normal-CDF/quantile
implementations, not approximate placeholders or numerical-integration stubs.
A small domain adapter preserves R's zero-beta-shape limits and out-of-range CDF
saturation, which differ from r_regularized_beta's strict argument domain.

No new stats, sorting, normal-quantile, beta-integration, BLAS/LAPACK, or general
R-compatibility implementation is embedded. No rfortran-compat dependency is used.
The sibling path dependency is intentional: FPM resolves the supplied repository
layout. To reproduce the tested shared sources, use the commit above; a path
manifest cannot itself enforce that Git revision.
