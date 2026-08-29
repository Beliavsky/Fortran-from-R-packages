# TruncatedNormal-fortran

A modern free-form Fortran/FPM translation of the computational core of the R
package **TruncatedNormal 2.3**.  It implements univariate and multivariate
truncated Normal and Student-t probability estimation and exact simulation,
including the minimax-tilting algorithms of Botev (2017) and Botev &
L'Ecuyer (2015).

## Included computational functionality

The umbrella module is:

```fortran
use truncated_normal
```

The translated API includes:

- robust `lnNpr`, `phinv`, `norminvp`, and `trandn` kernels;
- vectorized univariate `qtnorm_vec` and `rtnorm_vec`;
- both upstream Cholesky permutation strategies through
  `cholperm(..., method='GGE')` and `cholperm(..., method='GB')`;
- multivariate Normal and Student-t log/densities: `dmvnorm`, `dmvt`;
- Monte Carlo probability estimators: `mvncdf`, `mvtcdf`;
- randomized QMC estimators: `mvnqmc`, `mvtqmc`;
- lower-level importance estimators `mvnpr`, `mvnprqmc`, `mvtpr`, `mvtprqmc`;
- exact minimax-tilting generators `mvrandn`, `mvrandt`;
- R-style computational wrappers `pmvnorm`, `pmvt`, `dtmvnorm`, `dtmvt`,
  `ptmvnorm`, `ptmvt`, `rtmvnorm`, and `rtmvt`;
- `tregress`, returning the accepted Student scale variable and Gaussian
  component so that `sqrt(df)*Z/R` has the required truncated Student law;
- conditional sampling for effectively fixed coordinates, matching the
  upstream high-level `rtmvnorm`/`rtmvt` behavior.

Fortran matrices of observations and generated samples use shape `(n,d)`, so
one sample occupies one row.  `tregress_result%z` follows the same convention.

## Dependencies

The attached translated dependencies are vendored under `dependencies/` and
referenced by path from `fpm.toml`:

- `nleqslv-fortran` for the nonlinear minimax-tilting saddlepoint equations;
- `alabama-fortran` as the constrained fallback used by the Normal solver;
- `spacefillr-fortran` for Owen-scrambled Sobol point sets in the Normal QMC
  estimator;
- `qrng-fortran` for randomized Sobol points in the Student-t QMC estimator;
- the MIT-licensed `r_mod.F90` supplied with the earlier work, reused through
  the qrng dependency for Normal/Student helpers, RNGs, and linear solves.

BLAS and LAPACK are linked by the FPM manifest.  R, Rcpp, and RcppArmadillo are
not required.

## Build

With FPM and BLAS/LAPACK installed:

```text
fpm build
fpm test
fpm run --example basic_usage
```

The translation uses standard free-form line lengths and does not require
`-ffree-line-length-none`.

## Tests retained in this port

`test/` contains tests for:

- exact one-dimensional Normal and Student probabilities;
- extreme-tail `lnNpr` stability;
- exact untruncated density normalization;
- both GGE and GB Cholesky permutation paths;
- the independent bivariate Normal orthant probability `1/4`;
- a nontrivial bivariate Normal rectangle probability;
- a nontrivial bivariate Student-t rectangle probability;
- the upstream 15-dimensional equicorrelation identity: for pairwise
  correlation 0.5, the positive orthant probability is `1/(d+1)=1/16`;
- bounded Normal and Student exact sampling;
- fixed-coordinate conditional Normal and Student sampling;
- reconstruction of truncated Student draws from `tregress`.

The fixed two-dimensional reference probabilities used in the tests were also
checked independently with SciPy 1.17 during validation; SciPy is not a build
or runtime dependency.

## Scope

The R package's plotting/documentation/data plumbing, R argument recycling,
S3 presentation, dataset declarations, Rcpp registration, and other
non-computational interface machinery are intentionally omitted.  See
`PORTING_NOTES.md` for the detailed API mapping and implementation decisions.

## License

TruncatedNormal-derived Fortran code is GPL-3.0-only, matching upstream
`License: GPL-3`.  Vendored dependency ports retain their own license files
and notices.  The reused `r_mod.F90` remains MIT-licensed.  See `NOTICE.md`,
`LICENSE-GPL-3.txt`, and `dependencies/*/LICENSE*`.
