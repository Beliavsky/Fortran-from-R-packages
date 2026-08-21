# Translation status

Source package: `bzinb` 1.0.8.

## Exported R API mapping

| R routine | Fortran API | Status |
|---|---|---|
| `bp` | `fit_bp` | Complete numerical equivalent |
| `rbp` | `rbp_sample` | Complete |
| `lik.bp` | `loglik_bp` | Complete |
| `bzip.a` | `fit_bzip_a` | Upstream EM translated |
| `rbzip.a` | `rbzip_a_sample` | Complete |
| `lik.bzip.a` | `loglik_bzip_a` | Complete |
| `bzip.b` | `fit_bzip_b` | Upstream EM translated |
| `rbzip.b` | `rbzip_b_sample` | Complete |
| `lik.bzip.b` | `loglik_bzip_b` | Complete |
| `bnb` | `fit_bnb`, `fit_bnb_em` | Specialized upstream BNB EM translated; direct optimizer retained as `fit_bnb_direct` |
| `rbnb` | `rbnb_sample` | Complete |
| `lik.bnb` | `loglik_bnb` | Complete |
| `bzinb` | `fit_bzinb`, `fit_bzinb_em` | Specialized upstream BZINB EM translated; direct optimizer retained as `fit_bzinb_direct` |
| `rbzinb` | `rbzinb_sample` | Complete |
| `lik.bzinb` | `loglik_bzinb` | Complete |
| `bzinb.se` | `bzinb_standard_errors` | Upstream score outer-product information translated |
| `idigamma` | `idigamma` | Source-compatible inverse-digamma Newton iteration |
| `weighted.pc` | `weighted_pearson_correlation` | Complete |
| `pairwise.bzinb` | `pairwise_bzinb_full` | Full numerical result; `pairwise_bzinb` retained as compact rho/SE matrix interface |

The non-exported density functions (`dbp`, `dbzip.a`, `dbzip.b`, `dbnb`,
`dbzinb`) remain public Fortran PMF/log-PMF routines because they are useful
building blocks in native applications.

## Specialized BNB/BZINB EM parity

Version 0.2.0 translates the three C++ sources that were the main v0.1 parity
gap:

- `src/expt.cpp`: latent expectation sums, overflow/underflow rescaling, score
  contributions, and score outer-product information;
- `src/opt.cpp`: inverse-digamma updates and the coupled Newton update for
  `log(b1)`;
- `src/em.cpp`: mixture-probability M-step, `b2` latent-moment update,
  historical-maximum selection, convergence tolerance, and the 100-iteration
  no-improvement safeguard.

The upstream implementation parameterizes the information matrix by
`(a0,a1,a2,b1,b2,p1,p2,p3)`, with `p4=1-p1-p2-p3`. The Fortran fit result stores
that exact 8 x 8 information matrix and additionally expands its inverse to a
9 x 9 covariance matrix by the linear `p4` constraint. Thus the first 8 x 8
covariance block agrees with the upstream parameterization while the ninth row
and column are a native convenience extension.

The C++ mixture-score formulas have a somewhat unusual convention (for example,
for a positive-positive observation the sixth score component is 1, rather
than the ordinary derivative `1/p1`). The Fortran information engine preserves
that source convention exactly instead of replacing it with a textbook score.

The C++ code contains an undefined-behavior corner case if all three
inverse-digamma targets enter its `>=600` approximation branch before the local
`idgam` array has ever been initialized. The Fortran implementation uses the
same large-value approximation but supplies finite shape values rather than
emulating uninitialized memory. This is a safety fix for an otherwise
pathological numerical branch.

## Initial values and edge cases

The default BNB/BZINB moment/profile initial values follow the R code. The R
`bnb()` wrapper appears to overwrite a user-supplied `initial` value while
constructing its default start. The Fortran API intentionally honors an explicit
`initial`, which is both useful and consistent with the documented R argument.

The upstream all-zero shortcut is represented by a typed Fortran result:
positive parameters are `1e-10`, BZINB mixture probabilities are
`(1,0,0,0)`, log likelihood is zero, and the fit is returned immediately.

## Pairwise API

`pairwise_bzinb_full` supplies the numerical content of `pairwise.bzinb`:

- feature-pair indices;
- rho and its standard error;
- optional full BZINB parameter and standard-error arrays;
- log likelihood and best EM iteration;
- convergence flag;
- nonzero proportions and their pairwise minimum;
- optional random subsampling of feature pairs.

R data-frame labels and formatted `pair` strings are intentionally not part of
the Fortran result type.

## Excluded infrastructure

- Rcpp registration and Boost wrappers;
- R `Vectorize`, data-frame/table construction, names/dimnames;
- R `try`, warning, and printing semantics;
- R console formatting and documentation-only examples.

There is no remaining intentionally substituted statistical fitting algorithm in
the exported BNB/BZINB path. Remaining differences are interface/runtime
semantics and the pathological uninitialized-memory C++ branch noted above.
