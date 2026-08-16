# Porting notes

## Architecture

The R package combines probability functions with `gamlss.family` closures.
The Fortran port separates these concerns:

- `discretedists_distributions` contains the standalone probability laws.
- `discretedists_estimators` contains exported likelihood/start-value helpers.
- `discretedists_family` provides a Fortran derived-type analogue of the
  computational parts of the GAMLSS family objects.
- `discretedists_numerics` supplies local summation, link, Lambert-W(-1), and
  Nelder-Mead helpers.

No R, Rcpp, `gamlss`, `gamlss.dist`, `pracma`, or `nleqslv` runtime is needed.
COMPO/COMPO2 use the vendored `COMPoissonReg-fortran` dependency.

## Intentional compatibility corrections

### DLD CDF indexing

The supplied R `pDLD` returns zero at `q=0` even though `dDLD(0)` is positive.
Direct summation of the PMF shows that the formula in `pDLD` is shifted by one
support point.  The Fortran CDF uses the mathematically consistent `q+1`
indexing, and a regression test checks CDF = cumulative PMF.

### Quantile probability semantics

The Fortran quantile routines use conventional `lower_tail` and `log_p`
semantics: `log_p` applies to the input probability.  This corrects upstream
branches such as `qDIKUM` that apply logarithms to the resulting quantile.

### Closed-form CDFs

DGEII and GGEO PMFs telescope.  Their Fortran CDFs use the resulting exact
closed forms instead of repeated summation.  Quantiles and RNGs therefore avoid
unnecessary accumulation error.

## GAMLSS family representation

`discrete_family_t` implements the numerical behavior needed by a Fortran
caller: link/inverse-link/derivative, PMF/CDF dispatch, score and curvature
contributions, deviance increments, parameter validity, initialization,
moments and randomized quantile residuals.  R formula parsing, environments,
S3 classes and plotting are not reproduced.

For families whose upstream code uses numerical differentiation, the Fortran
family object likewise uses stable central finite differences.  Where upstream
manual score formulas are available, they are retained.

Custom GAMLSS `"own"` link closures cannot be represented by an R-style string;
a Fortran application needing a custom link should wrap/extend the descriptor.

## COMPoissonReg numerical tolerance

The vendored CMP implementation uses its documented hybrid/truncation
normalizer.  Consequently a `nu=1` CMP probability can differ from the exact
Poisson value by the normalizer truncation tolerance (about 1e-7 in the tested
case), even though CMP converges algebraically to Poisson.

## RNG

Random streams use Fortran `random_number`, so generated sequences are not
intended to match R's RNG stream bit-for-bit.  Distributional behavior is
checked instead.
