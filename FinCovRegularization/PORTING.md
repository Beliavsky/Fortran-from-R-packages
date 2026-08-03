# Porting notes

## Scope

All 16 exported numerical routines in FinCovRegularization 1.1.0 were
translated. The following R-only presentation infrastructure was not ported:

- `plot.CovCv`
- `print.CovCv`
- `summary.CovCv`

The original `.rda` sample dataset is preserved under `original/data/`, but no
R-data parser is bundled. The Fortran examples use deterministic synthetic
returns.

## Numerical implementation

The translation is self-contained and does not require BLAS, LAPACK, or a
quadratic-programming package. Dense Gaussian elimination with partial
pivoting, a symmetric Jacobi eigensolver, active-set long-only GMVP, and
Nelder-Mead risk-parity optimization are implemented in the source tree.

Cross-validation uses a portable Park-Miller generator rather than the R random
number generator. Therefore the same seed is reproducible within this Fortran
port, but it does not select the same rows as R's `set.seed` and `sample`.

## Intentional corrections and clarifications

The multi-factor branch of the original `MacroFactor.Cov` contains the
expression

```text
beta[-1,] %*% cov(factor) %*% t(beta[-1,])
```

which is dimensionally invalid for the documented `(n,q)` factor matrix when
`q > 1`. The Fortran translation uses the standard and dimensionally correct
form

```text
transpose(beta_slopes) * covariance_factors * beta_slopes
```

The one-factor result is unchanged.

The original `Ind.Cov` does not simply retain the diagonal of its input. It
first applies `cov(sigma)` and then retains that covariance matrix's diagonal.
The Fortran routine preserves this implemented behavior.

`StatFactor.Cov` calls singular values "eigenvalues" and applies the Kaiser
cutoff directly to singular values. The Fortran implementation preserves that
selection rule while obtaining equivalent singular values from the sample
covariance eigenvalues.

The R portfolio routines round returned weights to four digits. The Fortran
routines retain full precision; callers can round for presentation.
