# Porting notes

## Upstream

- Package: `betafunctions`
- Upstream version: 1.9.0
- Author: Haakon Eidem Haakstad
- Upstream license: CC0
- Upstream snapshot: `upstream/betafunctions-master/`

This is a clean Fortran translation of computational behavior. The original R
source is retained in the distribution for traceability.

## Numerical replacements

Base-R numerical services have been replaced with self-contained Fortran code:

- regularized incomplete beta: continued fraction with symmetry transform;
- beta inverse CDF: safeguarded bisection;
- regularized incomplete gamma / chi-square upper tail: series and continued fraction;
- adaptive integration: 15-point Gauss-Kronrod recursion;
- gamma RNG: Marsaglia-Tsang;
- beta RNG: ratio of independent gamma variates;
- binomial RNG: Bernoulli summation;
- discrete mixture RNGs: inverse cumulative sampling.

The four-parameter beta CDF and moments use analytic transforms rather than
numerically integrating the density as the R package does. This is mathematically
equivalent and generally more accurate.

## Compatibility choices

Several nonstandard upstream conventions are deliberately retained:

1. `pBetaBinom()` and `pcBinom()` use the lower-tail convention `P(X < q)` rather
   than the usual R discrete-CDF convention `P(X <= q)`.
2. The upstream `qBeta.4P(..., lower.tail = FALSE)` computes
   `1 - qbeta(p, alpha, beta)` before rescaling. The Fortran upper-tail branch
   preserves this behavior rather than replacing it with the more conventional
   `qbeta(1-p, ...)`.
3. `dGammaBinom(..., nc = FALSE)` is not normalized as a continuous density.
   `gamma_binomial_pdf` therefore defaults to the same unnormalized kernel;
   pass `normalized=.true.` for the normalized form.
4. Livingston-Lewis and Hanson-Brennan binary and multicategory routines are
   implemented by the same general category-matrix engine. This removes duplicated
   R control flow without changing the underlying probability calculations.

## Corrections / intentional deviations

The R model-fit routines compute and report a conditional degrees-of-freedom value,
but then calculate the chi-square p-value using `n_bins - 4` unconditionally.
The Fortran code uses the reported degrees of freedom for the p-value. This is an
intentional correction.

The Livingston-Lewis model-fit binning code in R mixes noninteger `pbinom` cut
semantics with rounded observed-score bins and is described upstream as an
experimental approximation. The Fortran helper uses coherent integer score bins
before applying the same expected-count merging rule. Distribution fitting and
classification calculations themselves do not depend on this model-fit helper.

No plotting or graphics-state behavior is translated.

## Classification orientation

For multicategory results, `accuracy_matrix(row, column)` is interpreted as

- row: observed category;
- column: true category.

`consistency_matrix(row, column)` is the joint category probability on two
independent administrations conditional on the same latent true score and then
integrated over the fitted true-score distribution.

For each category the one-vs-rest statistics reproduce the R package's
multicategory calculations.

## Validation

The source was compiled with GNU Fortran using:

```text
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all
```

Tests include:

- incomplete beta, inverse beta, binomial CDF, and chi-square tail references;
- four-parameter beta PDF/CDF/quantile references;
- beta-binomial closed-form reference when `l=0, u=1`;
- normalized Gamma-Binomial numerical-integration references;
- compound-binomial normalization;
- two- and four-parameter moment fitting;
- independent Livingston-Lewis accuracy and consistency matrix references;
- independent Hanson-Brennan matrix references;
- Cronbach alpha and McDonald omega sanity checks;
- RNG support and sample-mean checks.
