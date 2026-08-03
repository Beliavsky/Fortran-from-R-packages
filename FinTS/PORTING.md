# Porting notes

## General representation

R vectors, matrices, lists, S3 objects, and data frames are represented by Fortran arrays and derived types. Dots in R names are written as underscores where necessary, for example `FinTS.stats` becomes `FinTS_stats` and `as.yearmon2` becomes `as_yearmon2`.

Plotting and display code was not translated. `plotArmaTrueacf` remains available as a compatibility name, but it only computes and returns roots, ACF/PACF values, damping, and periods.

## ARIMA

The original `ARIMA` function is mostly a wrapper around `stats::arima`. Reproducing R's exact state-space likelihood, diffuse initialization, missing-value processing, coefficient covariance matrix, and optimizer behavior would require a much larger runtime.

The Fortran implementation instead uses:

1. regression removal on the original series;
2. seasonal and ordinary differencing;
3. multiplicative seasonal/nonseasonal AR and MA polynomials;
4. conditional innovations recursion;
5. profiled conditional Gaussian likelihood;
6. Nelder-Mead optimization;
7. optional reflection-coefficient transformations for stationarity and invertibility.

`CSS`, `ML`, and `CSS-ML` are accepted for source-level familiarity, but all currently select this conditional Gaussian implementation. The method name is retained in the result.

The R arguments `fixed`, `optim.control`, and `kappa` are not represented. `initial` is supported, but its AR and MA entries are on the unconstrained optimizer scale when `transform_pars=.true.`.

For differenced models, residuals are aligned to the final `n-d-D*period` observations. Earlier positions are IEEE NaNs. Fitted values for the aligned positions are reported as `x-residual`; they are diagnostic alignments rather than a full inverse-differenced state-space reconstruction.

## ACF and tests

The R `acf` implementation handles several missing-data policies. The translated `acf`, `cross_acf`, `AutocorTest`, and `ArchTest` require finite observations, matching the original default `na.fail` behavior for `Acf`.

Chi-square p-values use an internal regularized incomplete-gamma implementation.

## APCA

The APCA implementation follows the two eigenvalue/regression passes in the R source. Symmetric eigendecomposition uses an internal Jacobi method. Factor signs can differ from R because eigenvector signs are not unique.

## ARMA roots and theoretical correlations

Polynomial roots use a Durand-Kerner iteration. Theoretical ARMA autocovariances use an adaptively truncated infinite-MA representation. This is accurate for ordinary stationary models but may converge more slowly near a unit root.

The original function obtains polynomial roots through the optional R package `polynom`; this port has no external dependency.

## Summary moments

The skewness and excess kurtosis formulas match the normalized central-moment convention used by the default `e1071` routines referenced by FinTS. Standard deviation uses the sample denominator `n-1`.

## Compound interest ambiguity

The original R source computes a conditional net-value expression but then discards it and returns the gross value unconditionally. Its documentation and examples also disagree about the meaning of `net.value`.

This port uses the apparent intended executable convention suggested by the example and source expression:

- omitted or `.false.`: gross value;
- `.true.`: gross value minus one.

This difference is documented rather than silently reproducing the discarded-expression bug.

## Dates and data access

`as.yearmon2` is represented for numeric `yyyy.mm` and integer `yyyymm` inputs. Arbitrary R date classes, format strings, names attributes, and `zoo` indices do not have direct Fortran equivalents.

`read.yearmon`, `url2data`, `package.dir`, and `runscript` are data-access or R-runtime utilities and are not part of the computational library.
