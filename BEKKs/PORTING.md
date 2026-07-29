# Porting notes

## Source package

This project translates the computational algorithms in BEKKs 1.4.7 from R/RcppArmadillo to modern Fortran. The complete source package is retained under `original/`.

## Directly represented numerical areas

- full, diagonal, and scalar BEKK likelihoods;
- symmetric and asymmetric covariance recursions;
- model-validity and stationarity checks;
- BHHH estimation;
- score, Hessian, OPG, and QML covariance calculations;
- simulation and filtering;
- forecast paths and uncertainty bands;
- volatility impulse responses;
- VaR and coverage backtesting;
- multivariate portmanteau testing;
- Monte Carlo parameter-recovery evaluation;
- matrix construction and generalized-inverse utilities.

## Numerical score and Hessian

The original native code contains long generated analytical derivative routines. The Fortran port computes per-observation scores with adaptive central finite differences, falling back to one-sided differences when a perturbation violates the parameter domain. The likelihood Hessian is computed by second finite differences.

This preserves the public score, BHHH, covariance, t-ratio, and VIRF-inference workflows, but exact gradients, convergence paths, iteration counts, and final estimates can differ slightly from the R/C++ implementation.

## Linear algebra

RcppArmadillo operations are replaced by explicit BLAS/LAPACK calls and modern Fortran matrix operations. General inverses use an SVD-based Moore-Penrose pseudoinverse, allowing singular OPG and diagnostic matrices to be handled consistently.

## Starting values

The deterministic start reproduces the source package's `crossprod(data)/T` unconditional second moment and model-specific coefficient defaults. Randomized searches use a portable explicit Fortran RNG, so their streams are not identical to R's RNG.

## Forecast recursion

The first forecast uses the last realized return, as in the source. At horizons beyond one, the Fortran implementation propagates the conditional expectation `E[r r' | F] = H`. This is the standard analytical BEKK covariance forecast.

Some symmetric R methods instead construct a matrix square root of the prior forecast and use it as a pseudo-return. Consequently, multi-step symmetric forecasts may not be bit-for-bit identical, while the Fortran recursion has a direct conditional-expectation interpretation.

Forecast confidence paths follow the source package's endpoint-parameter approach: lower and upper parameter vectors are formed from fitted standard errors and then filtered/forecast separately. If an endpoint is invalid or nonstationary, the central path is used for that side rather than aborting the entire forecast.

## Volatility impulse responses

The Fortran VIRF propagates the covariance difference caused by the shock. It deliberately excludes a constant `C C'` term from the response recursion because a volatility impulse response is the difference between shocked and baseline covariance paths; the constant cancels.

Confidence bands use a numerical Jacobian of the complete VIRF path and the robust parameter covariance matrix when available. This follows the delta-method intent of the R implementation without relying on `numDeriv`.

Unlike the public R wrapper, the typed Fortran routine can also evaluate asymmetric models. The asymmetric indicator is applied to the initial shock, and its expected value is used in later propagation.

## Portmanteau test correction

The R source calculates the p-value by passing an already computed probability back into its `p_val_q` helper. The Fortran port instead computes the intended upper-tail chi-square probability directly from the test statistic and degrees of freedom.

## VaR and Student-t approximation

Normal and empirical residual quantiles follow the source design. The Student-t option infers degrees of freedom from residual excess kurtosis and standardizes the t quantile to unit variance. R interpolation metadata and arbitrary sample-size warnings are not part of the numerical API.

## R-only infrastructure omitted

- S3 dispatch and print/summary methods;
- plotting and `ggplot2` objects;
- `xts`, `zoo`, `ts`, and calendar metadata;
- progress bars and `future` parallel orchestration;
- automatic conversion of bundled `.RData` datasets;
- CRAN documentation-generation infrastructure.

The numerical arrays needed for plotting are available in result types. Original datasets and documentation remain under `original/`.

## Array orientation

R stores vectors and matrices in column-major order, as does Fortran, but the public time-series layout is made explicit:

- `data(T,N)`;
- `h(N,N,T)`;
- `residuals(T,N)`.

`C` is stored as a lower-triangular matrix and the intercept is `C C'`.
