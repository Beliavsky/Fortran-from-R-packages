# fBasics modern Fortran

This project is a modern Fortran numerical translation of the computational core of the R package **fBasics 4052.98**. Version 0.3.0 extends the translation with prewhitened HAC/GMM covariance estimation and the Andrews ARMA(1,1) plug-in bandwidth path, in addition to the probability, fitting, moment-estimation, smoothing, spatial-interpolation, and hypothesis-testing algorithms introduced in version 0.2.0.

Plotting, interactive GUI code, R classes, formulas, `timeSeries` metadata, and other R-specific infrastructure are intentionally excluded. `API_MAP.md` identifies direct implementations, numerical analogues, and remaining differences.

## Licensing

The source package declares `GPL (>= 2)`. This translation therefore uses **GPL-2.0-or-later**. The complete GNU GPL version 2 text is in `LICENSE`. Every Fortran source, application, example, and test file carries the GPL-2.0-or-later SPDX identifier and attribution to the fBasics package.

## Implemented and tested computational surface

### Matrix and time-series utilities

- LAPACK matrix inversion and SVD numerical rank
- Matrix one, spectral, and infinity norms
- Trace, `vec`, `vech`, Kronecker products, and triangular extraction
- Hilbert and Pascal matrices
- Positive-definiteness checks and eigenvalue-based repair
- Rectangular grid vectors
- Lag/lead matrices and polynomial distributed-lag regressors

### Descriptive and robust statistics

- Basic summaries and confidence intervals
- Mean, variance, standard deviation, quantiles, median, skewness, and excess kurtosis
- Row means, variances, standard deviations, skewness, kurtosis, extrema, products, and quantiles
- IQR, Bowley skewness, percentile kurtosis, and the first four sample L-moments
- Covariance and correlation matrices

### Random generation and utility functions

- Reproducible Park-Miller LCG
- Normal, Student-t, Gamma, chi-square, inverse-Gaussian, and generalized-inverse-Gaussian-related generation
- Heaviside, sign, approximate delta, boxcar, and ramp functions

### Stable distributions

- Nolan S1 characteristic function
- Density by Fourier inversion
- CDF by Gil-Pelaez inversion
- Quantiles by safeguarded inversion
- Chambers-Mallows-Stuck simulation
- Exact Normal and symmetric Cauchy special cases
- Empirical-characteristic-function fitting
- Bounded numerical maximum-likelihood fitting
- Numerical Hessian and covariance output

The stable implementation uses the S1 parameterization. It is a tested numerical implementation, not a claim of exact `stabledist` integration, optimization, or parameterization defaults.

### Generalized-lambda distributions

- Ramberg-Schmeiser density, CDF, quantile, simulation, mode, and quantile matching
- FMKL density, CDF, quantile, and simulation
- Five-parameter FM5 density, CDF, quantile, and simulation
- Robust quantile, maximum-likelihood, maximum-product-spacing, goodness-of-fit, and histogram fitting paths

### Generalized-hyperbolic family

- NIG, GH, hyperbolic, GHT, standardized GH, standardized NIG, and standardized GHT density/CDF/quantile/simulation paths
- Analytical and numerical moments
- Robust quantile-moment summaries
- Bounded fit wrappers for GH, HYP, GHT, SGH, SNIG, and SGHT
- Numerical Hessians and covariance matrices where the fit path supplies them

### Smoothing-density estimation

- Penalized cubic B-spline log-density estimation
- Automatic normalized density and CDF grids
- Density, CDF, quantile, and random-generation procedures

This is a self-contained numerical analogue of the package's smoothing-spline-density surface; it does not reproduce the external GSS optimizer iteration for iteration.

### GMM, GEL, and HAC calculations

- Identity, two-step, iterated, and continuously updated GMM
- EL, ET, CUE, and ETEL generalized empirical likelihood
- User-supplied moment callbacks over plain arrays
- IID and HAC moment covariance matrices
- Bartlett, Parzen, truncated, Tukey-Hanning, and Quadratic Spectral kernels
- Newey-West, Andrews AR(1), and Andrews ARMA(1,1) plug-in bandwidths
- VAR(p) prewhitening of moment conditions with exact long-run recoloring
- Automatic fallback to unprewhitened/Newey-West calculations for singular or degenerate fits
- Reported bandwidth, prewhitening order, and prewhitening status in `gmm_result`
- Numerical moment Jacobians, sandwich covariance matrices, J tests, and linear Wald restrictions
- Observation-level GEL weights

Example prewhitened GMM call:

```fortran
call fit_gmm(moment_fun, data, start, lower, upper, "two_step", result, &
  cov_type="HAC", kernel="Quadratic Spectral", bandwidth=-1, &
  prewhite_order=1, bandwidth_method="andrews_arma11")
```

The lower-level `moment_covariance` procedure accepts the same prewhitening and bandwidth controls. `prewhiten_var`, `fit_arma11_css`, and `andrews_bandwidth_value` are public for independent use and testing.

### Interpolation and spatial prediction

- One-dimensional linear interpolation
- Bilinear interpolation on rectangular grids
- Inverse-distance and local-plane scattered interpolation
- Exhaustive two-dimensional triangulated piecewise-linear interpolation
- Ordinary kriging with exponential, Gaussian, or spherical covariance models
- Estimated covariance range, sill, and nugget
- Kriging predictions, prediction variances, and unbiased weights

The triangulated interpolation is piecewise linear rather than the dependency's exact Akima cubic surface. The kriging implementation is ordinary kriging over plain coordinate arrays rather than an R formula/GLS wrapper.

### Hypothesis tests

- Pearson, Spearman, and Kendall correlation tests
- Pooled and Welch two-sample t-tests
- Two-sample variance F-test
- Jarque-Bera and Urzua finite-sample-adjusted Jarque-Bera tests
- D'Agostino skewness, kurtosis, and omnibus tests
- Shapiro-Wilk and Shapiro-Francia tests
- Ansari-Bradley, Mood, Bartlett, and Fligner-Killeen scale tests
- Wilcoxon rank-sum and two-sample Kruskal-Wallis tests
- Two-sample Kolmogorov-Smirnov test
- One-sample Normal KS and Pearson chi-square tests
- Cramer-von Mises, Anderson-Darling, and Lilliefors Normality tests

Several p-values use published asymptotic or response-surface approximations. The large historical LM/ALM lookup-table objects shipped with the R package are not duplicated.

### Maximum drawdown

- Brownian expected maximum drawdown using the package's Appendix-B interpolation values and asymptotic tails
- Trendless Brownian maximum-drawdown density and CDF series
- Brownian-path maximum-drawdown simulation
- Realized maximum drawdown from simple returns

## Build and validation

GNU Fortran, LAPACK, and BLAS are required.

```sh
make check
```

This runs the license audit, a runtime-checked debug build, five numerical suites and all applications, then repeats the build and tests with optimization and warnings treated as errors.

Individual commands:

```sh
make debug
make release
./scripts/run_tests.sh debug
./scripts/run_tests.sh release
```

## Applications

Run the demonstration:

```sh
build/release/demo_fbasics
```

Fit a one-column numeric CSV or a `Date,Value` CSV:

```sh
build/release/analyze_csv data/example_returns.csv normal
build/release/analyze_csv data/example_returns.csv student
build/release/analyze_csv data/example_returns.csv nig
build/release/analyze_csv data/example_returns.csv gld
build/release/analyze_csv data/example_returns.csv fmkl
build/release/analyze_csv data/example_returns.csv fm5
build/release/analyze_csv data/example_returns.csv gh
build/release/analyze_csv data/example_returns.csv hyp
build/release/analyze_csv data/example_returns.csv sgh
build/release/analyze_csv data/example_returns.csv snig
build/release/analyze_csv data/example_returns.csv stable
build/release/analyze_csv data/example_returns.csv stable-mle
build/release/analyze_csv data/example_returns.csv ssd
```

The CSV application reads the final comma-separated field from each line and skips nonnumeric headers.

## Important numerical differences

- Bounded Nelder-Mead replaces the R optimization stack.
- Stable density and CDF values use numerical Fourier integration in Nolan's S1 parameterization.
- NIG/GH-family CDFs use adaptive numerical integration.
- General GIG simulation uses a log-scale slice sampler; NIG simulation uses its inverse-Gaussian mixture directly.
- Long-memory or formula-based model infrastructure is not relevant to this package and is not introduced.
- The spline-density estimator uses penalized B-splines rather than exact GSS internals.
- Ordinary kriging and triangulated linear interpolation replace dependency-backed R spatial interfaces.
- Andrews ARMA(1,1) parameters are estimated by bounded conditional-sum-of-squares rather than R's `arima` maximum likelihood.
- VAR prewhitening uses no-intercept least squares on centered moment conditions and recolors with `(I - A_1 - ... - A_p)^{-1}`.
- Random streams are reproducible within this Fortran implementation but do not reproduce R streams.

## Remaining exclusions

The remaining omissions are predominantly nonnumerical infrastructure or exact legacy-data reproduction:

- Plotting, sliders, GUI code, palettes, symbols, and interactive locator utilities
- S3/S4 containers such as `fDISTFIT` and `fHTEST`
- R formulas, model frames, printing methods, `timeSeries`/`timeDate` metadata, and package datasets
- Exact Akima cubic interpolation internals
- Exact `stabledist`, GSS, `gmm`, `sandwich`, and R optimizer/random-stream equivalence
- GMM formula parsing and exact R `ar`/`arima` optimizer behavior
- Large historical finite-sample Jarque-Bera LM/ALM lookup tables
- Character-table, color, and package-introspection helpers that have no meaningful plain-array numerical Fortran counterpart

`fpm.toml` is supplied, but `fpm` was unavailable in the validation environment and is not claimed as tested.
