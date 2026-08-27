# Fortran from R packages

This repository collects experimental modern Fortran translations and
ports by ChatGPT on High mode of computational code in 445 package
directories derived from R packages. Each subdirectory is an
independent Fortran Package Manager (fpm) project with its own
documentation, tests, provenance record, and license.

These projects are unofficial and are not endorsed by the original package
authors, CRAN, or the R Foundation. They have not been validated by a human
unless the documentation says otherwise.

Browse packages through the [domain-oriented task views](TASK_VIEWS/), the
[category and keyword index](PROJECT_INDEX.md), or the
[NIST GAMS mathematical-software index](GAMS_INDEX.md).

The [R versus Fortran comparison suite](comparisons/) runs deterministic
correctness checks and exploratory kernel timings against installed R packages.
It currently covers probability distributions, fractional time-series methods,
fixed-parameter GARCH filters, financial performance measures, time-series
diagnostics, covariance shrinkage, clustering, and the separate `mixtools`
mixture-model comparison. Fifty-eight cases use a shared dated asset-price fixture,
including fitted sGARCH-t, GJR-GARCH-t, and NAGARCH-t models.
These checks
are evidence for selected operations, not validation of every translated API.

## Downloading one package

Because the complete repository is large, Git users can download a single
package directory and its subdirectories with a shallow sparse clone. For
example, these Command Prompt commands download `tsa` without downloading the
file contents of the other packages:

```bat
git clone --depth 1 --filter=blob:none --sparse https://github.com/Beliavsky/Fortran-from-R-packages.git tsa-download
cd tsa-download
git sparse-checkout set tsa
```

Replace `tsa` in the final command with the desired package directory. Git's
cone-mode sparse checkout also retains a few files from the repository root,
such as this README.

Alternatively, download [`download_build_package.py`](download_build_package.py)
and give it a package directory name. The script performs the sparse clone,
includes any repository-local FPM dependencies, checks that the package exists,
and runs both `fpm build` and `fpm run`:

```bat
python download_build_package.py rugarch
```

It creates `rugarch-download` in the current directory and refuses to overwrite
an existing destination.

## Building on Windows

Most packages build normally with `fpm build`. Some numerical packages declare
system BLAS, LAPACK, or ARPACK libraries in `fpm.toml`. On Windows, run
`setup_windows_linalg.bat` once in the same Command Prompt before building such
packages. The script detects compatible import libraries supplied by GNU Octave
or JAGS, adds their `lib` directory to `FPM_LDFLAGS`, and adds their `bin`
directory to `PATH`. GNU Octave is preferred because it also supplies ARPACK.

The configuration affects only the current Command Prompt. It does not copy
binary libraries into this repository or commit machine-specific paths.

## Why Fortran?

Modern Fortran is well suited to the numerical kernels found in statistical R
packages. It provides native multidimensional arrays, column-major storage
compatible with R, concise mathematical code, and strong performance for
numerical computation.

C++ has stronger language-binding and distribution ecosystems, so the
translations can be exposed through portable `bind(C)` wrappers and thin
adapters for R, Python, MATLAB, and Octave without rewriting the working
Fortran implementations. See [Why Fortran rather than
C++?](docs/FORTRAN_VS_CPP.md) for the tradeoffs and project rationale.

## Shared numerical modules

The [`rfortran-core`](rfortran-core/) library is the first step toward removing
duplicate R-like numerical helpers from the translations. Its dependency-free
modules currently provide kinds, missing-value policies, weighted and unweighted
descriptive statistics, normal and central Student-t/chi-square/F distribution
functions, stable
ordering, type-7 quantiles and medians, median absolute deviation,
weighted-ECDF quantiles, several explicitly named weighted-quantile
interpolation conventions,
stable log reductions and elementary transforms, positive-domain and
integer-combinatorial special functions, regularized gamma/beta functions,
paired covariance/correlation,
trailing rolling sums, extrema and moments, differencing, and time-series covariance functions. The core now has a
direct deterministic comparison suite against R reference implementations;
`FinTS`, `fracdiff`, and
`rugarch` retain package-level R comparisons, while additional migrated
packages exercise compatibility wrappers. See the
[shared-module design and migration
notes](docs/SHARED_MODULES.md).

The optional [`rfortran-linalg`](rfortran-linalg/) library centralizes checked
linear solves, general matrix inverses, symmetric eigendecompositions,
Cholesky factors, SPD inverse/log-determinant calculations, real and complex
thin SVD, economy-size QR, singular values, numerical rank, rank-revealing
pivoted QR, QR least squares, SVD least squares, full SVD, and spectral radius.
Pivoted-QR least squares is also available when callers need a numerical rank
and a selectable rank threshold.
It also provides general real and complex eigendecomposition and real Schur
decomposition, complex Schur decomposition, complex linear solves, and real
and complex matrix balancing.
Triangular inversion, general determinants, and signed log-absolute-
determinants are also provided.
It uses a pinned pure-Fortran
LAPACK backend, so migrated packages do not require system `-llapack` or
`-lblas` libraries. Basic operations already supplied by the language, such
as `norm2`, remain intrinsic calls rather than shared wrappers. Current
migrations cover `CEoptim`, `cmaes`, `cccp`, `bayesianOU`, `CLA`,
`RiskPortfolios`, `frapo`, `riskParityPortfolio`, `stochfactor`, `ks`, `nmof`,
`mixsqp`, `garchx`, `tvgarch`, `fbonds`, `fnonlinear`, `rquantlib`,
`MultiATSM`, `BEKKs`, `compositions`, `fmultivar`, `fcopulae`, `fbasics`,
`irlba`, `msm`, `matrixdist`, `etrm`, `esback`, `apt`, `matchingMarkets`, `gmm`,
`nnet`, `Rmalschains`, `statmod`, `robustbase`, `fastmatrix`, `L1pack`,
`lgarch`, `tsdyn`, `expm`, `pa`, `fportfolio`, `Rcsdp`, `Rdsdp`, `gogarch`,
and `mclust`.

`nleqslv` requires low-level solver-specific BLAS/LAPACK kernels rather than
the checked high-level API. It uses the same pinned pure-Fortran LAPACK backend
directly and therefore also avoids system `-llapack` and `-lblas` libraries.

The internal [`rfortran-arpack`](rfortran-arpack/) package supplies the
double-precision ARPACK-NG 3.9.1 iterative eigensolvers used by `RSpectra` and
`bigstatsr`. Its upstream sources have been converted to free source form and
use the same pinned pure-Fortran LAPACK backend, avoiding system `-larpack`,
`-llapack`, and `-lblas` libraries in those packages.

## Language interfaces

A [pilot interface layer](interfaces/) exposes selected procedures from
`pbivnorm`, `normalp`, and `bondAnalyst` through explicit-shape `bind(C)`
wrappers, thin R `.C()` functions, a small Python/NumPy package, and a shared
Octave/MATLAB MEX gateway. The gateway is tested with Octave; MATLAB build and
convenience functions are supplied but remain unverified without a MATLAB
installation. Cross-language tests cover numerical results and shape-error
handling while leaving the existing computational implementations intact.

## Packages and licenses

| Package | What it does | Upstream version · License |
| --- | --- | --- |
| [`ABCoptim`](ABCoptim/) | Performs artificial bee colony optimization. | `0.15.0` · `MIT` |
| [`ACDm`](ACDm/) | Estimates and simulates autoregressive conditional duration models. | `1.1.0` · `GPL-3.0-or-later` |
| [`actuar`](actuar/) | Provides actuarial distributions, aggregate-loss and ruin models, coverage transformations, credibility methods, and minimum-distance estimation. | `3.3-7` · `GPL-2.0-or-later` |
| [`adagio`](adagio/) | Provides discrete and global optimization routines. | `0.9.2` · `GPL-3.0-or-later` |
| [`AdequacyModel`](AdequacyModel/) | Assesses probabilistic-model adequacy and provides general-purpose optimization routines. | `2.0.0` · `GPL-2.0-or-later` |
| [`ADGofTest`](ADGofTest/) | Performs the Anderson-Darling goodness-of-fit test. | `0.3` · `GPL` |
| [`AEP`](AEP/) | Models asymmetric exponential-power distributions and performs robust regression. | `0.1.4` · `GPL-2.0-or-later` |
| [`alabama`](alabama/) | Performs constrained nonlinear optimization. | `2025.1.0` · `GPL-2.0-or-later` |
| [`anMC`](anMC/) | Computes high-dimensional orthant probabilities. | `0.2.5` · `GPL-3.0-only` |
| [`ao`](ao/) | Performs alternating optimization. | `1.2.3` · `GPL-3.0-only` |
| [`apt`](apt/) | Models asymmetric price transmission with threshold cointegration and error-correction methods. | `4.0` · `GPL-2.0-or-later` |
| [`arfima`](arfima/) | Fits, simulates, filters, and forecasts long-memory ARFIMA models. | `1.8-2` · `MIT` |
| [`argus`](argus/) | Provides probability functions and random generation for the Argus distribution. | `0.1.1` · `GPL-2.0-or-later` |
| [`backtest`](backtest/) | Explores portfolio-based conjectures about financial instruments. | `0.3-4` · `GPL-2.0-or-later` |
| [`bayesgarch`](bayesgarch/) | Performs Bayesian estimation of GARCH models with Student-t innovations. | `2.1.10` · `GPL-2.0-or-later` |
| [`bayesianOU`](bayesianOU/) | Fits Bayesian nonlinear Ornstein-Uhlenbeck models. | `0.2.0` · `MIT` |
| [`bayesm`](bayesm/) | Provides Bayesian inference for marketing and micro-econometric models. | `3.1-7` · `GPL-2.0-or-later` |
| [`BB`](BB/) | Solves and optimizes large-scale nonlinear systems. | `2026.1.0` · `GPL-3.0-only` |
| [`bcc1997`](bcc1997/) | Prices European options with stochastic volatility, rates, and jumps. | `0.1.1` · `GPL-2.0-or-later` |
| [`BEKKs`](BEKKs/) | Estimates and analyzes BEKK multivariate conditional-volatility models. | `1.4.7` · `MIT` |
| [`BenfordTests`](BenfordTests/) | Tests whether data conform to Benford's law. | `1.2.0` · `GPL-3.0-only` |
| [`betafunctions`](betafunctions/) | Works with two- and four-parameter beta distributions and psychometric classification analysis. | `1.9.0` · `CC0-1.0` |
| [`betategarch`](betategarch/) | Estimates Beta-t-EGARCH volatility models. | `3.4` · `GPL-2.0-only` |
| [`BGFD`](BGFD/) | Provides Bell-G and complementary Bell-G families of probability distributions. | `0.1` · `GPL-2.0-or-later` |
| [`BiasedUrn`](BiasedUrn/) | Provides Fisher and Wallenius noncentral hypergeometric distributions. | `2.0.12` · `GPL-3.0-only` |
| [`bidask`](bidask/) | Estimates bid-ask spreads efficiently from OHLC prices. | `2.1.5` · `MIT` |
| [`bigstatsr`](bigstatsr/) | Provides statistical tools for large file-backed matrices. | `1.6.2` · `GPL-3.0-only` |
| [`BivGeo`](BivGeo/) | Provides the Basu-Dhar bivariate geometric distribution. | `2.1.1` · `GPL-2.0-or-later` |
| [`bivgeom`](bivgeom/) | Provides Roy's bivariate geometric distribution. | `1.0` · `GPL` |
| [`bivpois`](bivpois/) | Provides distribution, estimation, and simulation methods for the bivariate Poisson distribution. | `1.2` · `GPL-2.0-or-later` |
| [`blmodel`](blmodel/) | Computes Black-Litterman posterior distributions. | `1.0.2` · `GPL-3.0-only` |
| [`bondAnalyst`](bondAnalyst/) | Performs fixed-income valuation and yield, spread, and duration calculations. | `1.0.1` · `GPL-3.0-only` |
| [`BondValuation`](BondValuation/) | Values fixed-coupon bonds with odd coupon periods and multiple day-count conventions. | `0.1.1` · `GPL-3.0-only` |
| [`boot`](boot/) | Provides bootstrap methods and functions. | `1.3-32` · `Unlimited` |
| [`bridgedist`](bridgedist/) | Implements the bridge distribution with a logit link. | `0.1.3` · `GPL-2.0-or-later` |
| [`bvls`](bvls/) | Solves bounded-variable least-squares problems with the Stark-Parker algorithm. | `1.4` · `GPL-2.0-or-later` |
| [`bzinb`](bzinb/) | Estimates bivariate zero-inflated negative-binomial models. | `1.0.8` · `GPL-2.0-only` |
| [`calibrar`](calibrar/) | Automates parameter estimation for complex models. | `0.9.0` · `GPL-2.0-only` |
| [`caRamel`](caRamel/) | Performs multi-objective evolutionary optimization. | `1.5` · `GPL-3.0-only` |
| [`cbinom`](cbinom/) | Provides the continuous analogue of the binomial distribution. | `1.6` · `GPL-2.0-or-later` |
| [`cccp`](cccp/) | Solves cone-constrained convex optimization problems. | `0.3-3` · `GPL-3.0-or-later` |
| [`CCd`](CCd/) | Provides the Cauchy-Cacoullos discrete Cauchy distribution. | `1.1` · `GPL-2.0-or-later` |
| [`cec2005benchmark`](cec2005benchmark/) | Provides the CEC 2005 real-parameter optimization benchmark suite. | `1.0.4` · `GPL-3.0-or-later` |
| [`cec2013`](cec2013/) | Provides CEC 2013 optimization benchmark functions. | `0.1-5` · `GPL-3.0-or-later` |
| [`CEoptim`](CEoptim/) | Performs optimization using the cross-entropy method. | `1.3` · `GPL-2.0-or-later` |
| [`CGNM`](CGNM/) | Fits nonlinear models using the cluster Gauss-Newton method. | `0.9.3` · `MIT` |
| [`ChernoffDist`](ChernoffDist/) | Evaluates Chernoff's distribution. | `0.1.0` · `GPL-3.0-only` |
| [`chyper`](chyper/) | Provides conditional hypergeometric distributions. | `0.3.1` · `MIT` |
| [`CircStats`](CircStats/) | Provides methods for circular statistics. | `0.2-7` · `GPL-2.0-only` |
| [`CLA`](CLA/) | Implements the Markowitz critical line algorithm for portfolio optimization. | `0.96-3` · `GPL-3.0-or-later` |
| [`clarabel`](clarabel/) | Provides a Fortran interface to the Clarabel conic interior-point solver. | `0.11.2` · `Apache-2.0` |
| [`clue`](clue/) | Combines and analyzes cluster ensembles. | `0.3-68` · `GPL-2.0-only` |
| [`cluster`](cluster/) | Provides clustering methods for finding groups in data. | `2.1.8.3` · `GPL-2.0-or-later` |
| [`cmaes`](cmaes/) | Performs optimization using covariance-matrix adaptation evolution strategies. | `1.0-12` · `GPL-2.0-only` |
| [`coda`](coda/) | Analyzes and diagnoses Markov-chain Monte Carlo output. | `0.19-4.1` · `GPL-2.0-or-later` |
| [`COMPoissonReg`](COMPoissonReg/) | Fits Conway-Maxwell-Poisson regression models. | `0.8.2` · `GPL-2.0-only OR GPL-3.0-only` |
| [`compositions`](compositions/) | Analyzes compositional and positive data. | `2.0-9` · `GPL-2.0-or-later` |
| [`compound.Cox`](compound.Cox/) | Performs survival feature screening and compound-covariate prediction, including copula-based dependent-censoring analyses. | `3.33` · `GPL-2.0-only` |
| [`coneproj`](coneproj/) | Performs cone projections, quadratic programming, and shape-restricted regression. | `1.23` · `GPL-2.0-or-later` |
| [`contfrac`](contfrac/) | Evaluates and manipulates continued fractions. | `1.1-12` · `GPL-2` |
| [`copula`](copula/) | Models multivariate dependence with common copula families. | `1.1-7` · `GPL-3.0-or-later` |
| [`corpcor`](corpcor/) | Estimates covariance and partial-correlation matrices efficiently. | `1.6.10` · `GPL-3.0-or-later` |
| [`corpmetrics`](corpmetrics/) | Provides corporate valuation, financial-metric, and modelling calculations. | `1.0` · `GPL-2.0-or-later` |
| [`countDM`](countDM/) | Estimates statistical models for count data. | `0.1.0` · `GPL-2.0-or-later` |
| [`creditr`](creditr/) | Values credit-default swaps using the ISDA CDS Standard Model. | `0.6.2` · `GPL-3.0-only AND LicenseRef-ISDA-CDS-Standard-Model` |
| [`cubature`](cubature/) | Performs adaptive multivariate integration over hypercubes. | `2.1.4-1` · `GPL-3.0-or-later` |
| [`cvar`](cvar/) | Computes value at risk and expected shortfall from distributions or samples. | `0.6` · `GPL-2.0-or-later` |
| [`degreenet`](degreenet/) | Models and analyzes network degree distributions. | `1.3-7` · `GPL-3.0-or-later` |
| [`Delaporte`](Delaporte/) | Provides distribution, probability, quantile, and random-generation functions for the Delaporte distribution. | `8.4.3` · `BSD-2-Clause` |
| [`DEoptim`](DEoptim/) | Performs global optimization by differential evolution. | `2.2-8` · `GPL-2.0-or-later` |
| [`deoptimr`](deoptimr/) | Performs global optimization using differential evolution. | `1.2-0` · `GPL-2.0-or-later` |
| [`derivmkts`](derivmkts/) | Provides derivative pricing and financial-market calculations. | `0.2.5.1` · `MIT` |
| [`desirability`](desirability/) | Optimizes and ranks solutions using desirability functions. | `2.1` · `GPL-2.0-only` |
| [`deSolve`](deSolve/) | Solves initial-value problems for differential equations. | `1.42` · `GPL-2.0-or-later` |
| [`dfoptim`](dfoptim/) | Provides derivative-free optimization algorithms. | `2023.1.0` · `GPL-2.0-or-later` |
| [`DiceDesign`](DiceDesign/) | Constructs designs for computer experiments. | `1.10` · `GPL-3.0-only` |
| [`DiceKriging`](DiceKriging/) | Fits Gaussian-process models for computer experiments. | `1.6.1` · `GPL-2.0-only OR GPL-3.0-only` |
| [`Directional`](Directional/) | Provides efficient statistical and mathematical functions, including directional-data methods. | `2.1.5.2` · `GPL-2.0-or-later` |
| [`DirichletReg`](DirichletReg/) | Fits Dirichlet regression models. | `0.7-2` · `GPL-2.0-or-later` |
| [`dirmult`](dirmult/) | Estimates Dirichlet-multinomial distributions. | `0.1.3-5` · `GPL-2.0-or-later` |
| [`DiscreteDists`](DiscreteDists/) | Provides discrete statistical distributions. | `1.1.2` · `MIT` |
| [`DiscreteInverseWeibull`](DiscreteInverseWeibull/) | Provides the discrete inverse Weibull distribution. | `1.0.2` · `GPL-2.0-only` |
| [`DiscreteLaplace`](DiscreteLaplace/) | Provides discrete Laplace distributions. | `1.1.1` · `GPL` |
| [`DiscreteWeibull`](DiscreteWeibull/) | Provides type 1 and type 3 discrete Weibull distributions. | `1.1` · `GPL-2.0-only` |
| [`distr`](distr/) | Provides composable probability-distribution objects and calculations. | `2.9.7` · `LGPL-3.0-only` |
| [`Dowd`](Dowd/) | Provides quantitative financial risk-management calculations. | `0.12` · `GPL-2.0-only OR GPL-3.0-only` |
| [`Dykstra`](Dykstra/) | Solves quadratic programs using cyclic projections. | `1.0-0` · `GPL-2.0-or-later` |
| [`ecd`](ecd/) | Models elliptic lambda distributions and prices options. | `0.9.2.4` · `Artistic-2.0` |
| [`ECOSolveR`](ECOSolveR/) | Solves conic optimization problems with an embedded conic solver. | `0.6.1` · `GPL-3.0-or-later` |
| [`ecpdist`](ecpdist/) | Provides the extended Chen-Poisson lifetime distribution. | `0.2.1` · `GPL-3.0-only` |
| [`elliptic`](elliptic/) | Computes Weierstrass and Jacobi elliptic functions. | `1.5-1` · `GPL-2.0-only` |
| [`epo`](epo/) | Performs enhanced portfolio optimization with correlation shrinkage. | `0.1.0.9000` · `MIT` |
| [`esback`](esback/) | Backtests expected-shortfall forecasts. | `0.3.1` · `GPL-3.0-only` |
| [`etrm`](etrm/) | Values and analyzes energy-trading and risk-management instruments. | `1.0.2` · `MIT` |
| [`evir`](evir/) | Performs extreme-value analysis for financial risk. | `1.7-4` · `GPL-2.0-or-later` |
| [`ewens`](ewens/) | Provides the Ewens sampling distribution. | `0.1.0` · `GPL-3.0-or-later` |
| [`expint`](expint/) | Computes exponential integrals and incomplete gamma functions. | `0.2-1` · `GPL-3.0-or-later` |
| [`expm`](expm/) | Computes matrix exponentials, logarithms, square roots, and related functions. | `1.0-0` · `GPL-3.0-or-later` |
| [`extraDistr`](extraDistr/) | Provides additional univariate and multivariate probability distributions. | `1.10.0.5` · `GPL-2.0-only` |
| [`fastcluster`](fastcluster/) | Performs fast hierarchical agglomerative clustering. | `1.3.0` · `BSD-2-Clause` |
| [`fastmatrix`](fastmatrix/) | Computes matrices used in statistical methods efficiently. | `0.6-6` · `GPL-3.0-only` |
| [`fattailsr`](fattailsr/) | Provides Kiener distributions and fat-tail analytics. | `2.0.1` · `GPL-2.0-only` |
| [`fbasics`](fbasics/) | Provides financial-market statistics, distributions, and utilities. | `4052.98` · `GPL-2.0-or-later` |
| [`fbonds`](fbonds/) | Prices bonds and fits Nelson-Siegel family term structures. | `3042.78` · `GPL-2.0-or-later` |
| [`fcl`](fcl/) | Provides dated cash-flow, bond, yield, duration, and financial-calendar calculations. | `0.1.4` · `MIT` |
| [`fcopulae`](fcopulae/) | Models dependence with elliptical, Archimedean, and empirical copulas. | `4052.86` · `GPL-2.0-or-later` |
| [`FER`](FER/) | Provides financial-engineering option-pricing formulas. | `0.94` · `GPL-2.0-or-later` |
| [`fextremes`](fextremes/) | Models extreme values and financial tail risk. | `4032.84` · `GPL-2.0-or-later` |
| [`ffp`](ffp/) | Computes fully flexible probabilities for stress testing and portfolio construction. | `0.2.2.9000` · `MIT` |
| [`fGarch`](fGarch/) | Fits and analyzes univariate GARCH and APARCH models. | `4052.93` · `GPL-2.0-or-later` |
| [`fhmm`](fhmm/) | Fits hidden Markov models to financial data. | `1.4.3` · `GPL-3.0-only` |
| [`financialmath`](financialmath/) | Provides financial mathematics for actuaries. | `0.1.1` · `GPL-2.0-only` |
| [`fincal`](fincal/) | Provides time-value-of-money and computational-finance calculations. | `0.6.3` · `GPL-2.0-or-later` |
| [`FinCovRegularization`](FinCovRegularization/) | Estimates and regularizes covariance matrices for finance. | `1.1.0` · `GPL-2.0-only` |
| [`fingraph`](fingraph/) | Learns graph structures for financial markets. | `0.1.0` · `GPL-3.0-only` |
| [`FinTS`](FinTS/) | Provides methods accompanying analysis of financial time series. | `0.4-9` · `GPL-2.0-or-later` |
| [`fitdistrplus`](fitdistrplus/) | Fits parametric distributions to censored and uncensored data. | `1.2-6` · `GPL-2.0-or-later` |
| [`fitHeavyTail`](fitHeavyTail/) | Estimates means and covariance matrices under heavy tails. | `0.2.0.9000` · `GPL-3.0-only` |
| [`FKF`](FKF/) | Performs fast multivariate Kalman filtering and smoothing. | `0.2.6` · `GPL-2.0-or-later` |
| [`flexsurv`](flexsurv/) | Fits flexible parametric survival and multi-state models. | `2.3.2` · `GPL-2.0-or-later` |
| [`FLSSS`](FLSSS/) | Solves subset-sum, multidimensional-knapsack, and generalized-assignment problems. | `9.2.8` · `GPL-3.0-only` |
| [`fmbasics`](fmbasics/) | Provides foundational financial-market calculations. | `0.3.99` · `GPL-2.0-only` |
| [`fmultivar`](fmultivar/) | Provides multivariate distributions and financial-data analysis. | `4031.84` · `GPL-2.0-or-later` |
| [`FNN`](FNN/) | Provides fast nearest-neighbor search algorithms. | `1.1.4.1` · `GPL-2.0-or-later` |
| [`fnonlinear`](fnonlinear/) | Models nonlinear and chaotic time series. | `4052.83` · `GPL-2.0-or-later` |
| [`fportfolio`](fportfolio/) | Performs portfolio selection, optimization, risk analysis, and backtesting. | `4023.84` · `GPL-2.0-or-later` |
| [`fracdiff`](fracdiff/) | Estimates, simulates, and analyzes fractionally differenced time-series models. | `1.5-4` · `GPL-2.0-or-later` |
| [`frapo`](frapo/) | Provides financial risk modelling and portfolio optimization methods. | `0.4-2` · `GPL-3.0-or-later` |
| [`frbinom`](frbinom/) | Provides fractional binomial distributions. | `1.0.0` · `MIT` |
| [`GA`](GA/) | Performs optimization with genetic algorithms. | `3.2.5` · `GPL-2.0-or-later` |
| [`gamlss`](gamlss/) | Fits generalized additive models for location, scale, and shape. | `5.5-0` · `GPL-3.0-only` |
| [`gamlss.dist`](gamlss.dist/) | Provides distributions for generalized additive location-scale-shape models. | `6.1-1` · `GPL-3.0-only` |
| [`GARCHIto`](GARCHIto/) | Estimates unified and realized GARCH-Ito volatility models. | `0.1.0` · `GPL-3.0-only` |
| [`garchsk`](garchsk/) | Estimates GARCH models with conditional skewness and kurtosis. | `0.1.0` · `GPL-2.0-or-later` |
| [`garchx`](garchx/) | Fits GARCH models with exogenous covariates. | `1.7` · `GPL-2.0-or-later` |
| [`GB2`](GB2/) | Provides properties, likelihoods, and estimation for generalized beta distributions of the second kind. | `2.1.2` · `GPL-2.0-or-later` |
| [`gcpm`](gcpm/) | Models credit-portfolio risk analytically and by Monte Carlo simulation. | `1.2.2` · `GPL-2.0-only` |
| [`genalg`](genalg/) | Performs optimization with genetic algorithms. | `0.2.1` · `GPL-2.0-only` |
| [`GenBinomApps`](GenBinomApps/) | Computes generalized-binomial probabilities and Clopper-Pearson confidence intervals. | `1.2.1` · `GPL-3.0-only` |
| [`GeneralizedHyperbolic`](GeneralizedHyperbolic/) | Provides generalized hyperbolic and related probability distributions. | `0.8-7` · `GPL-2.0-or-later` |
| [`GenSA`](GenSA/) | Performs global optimization using generalized simulated annealing. | `1.1.15` · `GPL-2.0-only` |
| [`ghyp`](ghyp/) | Evaluates, fits, and simulates generalized hyperbolic distributions. | `1.6.5` · `GPL-2.0-or-later` |
| [`gkwdist`](gkwdist/) | Provides the generalized Kumaraswamy distribution family. | `1.1.4` · `MIT` |
| [`glmnet`](glmnet/) | Fits lasso and elastic-net regularized generalized linear models. | `5.0` · `GPL-2.0-only` |
| [`globalOptTests`](globalOptTests/) | Provides objective functions for benchmarking global optimization methods. | `1.1` · `GPL-3.0-or-later` |
| [`gnorm`](gnorm/) | Evaluates and simulates generalized normal distributions. | `1.0.2` · `GPL-2.0-or-later` |
| [`goftest`](goftest/) | Performs classical goodness-of-fit tests for univariate distributions. | `1.2-3` · `GPL-2.0-or-later` |
| [`gogarch`](gogarch/) | Fits generalized orthogonal GARCH models. | `0.7-6` · `GPL-2.0-or-later` |
| [`good`](good/) | Fits Good regression models. | `1.0.2` · `GPL-2.0-or-later` |
| [`GPareto`](GPareto/) | Uses Gaussian processes for Pareto-front estimation and optimization. | `1.1.9` · `GPL-3.0-only` |
| [`GPArotation`](GPArotation/) | Provides gradient-projection methods for orthogonal and oblique factor rotation. | `2026.8-2` · `GPL-2.0-or-later` |
| [`graDiEnt`](graDiEnt/) | Performs stochastic quasi-gradient differential-evolution optimization. | `1.0.1` · `MIT` |
| [`greeks`](greeks/) | Computes option sensitivities, implied volatility, and Monte Carlo Greeks. | `1.5.6` · `MIT` |
| [`greeks1`](greeks1/) | Computes option sensitivities and implied volatilities. | `1.5.6` · `MIT` |
| [`greybox`](greybox/) | Provides regression model-building and forecasting tools. | `2.0.8` · `LGPL-2.1-only` |
| [`gsl`](gsl/) | Wraps GNU GSL special functions and numerical routines. | `2.1-9` · `GPL-3.0-only` |
| [`gslnls`](gslnls/) | Performs multistart nonlinear least-squares fitting. | `1.4.2` · `LGPL-3.0-only` |
| [`hdshop`](hdshop/) | Constructs high-dimensional shrinkage optimal portfolios. | `0.1.7` · `GPL-3.0-only` |
| [`hermite`](hermite/) | Provides the generalized Hermite distribution. | `1.2.1` · `GPL-2.0-or-later` |
| [`HierPortfolios`](HierPortfolios/) | Constructs portfolios using hierarchical risk clustering. | `1.0.2` · `GPL-2.0-only` |
| [`highfrequency`](highfrequency/) | Analyzes high-frequency trade and quote data. | `1.0.2` · `GPL-2.0-or-later` |
| [`highOrderPortfolios`](highOrderPortfolios/) | Designs portfolios using mean, variance, skewness, and kurtosis. | `0.1.1` · `GPL-3.0-only` |
| [`highs`](highs/) | Provides a Fortran interface to the HiGHS optimization solver. | `1.14.0-2` · `GPL-2.0-or-later` |
| [`hyper2`](hyper2/) | Provides hyperdirichlet distributions and likelihood calculations. | `3.2-3` · `GPL-3.0-or-later` |
| [`hypergeo`](hypergeo/) | Computes the Gauss hypergeometric function. | `1.2-14` · `GPL-2.0-only` |
| [`ICSNP`](ICSNP/) | Provides tools for multivariate nonparametric statistics. | `1.1-3` · `GPL-2.0-or-later` |
| [`igraph`](igraph/) | Provides graph and network analysis algorithms. | `2.3.3` · `GPL-2.0-or-later` |
| [`imputeFin`](imputeFin/) | Imputes missing values and outliers in financial time series. | `0.1.2.9000` · `GPL-3.0-only` |
| [`IndGenErrors`](IndGenErrors/) | Tests independence between innovations of generalized-error models. | `0.1.6` · `GPL-3.0-only` |
| [`INFOSET`](INFOSET/) | Computes informative distribution sets for asset returns. | `4.1.1` · `GPL-2.0-or-later` |
| [`intradayModel`](intradayModel/) | Models and forecasts periodic financial intraday signals. | `0.0.1` · `Apache-2.0` |
| [`intrinsicFRP`](intrinsicFRP/) | Estimates and tests factor-model asset-pricing relationships. | `2.1.0` · `GPL-3.0-or-later` |
| [`invgamstochvol`](invgamstochvol/) | Computes likelihoods and posterior smoothing for inverse-gamma stochastic-volatility models. | `1.0.0` · `MIT` |
| [`irlba`](irlba/) | Computes fast truncated singular-value decompositions. | `2.3.7` · `GPL-3.0-or-later` |
| [`isotone`](isotone/) | Performs isotone optimization and inequality-restricted regression. | `1.1-2` · `GPL-2.0-only` |
| [`Jdmbs`](Jdmbs/) | Prices options by Monte Carlo under geometric Brownian and jump-diffusion models. | `1.4` · `GPL-2.0-or-later` |
| [`jfe`](jfe/) | Analyzes financial and economic time-series data. | `2.5.11` · `GPL-2.0-or-later` |
| [`joker`](joker/) | Provides probability distributions and parameter-estimation methods. | `0.14.2` · `GPL-3.0-or-later` |
| [`jrvFinance`](jrvFinance/) | Provides NPV, IRR, annuity, bond-pricing, and Black-Scholes calculations. | `1.4.3` · `GPL-2.0-or-later` |
| [`JumpTest`](JumpTest/) | Detects jumps in financial time series. | `1.1` · `MIT` |
| [`kdensity`](kdensity/) | Performs parametrically guided kernel-density estimation with asymmetric kernels. | `1.2.0` · `MIT` |
| [`kernlab`](kernlab/) | Provides kernel-based machine-learning algorithms. | `0.9-33` · `GPL-2.0-only` |
| [`KernSmooth`](KernSmooth/) | Performs kernel density estimation and local polynomial regression. | `2.23-27` · `Unlimited` |
| [`kofnGA`](kofnGA/) | Uses a genetic algorithm for fixed-size subset selection. | `1.3` · `GPL-2.0-only` |
| [`KrigInv`](KrigInv/) | Performs Gaussian-process-based inversion and contour estimation. | `1.4.2` · `GPL-3.0-only` |
| [`ks`](ks/) | Provides multivariate kernel-smoothing methods. | `1.15.3` · `GPL-2.0-only` |
| [`L1pack`](L1pack/) | Provides routines for L1 estimation. | `0.62-4` · `GPL-3.0-only` |
| [`LaplacesDemon`](LaplacesDemon/) | Provides Bayesian inference and Markov-chain Monte Carlo algorithms. | `16.1.8` · `MIT` |
| [`lbfgs`](lbfgs/) | Performs limited-memory BFGS and orthant-wise optimization. | `1.2.1.2` · `GPL-2.0-or-later` |
| [`ldhmm`](ldhmm/) | Fits hidden Markov models based on lambda distributions. | `0.6.1` · `Artistic-2.0` |
| [`leaps`](leaps/) | Performs regression subset selection. | `3.2` · `GPL-2.0-or-later` |
| [`lgarch`](lgarch/) | Simulates and estimates log-GARCH models. | `0.7` · `GPL-2.0-only` |
| [`lifeinsurer`](lifeinsurer/) | Models traditional life-insurance contracts. | `1.0.1` · `GPL-2.0-or-later` |
| [`limSolve`](limSolve/) | Solves linear inverse models with equality and inequality constraints. | `2.0.3` · `GPL` |
| [`LindleyPowerSeries`](LindleyPowerSeries/) | Provides Lindley power-series distributions. | `1.0.1` · `GPL-2.0-or-later` |
| [`linprog`](linprog/) | Solves linear programming problems. | `0.9-6` · `GPL-2.0-or-later` |
| [`lme4`](lme4/) | Fits linear, generalized linear, and nonlinear mixed-effects models. | `2.1-0` · `GPL-2.0-or-later` |
| [`lmomco`](lmomco/) | Fits and evaluates probability distributions using L-moments. | `2.5.7` · `GPL` |
| [`Lmoments`](Lmoments/) | Computes L-moments and fits probability distributions using them. | `1.3-2` · `GPL-2.0-only` |
| [`locfit`](locfit/) | Performs local regression, likelihood, and density estimation. | `1.5-9.12` · `GPL-2.0-or-later` |
| [`longmemo`](longmemo/) | Provides statistical methods for long-memory processes. | `1.1-4` · `GPL-2.0-or-later` |
| [`LowRankQP`](LowRankQP/) | Solves low-rank quadratic programming problems. | `1.0.6` · `GPL-2.0-or-later` |
| [`lpSolve`](lpSolve/) | Solves linear and integer programming problems. | `5.6.23.9000` · `LGPL-2.0-only` |
| [`lsei`](lsei/) | Solves least-squares and quadratic-programming problems under constraints. | `1.3-1` · `GPL-2.0-or-later` |
| [`lsmontecarlo`](lsmontecarlo/) | Prices American options using least-squares Monte Carlo. | `1.0` · `GPL-3.0-only` |
| [`ltsa`](ltsa/) | Provides methods for linear time-series analysis. | `1.4.6.1` · `GPL-2.0-or-later` |
| [`magic`](magic/) | Creates and analyzes magic squares, hypercubes, and Latin squares. | `1.6-1-1` · `GPL-2.0-only` |
| [`ManifoldOptim`](ManifoldOptim/) | Performs optimization on Riemannian manifolds. | `1.0.2` · `GPL-2.0-or-later` |
| [`markowitzr`](markowitzr/) | Performs statistical inference for Markowitz portfolios. | `1.0.2.0002` · `LGPL-3.0-or-later` |
| [`marqLevAlg`](marqLevAlg/) | Performs parallelized Marquardt-Levenberg optimization. | `2.0.8` · `GPL-2.0-or-later` |
| [`MASS`](MASS/) | Provides statistical methods from Venables and Ripley's MASS. | `7.3-66` · `GPL-3.0-only` |
| [`matchingMarkets`](matchingMarkets/) | Analyzes stable matching markets. | `1.0-5` · `GPL-2.0-or-later` |
| [`matchingR`](matchingR/) | Provides algorithms for stable and optimal matching. | `2.0.0` · `GPL-2.0-or-later` |
| [`matlab`](matlab/) | Provides numerical and matrix utilities modelled after MATLAB functions. | `1.0.4.1` · `Artistic-2.0` |
| [`Matrix`](Matrix/) | Provides dense and sparse matrix algorithms. | `1.7-6` · `GPL-3.0-only` |
| [`MatrixExtra`](MatrixExtra/) | Provides additional methods for sparse matrices. | `0.1.15` · `GPL-3.0-only` |
| [`maxLik`](maxLik/) | Performs maximum-likelihood estimation and related calculations. | `1.5-2.2` · `GPL-2.0-or-later` |
| [`mbbefd`](mbbefd/) | Provides Maxwell-Boltzmann, Bose-Einstein, and Fermi-Dirac distributions. | `0.8.14` · `GPL-2.0-only` |
| [`mc2d`](mc2d/) | Provides tools for two-dimensional Monte Carlo simulations. | `0.2.2` · `GPL-2.0-or-later` |
| [`mcga`](mcga/) | Performs real-valued optimization with machine-coded genetic algorithms. | `3.0.9` · `GPL-2.0-or-later` |
| [`mclust`](mclust/) | Performs Gaussian-mixture modelling and model-based clustering. | `6.1.3` · `GPL-2.0-or-later` |
| [`mcmc`](mcmc/) | Provides Markov-chain Monte Carlo algorithms. | `0.9-8` · `MIT` |
| [`MCMCpack`](MCMCpack/) | Provides Markov-chain Monte Carlo methods for statistical models. | `1.7-1` · `GPL-3.0-only` |
| [`mco`](mco/) | Provides multiple-criteria optimization algorithms. | `1.17` · `GPL-2.0-only` |
| [`mcrp`](mcrp/) | Constructs multiple-criteria risk-parity portfolios. | `0.0-1` · `GPL-3.0-only` |
| [`metaheuristicOpt`](metaheuristicOpt/) | Provides population-based metaheuristic optimization algorithms. | `2.0.0` · `GPL-2.0-or-later` |
| [`Metrics`](Metrics/) | Computes evaluation metrics for machine learning. | `0.1.4` · `BSD-3-Clause` |
| [`mev`](mev/) | Models extreme values. | `2.2` · `GPL-3.0-only` |
| [`mfGARCH`](mfGARCH/) | Fits mixed-frequency GARCH models. | `0.2.2` · `MIT` |
| [`mgcv`](mgcv/) | Fits generalized additive models with automatic smoothness selection. | `1.9-4` · `GPL-2.0-or-later` |
| [`minqa`](minqa/) | Provides derivative-free optimization by quadratic approximation. | `1.2.8` · `GPL-2.0-only` |
| [`miscTools`](miscTools/) | Provides miscellaneous numerical and statistical utilities. | `0.6-30` · `GPL-2.0-or-later` |
| [`MixedIndTests`](MixedIndTests/) | Tests randomness and independence for discrete, continuous, and mixed data. | `1.2.0` · `GPL-3.0-only` |
| [`mixSPE`](mixSPE/) | Fits mixtures of power-exponential and skew power-exponential distributions. | `0.9.3` · `GPL-2.0-only` |
| [`mixsqp`](mixsqp/) | Estimates mixture proportions using sequential quadratic programming. | `0.3-54` · `MIT` |
| [`mixtools`](mixtools/) | Analyzes finite mixture models. | `2.0.0` · `GPL-2.0-or-later` |
| [`mize`](mize/) | Provides unconstrained numerical optimization algorithms. | `0.2.5.9000` · `BSD-2-Clause` |
| [`mlr`](mlr/) | Provides computational infrastructure and native learners for machine learning. | `2.19.3` · `GPL-2.0-or-later` |
| [`mlrMBO`](mlrMBO/) | Performs model-based and Bayesian optimization. | `1.1.6` · `GPL-3.0-only` |
| [`MM`](MM/) | Provides the multiplicative multinomial distribution. | `1.7-0` · `GPL-2.0-only` |
| [`MNB`](MNB/) | Provides diagnostics for multivariate negative-binomial regression models. | `1.2.0` · `GPL-2.0-or-later` |
| [`mnormt`](mnormt/) | Provides multivariate normal and Student-t distributions, including truncated probabilities. | `2.1.2` · `GPL-2.0-only OR GPL-3.0-only` |
| [`moments`](moments/) | Computes moments, cumulants, skewness, kurtosis, and related tests. | `0.14.1` · `GPL-2.0-or-later` |
| [`msgarch`](msgarch/) | Fits Markov-switching GARCH models. | `2.51` · `GPL-2.0-or-later` |
| [`msm`](msm/) | Fits continuous-time multi-state and hidden Markov models. | `1.8.2` · `GPL-2.0-or-later` |
| [`mstate`](mstate/) | Prepares data and estimates and predicts multi-state models. | `0.3.3` · `GPL-2.0-or-later` |
| [`MTS`](MTS/) | Provides a toolkit for multivariate time-series analysis. | `1.2.1` · `Artistic-2.0` |
| [`multiAssetOptions`](multiAssetOptions/) | Values European and American multi-asset options by finite differences. | `0.1-2` · `GPL-2.0-only OR GPL-3.0-only` |
| [`MultiATSM`](MultiATSM/) | Models multicountry term structures of interest rates. | `1.5.1-2` · `GPL-2.0-only OR GPL-3.0-only` |
| [`MultiRNG`](MultiRNG/) | Generates multivariate pseudorandom samples. | `1.2.4` · `GPL-2.0-only OR GPL-3.0-only` |
| [`mvtnorm`](mvtnorm/) | Computes multivariate normal and Student-t probabilities, densities, and random samples. | `1.4-2` · `GPL-2.0-only` |
| [`n1qn1`](n1qn1/) | Performs unconstrained optimization using full-memory BFGS. | `6.0.1-14` · `CeCILL-2.1` |
| [`nbconv`](nbconv/) | Evaluates arbitrary convolutions of negative-binomial distributions. | `1.0.1` · `GPL-3.0-or-later` |
| [`neighbours`](neighbours/) | Provides neighbourhood functions for local-search algorithms. | `0.1-5` · `GPL-3.0-only` |
| [`neldermead`](neldermead/) | Performs derivative-free optimization using Nelder-Mead methods. | `1.0-13` · `CeCILL-2.0` |
| [`new.dist`](new.dist/) | Provides alternative continuous and discrete probability distributions. | `0.1.2` · `GPL-3.0-only` |
| [`NFCP`](NFCP/) | Estimates term structures with N-factor commodity-pricing models. | `1.2.2` · `GPL-3.0-only` |
| [`nilde`](nilde/) | Finds nonnegative integer solutions of linear Diophantine equations. | `1.1-7` · `GPL-2.0-or-later` |
| [`NlcOptim`](NlcOptim/) | Solves nonlinear optimization problems with nonlinear constraints. | `0.6` · `GPL-3.0-only` |
| [`nleqslv`](nleqslv/) | Solves systems of nonlinear equations. | `3.3.7` · `GPL-2.0-or-later` |
| [`nlme`](nlme/) | Fits linear and nonlinear mixed-effects models. | `3.1-170` · `GPL-2.0-or-later` |
| [`nloptr`](nloptr/) | Provides nonlinear optimization methods inspired by the NLopt interface. | `2.2.1.9000` · `LGPL-3.0-or-later` |
| [`nls2`](nls2/) | Fits nonlinear regression models using brute-force starting-value searches. | `0.3-4` · `GPL-2.0-only` |
| [`nlsic`](nlsic/) | Solves nonlinear least-squares problems with inequality constraints. | `1.2.0` · `GPL-2.0-only` |
| [`nlsr`](nlsr/) | Provides methods for solving nonlinear least-squares problems. | `2026.4.29` · `GPL-2.0-only` |
| [`nmof`](nmof/) | Provides numerical optimization methods for finance and economics. | `2.12-0` · `GPL-3.0-only` |
| [`nnls`](nnls/) | Solves nonnegative and mixed-sign least-squares problems. | `1.6` · `GPL-2.0-or-later` |
| [`nonneg-cg`](nonneg-cg/) | Performs nonnegative conjugate-gradient minimization. | `0.1.6-1` · `BSD-2-Clause` |
| [`normalp`](normalp/) | Provides exponential-power distributions, estimation, and regression. | `0.7.2.1` · `GPL` |
| [`nspmix`](nspmix/) | Fits nonparametric and semiparametric mixture models. | `2.0-0` · `GPL-2.0-or-later` |
| [`numDeriv`](numDeriv/) | Computes accurate numerical gradients, Jacobians, and Hessians. | `2016.8-1.1` · `GPL-2.0-or-later` |
| [`nvmix`](nvmix/) | Computes and simulates multivariate normal variance-mixture distributions. | `0.1-2` · `GPL-3.0-or-later` |
| [`ob-analytics`](ob-analytics/) | Analyzes limit order books and liquidity. | `0.1.2` · `GPL-2.0-or-later` |
| [`onls`](onls/) | Performs orthogonal nonlinear least-squares regression. | `0.1-4` · `GPL-2.0-or-later` |
| [`OOR`](OOR/) | Performs global optimization using optimistic optimization. | `0.1.4` · `LGPL` |
| [`opthedging`](opthedging/) | Values and optimally hedges call and put options. | `1.0` · `GPL-2.0-or-later` |
| [`optimflex`](optimflex/) | Performs derivative-based optimization with user-defined convergence criteria. | `0.1.8` · `MIT` |
| [`optimx`](optimx/) | Provides an expanded toolkit for general-purpose optimization. | `2025-4.9` · `GPL-2.0-only` |
| [`optionpricing`](optionpricing/) | Prices options with efficient simulation algorithms. | `0.1.2` · `GPL-2.0-only OR GPL-3.0-only` |
| [`optmatch`](optmatch/) | Performs optimal distance-based matching for observational studies. | `0.10.8` · `MIT` |
| [`orthopolynom`](orthopolynom/) | Provides orthogonal and orthonormal polynomial families. | `1.0-6.1` · `GPL-2.0-only` |
| [`osqp`](osqp/) | Solves convex quadratic programs using the OSQP algorithm. | `1.0.0` · `Apache-2.0` |
| [`pa`](pa/) | Performs equity-portfolio performance attribution. | `1.2-4` · `GPL-2.0-only` |
| [`parma`](parma/) | Provides portfolio allocation and risk-management applications. | `1.7` · `GPL-3.0-or-later` |
| [`partitions`](partitions/) | Generates and analyzes additive integer partitions. | `1.10-9` · `GPL` |
| [`pbivnorm`](pbivnorm/) | Evaluates the vectorized standard bivariate normal cumulative distribution function. | `0.6.0` · `GPL-2.0-or-later` |
| [`pbo`](pbo/) | Estimates the probability of backtest overfitting. | `1.3.5` · `MIT` |
| [`pdqutils`](pdqutils/) | Provides distribution approximations using Gram-Charlier, Edgeworth, and Cornish-Fisher expansions. | `0.1.6` · `LGPL-3.0-or-later` |
| [`peerperformance`](peerperformance/) | Performs luck-corrected peer-performance analysis. | `2.4.0` · `GPL-2.0-or-later` |
| [`performanceanalytics`](performanceanalytics/) | Provides econometric tools for performance and risk analysis. | `2.1.0` · `GPL-2.0-or-later` |
| [`pgnorm`](pgnorm/) | Evaluates and simulates p-generalized normal distributions. | `2.0.1` · `GPL-2.0-or-later` |
| [`PINstimation`](PINstimation/) | Estimates the probability of informed trading. | `0.2.0` · `GPL-3.0-or-later` |
| [`piqp`](piqp/) | Solves quadratic programs with a proximal interior-point method. | `0.6.2` · `GPL-3.0-only` |
| [`pmultinom`](pmultinom/) | Computes one-sided multinomial probabilities. | `1.0.0` · `AGPL-3.0-only` |
| [`PMwR`](PMwR/) | Provides portfolio accounting, return analysis, attribution, and backtesting tools. | `1.2-0` · `GPL-3.0-only` |
| [`poibin`](poibin/) | Provides the Poisson-binomial distribution. | `1.6` · `GPL-2.0-only` |
| [`poilog`](poilog/) | Provides Poisson-lognormal and bivariate Poisson-lognormal distributions. | `0.4.2.1` · `GPL-3.0-only` |
| [`PoissonBinomial`](PoissonBinomial/) | Computes ordinary and generalized Poisson-binomial distributions. | `1.2.8` · `GPL-3.0-only` |
| [`polyaAeppli`](polyaAeppli/) | Provides the Polya-Aeppli distribution. | `2.0.2` · `GPL-2.0-or-later` |
| [`polynom`](polynom/) | Provides arithmetic and calculations for univariate polynomials. | `1.4-1` · `GPL-2.0-only` |
| [`PortfolioAnalytics`](PortfolioAnalytics/) | Performs portfolio analysis and constrained optimization. | `2.1.2` · `GPL-3.0-only` |
| [`portfoliooptim`](portfoliooptim/) | Performs small- and large-sample portfolio optimization. | `1.1.1` · `GPL-3.0-only` |
| [`PortfolioTesteR`](PortfolioTesteR/) | Tests investment strategies using an English-like specification. | `0.1.4` · `MIT` |
| [`portvine`](portvine/) | Estimates portfolio risk with rolling ARMA-GARCH and vine-copula models. | `1.0.3.9000` · `GPL-3.0-only` |
| [`poweRlaw`](poweRlaw/) | Analyzes power-law and other heavy-tailed distributions. | `1.0.0` · `GPL-3.0-only` |
| [`ppcor`](ppcor/) | Computes partial and semi-partial correlations. | `1.1` · `GPL-2.0-only` |
| [`ppso`](ppso/) | Performs particle-swarm optimization and dynamically dimensioned search. | `0.9-99994` · `Unlimited` |
| [`pracma`](pracma/) | Provides practical numerical mathematics functions. | `2.4.6` · `GPL-3.0-or-later` |
| [`PSDistr`](PSDistr/) | Provides practical numerical mathematics functions. | `2.4.6` · `GPL-3.0-or-later` |
| [`pso`](pso/) | Performs particle-swarm optimization. | `1.0.4` · `LGPL-3.0-only` |
| [`psoptim`](psoptim/) | Performs particle-swarm optimization. | `1.0` · `GPL-2.0-or-later` |
| [`psqn`](psqn/) | Performs partially separable quasi-Newton optimization. | `0.3.2` · `Apache-2.0` |
| [`PWEV`](PWEV/) | Builds weighted ensembles for volatility modelling using particle-swarm optimization. | `0.1.0` · `GPL-3.0-only` |
| [`qap`](qap/) | Provides heuristics for the quadratic assignment problem. | `0.1-2` · `GPL-3.0-only` |
| [`QCSIS`](QCSIS/) | Performs sure-independence screening using quantile correlation. | `0.1` · `GPL-2.0-only` |
| [`qrmtools`](qrmtools/) | Provides quantitative risk-management tools and distribution calculations. | `0.0-19` · `GPL-3.0-or-later` |
| [`quadform`](quadform/) | Evaluates quadratic forms efficiently. | `0.0-4` · `GPL` |
| [`quadprog`](quadprog/) | Solves quadratic programs with the Goldfarb-Idnani method. | `1.5-8` · `GPL-2.0-or-later` |
| [`quadprogXT`](quadprogXT/) | Solves quadratic programs with absolute-value constraints. | `0.0.6` · `GPL-2.0-or-later` |
| [`QuantBondCurves`](QuantBondCurves/) | Values bonds and swaps and calibrates interest-rate curves. | `0.3.3` · `GPL-3.0-or-later` |
| [`quantreg`](quantreg/) | Fits quantile-regression models. | `6.1` · `GPL-2.0-or-later` |
| [`quarks`](quarks/) | Calculates and backtests value at risk and expected shortfall. | `1.1.6` · `GPL-3.0-only` |
| [`R4GoodPersonalFinances`](R4GoodPersonalFinances/) | Supports optimization of personal financial decisions. | `1.2.0.9000` · `MIT` |
| [`ragtop`](ragtop/) | Prices equity derivatives with extensions of Black-Scholes. | `2.0.0` · `GPL-2.0-or-later` |
| [`randtoolbox`](randtoolbox/) | Provides pseudo-random and quasi-random number generators. | `2.0.5` · `BSD-3-Clause AND LicenseRef-WELL-Upstream` |
| [`rangen`](rangen/) | Provides fast random-number and matrix generators and sampling utilities. | `0.0.1` · `GPL-3.0-only` |
| [`RCEIM`](RCEIM/) | Performs optimization using a cross-entropy-inspired method. | `0.3` · `GPL-2.0-or-later` |
| [`RcppNumerical`](RcppNumerical/) | Provides numerical optimization, integration, and linear-algebra algorithms. | `0.7-0` · `GPL-2.0-or-later` |
| [`Rcsdp`](Rcsdp/) | Solves semidefinite programs using CSDP algorithms. | `0.1.57.6` · `CPL-1.0` |
| [`Rdsdp`](Rdsdp/) | Solves semidefinite programs using DSDP algorithms. | `1.0.6` · `GPL-3.0-only AND LicenseRef-DSDP AND Apache-2.0` |
| [`REBayes`](REBayes/) | Provides empirical-Bayes and nonparametric maximum-likelihood methods. | `2.60` · `GPL-2.0-or-later` |
| [`relsurv`](relsurv/) | Performs relative-survival analysis. | `2.3-3` · `GPL-2.0-or-later` |
| [`REN`](REN/) | Uses regularization ensembles for robust portfolio optimization. | `0.1.0` · `AGPL-3.0-or-later` |
| [`Rfast`](Rfast/) | Provides efficient statistical, mathematical, and data-analysis routines. | `2.1.5.2` · `GPL-2.0-or-later` |
| [`Rfast2`](Rfast2/) | Provides a second collection of efficient statistical and mathematical functions. | `0.1.5.6` · `GPL-2.0-or-later` |
| [`R-fixedincome`](R-fixedincome/) | Provides fixed-income rates, curves, interpolation, and Nelson-Siegel models. | `0.0.5` · `MIT` |
| [`rgenoud`](rgenoud/) | Performs genetic optimization with optional derivatives. | `5.9-0.3` · `GPL-3.0-only` |
| [`Risk`](Risk/) | Computes financial risk measures for continuous distributions. | `1.0` · `GPL-2.0-or-later` |
| [`riskParityPortfolio`](riskParityPortfolio/) | Designs risk-parity and risk-budgeting portfolios. | `0.2.2.9000` · `GPL-3.0-only` |
| [`RiskPortfolios`](RiskPortfolios/) | Constructs portfolios using risk-based allocation methods. | `2.1.7` · `GPL-2.0-or-later` |
| [`risksimul`](risksimul/) | Simulates rare portfolio losses under t-copulas with t or generalized-hyperbolic marginals. | `0.1.2` · `GPL-2.0-only OR GPL-3.0-only` |
| [`RM2006`](RM2006/) | Estimates conditional covariance using the RiskMetrics 2006 methodology. | `0.1.1` · `GPL-2.0-or-later` |
| [`Rmalschains`](Rmalschains/) | Performs continuous optimization with memetic algorithms and local search chains. | `0.2-11` · `GPL-3.0-only` |
| [`rmgarch`](rmgarch/) | Fits and analyzes multivariate GARCH models. | `1.4-2` · `GPL-3.0-only` |
| [`RMKdiscrete`](RMKdiscrete/) | Provides discrete probability distributions and helper functions. | `0.1` · `GPL-2.0-or-later` |
| [`rmoo`](rmoo/) | Performs multi-objective optimization. | `0.3.2` · `GPL-2.0-or-later` |
| [`rmutil`](rmutil/) | Provides utilities for nonlinear regression and repeated measurements. | `1.1.10` · `GPL-2.0-or-later` |
| [`Rnanoflann`](Rnanoflann/) | Performs nearest-neighbor searches with multiple distance measures. | `0.0.3` · `GPL-3.0-or-later` |
| [`rnd`](rnd/) | Extracts option-implied risk-neutral densities. | `1.2` · `GPL-2.0-or-later` |
| [`rngWELL`](rngWELL/) | Provides WELL-family pseudorandom-number generators. | `0.10-10` · `BSD-3-Clause AND LicenseRef-WELL-Upstream` |
| [`RobStatTM`](RobStatTM/) | Provides robust statistical estimators and tests. | `1.0.11` · `GPL-3.0-or-later` |
| [`robustbase`](robustbase/) | Provides fundamental robust statistical methods. | `0.99-7` · `GPL-2.0-or-later` |
| [`ROI.plugin.qpoases`](ROI.plugin.qpoases/) | Provides qpOASES quadratic programming methods for optimization infrastructure. | `1.0-3` · `GPL-3.0-only` |
| [`rootSolve`](rootSolve/) | Finds nonlinear roots, equilibria, and steady states of ordinary differential equations. | `1.8.2.4` · `GPL-2.0-or-later` |
| [`roptim`](roptim/) | Provides Nelder-Mead, BFGS, conjugate-gradient, L-BFGS-B, and simulated-annealing optimization. | `0.1.7` · `GPL-2.0-or-later` |
| [`RPEGLMEN`](RPEGLMEN/) | Fits Gamma and exponential generalized linear models with elastic-net regularization. | `1.1.4` · `GPL-2.0-or-later` |
| [`RPEIF`](RPEIF/) | Computes influence functions for risk and performance measures. | `1.2.5` · `GPL-3.0-or-later` |
| [`RPESE`](RPESE/) | Estimates standard errors for risk and performance measures. | `1.2.7` · `GPL-3.0-or-later` |
| [`rquantlib`](rquantlib/) | Implements option, bond, and fixed-income pricing algorithms. | `0.4.28` · `GPL-2.0-or-later` |
| [`rrcov`](rrcov/) | Provides robust multivariate estimators with high breakdown points. | `1.7-8` · `GPL-3.0-or-later` |
| [`RSDC`](RSDC/) | Fits regime-switching dynamic-correlation models. | `1.7-0` · `GPL-3.0-only` |
| [`Rsolnp`](Rsolnp/) | Performs constrained nonlinear optimization with an augmented Lagrangian. | `2.0.1` · `GPL-2.0-only` |
| [`RSpectra`](RSpectra/) | Computes selected eigenvalues, eigenvectors, and singular values of large matrices. | `0.16-2` · `MPL-2.0` |
| [`rtl`](rtl/) | Provides trading, risk, and analytics tools for commodities. | `1.3.9` · `MIT` |
| [`rugarch`](rugarch/) | Fits and analyzes univariate GARCH models. | `1.5-6` · `GPL-3.0-only` |
| [`rumidas`](rumidas/) | Fits univariate and double-asymmetric GARCH-MIDAS models. | `0.1.3` · `GPL-3.0-only` |
| [`Runuran`](Runuran/) | Generates random variates using UNU.RAN methods. | `0.41` · `GPL-2.0-or-later` |
| [`rvinecopulib`](rvinecopulib/) | Models pair copulas and D-vines. | `0.7.3.1.0` · `GPL-3.0-only` |
| [`rvmf`](rvmf/) | Generates von Mises-Fisher distributed random vectors. | `0.1.2` · `GPL-3.0-or-later` |
| [`SACCR`](SACCR/) | Calculates counterparty credit risk under the standardized approach. | `3.4` · `GPL-3.0-only` |
| [`sadists`](sadists/) | Provides additional probability distributions. | `0.2.6` · `LGPL-3.0-or-later` |
| [`sandwich`](sandwich/) | Computes robust sandwich covariance-matrix estimators. | `3.1-2` · `GPL-2.0-only OR GPL-3.0-only` |
| [`scip`](scip/) | Provides an interface to the SCIP optimization suite. | `1.10.0-3` · `Apache-2.0` |
| [`scs`](scs/) | Solves convex cone problems using the splitting conic solver. | `3.2.7` · `GPL-3.0-only` |
| [`sde`](sde/) | Simulates and performs inference for stochastic differential equations. | `2.0.21` · `GPL-2.0-or-later` |
| [`segmented`](segmented/) | Fits regression models with estimated breakpoints and changepoints. | `2.2-1` · `GPL-2.0-or-later` |
| [`sgt`](sgt/) | Provides the skewed generalized-t distribution and parameter fitting. | `2.0` · `GPL-3.0-or-later` |
| [`sharper`](sharper/) | Evaluates the statistical significance of Sharpe ratios. | `1.4.0` · `LGPL-3.0-or-later` |
| [`sharpeRratio`](sharpeRratio/) | Estimates Sharpe and signal-to-noise ratios without moment assumptions. | `1.4.3` · `GPL-3.0-only` |
| [`skellam`](skellam/) | Provides the Skellam distribution, estimation, and regression. | `0.2.4` · `GPL-2.0-or-later` |
| [`SkewHyperbolic`](SkewHyperbolic/) | Provides the skew hyperbolic Student-t distribution. | `0.4-2` · `GPL-2.0-or-later` |
| [`skewt`](skewt/) | Provides the skewed Student-t distribution. | `1.0` · `GPL` |
| [`skewunit`](skewunit/) | Estimates and analyzes skew-unit models. | `1.1` · `GPL-2.0-or-later` |
| [`slam`](slam/) | Provides sparse lightweight arrays and matrices. | `0.1-56` · `GPL-2.0-only` |
| [`SmithWilsonYieldCurve`](SmithWilsonYieldCurve/) | Constructs yield curves using the Smith-Wilson method. | `1.1.1` · `GPL-3.0-only` |
| [`smoof`](smoof/) | Provides single- and multi-objective optimization test functions. | `1.7.0` · `BSD-2-Clause` |
| [`smoots`](smoots/) | Smooths trends and forecasts equidistant time series. | `1.1.4` · `GPL-3.0-only` |
| [`sn`](sn/) | Provides skew-normal, skew-t, skew-Cauchy, multivariate, and SUN distributions. | `2.1.3` · `GPL-2.0-only OR GPL-3.0-only` |
| [`sna`](sna/) | Provides tools for social-network analysis. | `2.8` · `GPL-2.0-or-later` |
| [`soma`](soma/) | Performs general-purpose optimization with self-organizing migrating algorithms. | `1.2.0` · `GPL-2.0-only` |
| [`spacefillr`](spacefillr/) | Generates space-filling Halton and Sobol sequences. | `0.4.0` · `MIT` |
| [`spantest`](spantest/) | Performs mean-variance spanning tests. | `1.4-0` · `GPL-3.0-only` |
| [`sparseIndexTracking`](sparseIndexTracking/) | Designs sparse portfolios that track a financial index. | `0.1.1` · `GPL-3.0-only` |
| [`spectralGraphTopology`](spectralGraphTopology/) | Learns graph topologies from data using spectral constraints. | `0.2.3` · `GPL-3.0-only` |
| [`Spillover`](Spillover/) | Computes VAR-based generalized and orthogonalized connectedness measures. | `0.1.1` · `GPL-2.0-only` |
| [`splines`](splines/) | Provides B-spline and natural-spline basis calculations. | `2.0-7` · `GPL-2.0-or-later` |
| [`stochfactor`](stochfactor/) | Models univariate and factor stochastic volatility. | `stochvol 3.2.9 + factorstochvol 1.1.2` · `GPL-2.0-or-later` |
| [`stochQN`](stochQN/) | Provides stochastic limited-memory quasi-Newton optimizers. | `0.1.2-1` · `BSD-2-Clause` |
| [`stochvolTMB`](stochvolTMB/) | Estimates stochastic-volatility models by Laplace maximum likelihood. | `0.3.0` · `GPL-3.0-only` |
| [`stockAnalyst`](stockAnalyst/) | Provides equity-valuation, return, growth, and required-return calculations. | `1.0.1` · `GPL-3.0-only` |
| [`strand`](strand/) | Provides a framework for investment-strategy simulation. | `0.2.3` · `GPL-3.0-only` |
| [`StReg`](StReg/) | Fits static and dynamic Student-t regression models. | `1.1` · `GPL-2.0-only` |
| [`subplex`](subplex/) | Performs unconstrained optimization using the Subplex algorithm. | `1.9` · `GPL-3.0-only` |
| [`SuppDists`](SuppDists/) | Provides supplemental probability distributions and rank-statistic calculations. | `1.1-9.9` · `GPL-2.0-or-later` |
| [`survey`](survey/) | Analyzes complex survey samples. | `4.5` · `GPL-2.0-only OR GPL-3.0-only` |
| [`survival`](survival/) | Provides core methods for survival analysis. | `3.8-9` · `GPL-2.0-or-later` |
| [`svdnf`](svdnf/) | Performs discrete nonlinear filtering for stochastic-volatility models. | `0.1.11` · `GPL-3.0-only` |
| [`tabuSearch`](tabuSearch/) | Performs tabu search over binary configurations. | `1.2.0` · `GPL-2.0-or-later` |
| [`tensorA`](tensorA/) | Provides advanced tensor arithmetic with named indices. | `0.36.2.1` · `GPL-2.0-or-later` |
| [`timsac`](timsac/) | Provides time-series analysis, prediction, and control methods. | `1.3.8-6` · `GPL-2.0-or-later` |
| [`tolerance`](tolerance/) | Computes statistical tolerance intervals and regions. | `3.0.0` · `GPL-2.0-or-later` |
| [`Trading`](Trading/) | Provides trading, correlation, beta, and betting calculations. | `3.2` · `GPL-3.0-only` |
| [`trawl`](trawl/) | Estimates and simulates trawl processes. | `0.2.2` · `GPL-3.0-only` |
| [`treasuryTR`](treasuryTR/) | Generates Treasury total returns from yield data. | `0.1.6` · `MIT` |
| [`trust`](trust/) | Performs trust-region optimization. | `0.1-9` · `MIT` |
| [`trustOptim`](trustOptim/) | Performs trust-region optimization for functions with sparse Hessians. | `0.8.7.4` · `MPL-2.0` |
| [`tsa`](tsa/) | Provides time-series analysis methods. | `1.3.1` · `GPL-2.0-only OR GPL-3.0-only` |
| [`tscopula`](tscopula/) | Fits and simulates time-series copula models. | `0.3.9` · `GPL-3.0-only` |
| [`tsdistributions`](tsdistributions/) | Provides standardized distributions for time-series modelling. | `1.0.4` · `GPL-2.0-only` |
| [`tsdyn`](tsdyn/) | Fits nonlinear time-series models with regime switching. | `11.0.5.2` · `GPL-2.0-or-later` |
| [`tseries`](tseries/) | Provides time-series analysis and computational-finance methods. | `0.10-62` · `GPL-2.0-only OR GPL-3.0-only` |
| [`tserieschaos`](tserieschaos/) | Analyzes nonlinear and chaotic time series. | `0.1-13.1` · `GPL-2.0-only` |
| [`tsgarch`](tsgarch/) | Fits, simulates, forecasts, and diagnoses univariate GARCH models. | `1.0.4` · `GPL-2.0-only` |
| [`tsmarch`](tsmarch/) | Fits and forecasts multivariate ARCH models. | `1.0.3` · `GPL-2.0-only` |
| [`TSP`](TSP/) | Provides algorithms and infrastructure for traveling-salesperson problems. | `1.2.7` · `GPL-3.0-only` |
| [`tvgarch`](tvgarch/) | Fits time-varying GARCH models. | `2.4.3` · `GPL-2.0-or-later` |
| [`tvGarchKF`](tvGarchKF/) | Fits time-varying GARCH models through a state-space representation. | `0.0.1` · `GPL-3.0-or-later` |
| [`tvm`](tvm/) | Performs time-value-of-money and interest-rate calculations. | `0.5.2` · `MIT` |
| [`TVMVP`](TVMVP/) | Estimates time-varying covariance matrices for portfolio optimization. | `1.0.5` · `MIT` |
| [`ufRisk`](ufRisk/) | Estimates GARCH-based value at risk and expected shortfall and backtests the forecasts. | `1.0.7` · `GPL-3.0-only` |
| [`uncorbets`](uncorbets/) | Constructs uncorrelated bets using the minimum-torsion algorithm. | `0.1.2.9000` · `MIT` |
| [`vamc`](vamc/) | Values variable annuities by Monte Carlo simulation. | `0.2.1` · `GPL-2.0-only` |
| [`vares`](vares/) | Computes parametric value at risk and expected shortfall. | `1.0.2` · `GPL-2.0-or-later` |
| [`vasicekfit`](vasicekfit/) | Fits the extended Vasicek credit-loss model. | `0.2.0` · `MIT` |
| [`VGAM`](VGAM/) | Fits vector generalized linear and additive models. | `1.1-14` · `GPL-3.0-only` |
| [`vrtest`](vrtest/) | Performs variance-ratio and martingale-difference tests. | `1.2` · `GPL-2.0-only` |
| [`waveslim`](waveslim/) | Provides one-, two-, and three-dimensional wavelet methods. | `1.8.5` · `BSD-3-Clause` |
| [`wqc`](wqc/) | Performs wavelet quantile-correlation analysis. | `0.1.2` · `GPL-3.0-only` |
| [`xVA`](xVA/) | Computes credit-risk valuation adjustments. | `1.3` · `GPL-3.0-only` |
| [`ycevo`](ycevo/) | Estimates yield-curve evolution nonparametrically. | `0.2.1.9000` · `GPL-3.0-only` |
| [`YieldCurve`](YieldCurve/) | Models and estimates yield curves. | `5.1` · `GPL-2.0-or-later` |
| [`yieldcurves`](yieldcurves/) | Fits, analyzes, and decomposes yield curves. | `0.1.0` · `MIT` |
| [`yrnd`](yrnd/) | Extracts risk-neutral densities from financial prices and rates. | `0.1.5` · `GPL-3.0-only` |
| [`ZeroOneDists`](ZeroOneDists/) | Provides statistical distributions with support on the zero-one interval. | `1.0.0` · `MIT` |
| [`zigg`](zigg/) | Generates pseudorandom numbers using the Ziggurat method. | `0.0.2` · `GPL-2.0-or-later` |

The repository is an aggregate of separately licensed projects; there is no
additional repository-wide license. Consult the license/copying, notice,
origin, and README files within each package before redistributing or modifying
it.

## Attribution and provenance

Each project identifies the original package, version, authors or copyright
holders, and the nature of the Fortran work. Original package metadata is
retained under package `reference/` directories where available. The modern
Fortran files also carry SPDX license identifiers and dated modification
notices.
