# Fortran from R packages

This repository collects experimental modern Fortran translations and ports of
computational code in 89 package directories derived from R packages. Each
subdirectory is an independent
Fortran Package Manager (fpm) project with its own documentation, tests,
provenance record, and license.

These projects are unofficial and are not endorsed by the original package
authors, CRAN, or the R Foundation. They are experimental and incompletely
validated; independently verify numerical results before relying on them.

## Packages and licenses

| Package | What it does | Upstream version · License |
| --- | --- | --- |
| [`ACDm`](ACDm/) | Estimates and simulates autoregressive conditional duration models. | `1.1.0` · `GPL-3.0-or-later` |
| [`actuar`](actuar/) | Provides actuarial distributions, loss models, credibility methods, and risk calculations. | `3.3-7` · `GPL-2.0-or-later` |
| [`apt`](apt/) | Models asymmetric price transmission with threshold cointegration and error-correction methods. | `4.0` · `GPL-2.0-or-later` |
| [`backtest`](backtest/) | Explores portfolio-based conjectures about financial instruments. | `0.3-4` · `GPL-2.0-or-later` |
| [`bayesgarch`](bayesgarch/) | Performs Bayesian estimation of GARCH models with Student-t innovations. | `2.1.10` · `GPL-2.0-or-later` |
| [`bayesianOU`](bayesianOU/) | Fits Bayesian nonlinear Ornstein-Uhlenbeck models. | `0.2.0` · `MIT` |
| [`bcc1997`](bcc1997/) | Prices European options with stochastic volatility, rates, and jumps. | `0.1.1` · `GPL-2.0-or-later` |
| [`BEKKs`](BEKKs/) | Estimates and analyzes BEKK multivariate conditional-volatility models. | `1.4.7` · `MIT` |
| [`betategarch`](betategarch/) | Estimates Beta-t-EGARCH volatility models. | `3.4` · `GPL-2.0-only` |
| [`bidask`](bidask/) | Estimates bid-ask spreads efficiently from OHLC prices. | `2.1.5` · `MIT` |
| [`blmodel`](blmodel/) | Computes Black-Litterman posterior distributions. | `1.0.2` · `GPL-3.0-only` |
| [`bondAnalyst`](bondAnalyst/) | Performs fixed-income valuation and yield, spread, and duration calculations. | `1.0.1` · `GPL-3.0-only` |
| [`cccp`](cccp/) | Solves cone-constrained convex optimization problems. | `0.3-3` · `GPL-3.0-or-later` |
| [`copula`](copula/) | Models multivariate dependence with common copula families. | `1.1-7` · `GPL-3.0-or-later` |
| [`cvar`](cvar/) | Computes value at risk and expected shortfall from distributions or samples. | `0.6` · `GPL-2.0-or-later` |
| [`deoptimr`](deoptimr/) | Performs global optimization using differential evolution. | `1.2-0` · `GPL-2.0-or-later` |
| [`derivmkts`](derivmkts/) | Provides derivative pricing and financial-market calculations. | `0.2.5.1` · `MIT` |
| [`Dowd`](Dowd/) | Provides quantitative financial risk-management calculations. | `0.12` · `GPL-2.0-only OR GPL-3.0-only` |
| [`epo`](epo/) | Performs enhanced portfolio optimization with correlation shrinkage. | `0.1.0.9000` · `MIT` |
| [`etrm`](etrm/) | Values and analyzes energy-trading and risk-management instruments. | `1.0.2` · `MIT` |
| [`evir`](evir/) | Performs extreme-value analysis for financial risk. | `1.7-4` · `GPL-2.0-or-later` |
| [`fattailsr`](fattailsr/) | Provides Kiener distributions and fat-tail analytics. | `2.0.1` · `GPL-2.0-only` |
| [`fbasics`](fbasics/) | Provides financial-market statistics, distributions, and utilities. | `4052.98` · `GPL-2.0-or-later` |
| [`fbonds`](fbonds/) | Prices bonds and fits Nelson-Siegel family term structures. | `3042.78` · `GPL-2.0-or-later` |
| [`fcopulae`](fcopulae/) | Models dependence with elliptical, Archimedean, and empirical copulas. | `4052.86` · `GPL-2.0-or-later` |
| [`fextremes`](fextremes/) | Models extreme values and financial tail risk. | `4032.84` · `GPL-2.0-or-later` |
| [`fGarch`](fGarch/) | Fits and analyzes univariate GARCH and APARCH models. | `4052.93` · `GPL-2.0-or-later` |
| [`fhmm`](fhmm/) | Fits hidden Markov models to financial data. | `1.4.3` · `GPL-3.0-only` |
| [`financialmath`](financialmath/) | Provides financial mathematics for actuaries. | `0.1.1` · `GPL-2.0-only` |
| [`fmultivar`](fmultivar/) | Provides multivariate distributions and financial-data analysis. | `4031.84` · `GPL-2.0-or-later` |
| [`fnonlinear`](fnonlinear/) | Models nonlinear and chaotic time series. | `4052.83` · `GPL-2.0-or-later` |
| [`fportfolio`](fportfolio/) | Performs portfolio selection, optimization, risk analysis, and backtesting. | `4023.84` · `GPL-2.0-or-later` |
| [`fracdiff`](fracdiff/) | Estimates, simulates, and analyzes fractionally differenced time-series models. | `1.5-4` · `GPL-2.0-or-later` |
| [`frapo`](frapo/) | Provides financial risk modelling and portfolio optimization methods. | `0.4-2` · `GPL-3.0-or-later` |
| [`garchsk`](garchsk/) | Estimates GARCH models with conditional skewness and kurtosis. | `0.1.0` · `GPL-2.0-or-later` |
| [`garchx`](garchx/) | Fits GARCH models with exogenous covariates. | `1.7` · `GPL-2.0-or-later` |
| [`gcpm`](gcpm/) | Models credit-portfolio risk analytically and by Monte Carlo simulation. | `1.2.2` · `GPL-2.0-only` |
| [`ghyp`](ghyp/) | Evaluates, fits, and simulates generalized hyperbolic distributions. | `1.6.5` · `GPL-2.0-or-later` |
| [`gogarch`](gogarch/) | Fits generalized orthogonal GARCH models. | `0.7-6` · `GPL-2.0-or-later` |
| [`greeks`](greeks/) | Computes option sensitivities, implied volatility, and Monte Carlo Greeks. | `1.5.6` · `MIT` |
| [`hdshop`](hdshop/) | Constructs high-dimensional shrinkage optimal portfolios. | `0.1.7` · `GPL-3.0-only` |
| [`highfrequency`](highfrequency/) | Analyzes high-frequency trade and quote data. | `1.0.2` · `GPL-2.0-or-later` |
| [`lgarch`](lgarch/) | Simulates and estimates log-GARCH models. | `0.7` · `GPL-2.0-only` |
| [`longmemo`](longmemo/) | Provides statistical methods for long-memory processes. | `1.1-4` · `GPL-2.0-or-later` |
| [`lsmontecarlo`](lsmontecarlo/) | Prices American options using least-squares Monte Carlo. | `1.0` · `GPL-3.0-only` |
| [`markowitzr`](markowitzr/) | Performs statistical inference for Markowitz portfolios. | `1.0.2.0002` · `LGPL-3.0-or-later` |
| [`mcrp`](mcrp/) | Constructs multiple-criteria risk-parity portfolios. | `0.0-1` · `GPL-3.0-only` |
| [`msgarch`](msgarch/) | Fits Markov-switching GARCH models. | `2.51` · `GPL-2.0-or-later` |
| [`mvtnorm`](mvtnorm/) | Computes multivariate normal and Student-t probabilities, densities, and random samples. | `1.4-2` · `GPL-2.0-only` |
| [`nmof`](nmof/) | Provides numerical optimization methods for finance and economics. | `2.12-0` · `GPL-3.0-only` |
| [`nvmix`](nvmix/) | Computes and simulates multivariate normal variance-mixture distributions. | `0.1-2` · `GPL-3.0-or-later` |
| [`ob-analytics`](ob-analytics/) | Analyzes limit order books and liquidity. | `0.1.2` · `GPL-2.0-or-later` |
| [`opthedging`](opthedging/) | Values and optimally hedges call and put options. | `1.0` · `GPL-2.0-or-later` |
| [`optionpricing`](optionpricing/) | Prices options with efficient simulation algorithms. | `0.1.2` · `GPL-2.0-only OR GPL-3.0-only` |
| [`pa`](pa/) | Performs equity-portfolio performance attribution. | `1.2-4` · `GPL-2.0-only` |
| [`parma`](parma/) | Provides portfolio allocation and risk-management applications. | `1.7` · `GPL-3.0-or-later` |
| [`pbo`](pbo/) | Estimates the probability of backtest overfitting. | `1.3.5` · `MIT` |
| [`peerperformance`](peerperformance/) | Performs luck-corrected peer-performance analysis. | `2.4.0` · `GPL-2.0-or-later` |
| [`performanceanalytics`](performanceanalytics/) | Provides econometric tools for performance and risk analysis. | `2.1.0` · `GPL-2.0-or-later` |
| [`portfoliooptim`](portfoliooptim/) | Performs small- and large-sample portfolio optimization. | `1.1.1` · `GPL-3.0-only` |
| [`qrmtools`](qrmtools/) | Provides quantitative risk-management tools and distribution calculations. | `0.0-19` · `GPL-3.0-or-later` |
| [`Risk`](Risk/) | Computes financial risk measures for continuous distributions. | `1.0` · `GPL-2.0-or-later` |
| [`riskParityPortfolio`](riskParityPortfolio/) | Designs risk-parity and risk-budgeting portfolios. | `0.2.2.9000` · `GPL-3.0-only` |
| [`RiskPortfolios`](RiskPortfolios/) | Constructs portfolios using risk-based allocation methods. | `2.1.7` · `GPL-2.0-or-later` |
| [`risksimul`](risksimul/) | Simulates rare portfolio losses under t-copulas with t or generalized-hyperbolic marginals. | `0.1.2` · `GPL-2.0-only OR GPL-3.0-only` |
| [`RM2006`](RM2006/) | Estimates conditional covariance using the RiskMetrics 2006 methodology. | `0.1.1` · `GPL-2.0-or-later` |
| [`rmgarch`](rmgarch/) | Fits and analyzes multivariate GARCH models. | `1.4-2` · `GPL-3.0-only` |
| [`rnd`](rnd/) | Extracts option-implied risk-neutral densities. | `1.2` · `GPL-2.0-or-later` |
| [`robustbase`](robustbase/) | Provides fundamental robust statistical methods. | `0.99-7` · `GPL-2.0-or-later` |
| [`rquantlib`](rquantlib/) | Implements option, bond, and fixed-income pricing algorithms. | `0.4.28` · `GPL-2.0-or-later` |
| [`rtl`](rtl/) | Provides trading, risk, and analytics tools for commodities. | `1.3.9` · `MIT` |
| [`rugarch`](rugarch/) | Fits and analyzes univariate GARCH models. | `1.5-6` · `GPL-3.0-only` |
| [`sde`](sde/) | Simulates and performs inference for stochastic differential equations. | `2.0.21` · `GPL-2.0-or-later` |
| [`sharper`](sharper/) | Evaluates the statistical significance of Sharpe ratios. | `1.4.0` · `LGPL-3.0-or-later` |
| [`smoots`](smoots/) | Smooths trends and forecasts equidistant time series. | `1.1.4` · `GPL-3.0-only` |
| [`stochfactor`](stochfactor/) | Models univariate and factor stochastic volatility. | `stochvol 3.2.9 + factorstochvol 1.1.2` · `GPL-2.0-or-later` |
| [`stockAnalyst`](stockAnalyst/) | Provides equity-valuation, return, growth, and required-return calculations. | `1.0.1` · `GPL-3.0-only` |
| [`strand`](strand/) | Provides a framework for investment-strategy simulation. | `0.2.3` · `GPL-3.0-only` |
| [`svdnf`](svdnf/) | Performs discrete nonlinear filtering for stochastic-volatility models. | `0.1.11` · `GPL-3.0-only` |
| [`timsac`](timsac/) | Provides time-series analysis, prediction, and control methods. | `1.3.8-6` · `GPL-2.0-or-later` |
| [`tsdyn`](tsdyn/) | Fits nonlinear time-series models with regime switching. | `11.0.5.2` · `GPL-2.0-or-later` |
| [`tseries`](tseries/) | Provides time-series analysis and computational-finance methods. | `0.10-62` · `GPL-2.0-only OR GPL-3.0-only` |
| [`tserieschaos`](tserieschaos/) | Analyzes nonlinear and chaotic time series. | `0.1-13.1` · `GPL-2.0-only` |
| [`tvm`](tvm/) | Performs time-value-of-money and interest-rate calculations. | `0.5.2` · `MIT` |
| [`ufRisk`](ufRisk/) | Estimates GARCH-based value at risk and expected shortfall and backtests the forecasts. | `1.0.7` · `GPL-3.0-only` |
| [`vares`](vares/) | Computes parametric value at risk and expected shortfall. | `1.0.2` · `GPL-2.0-or-later` |
| [`vasicekfit`](vasicekfit/) | Fits the extended Vasicek credit-loss model. | `0.2.0` · `MIT` |
| [`vrtest`](vrtest/) | Performs variance-ratio and martingale-difference tests. | `1.2` · `GPL-2.0-only` |
| [`yieldcurves`](yieldcurves/) | Fits, analyzes, and decomposes yield curves. | `0.1.0` · `MIT` |

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
