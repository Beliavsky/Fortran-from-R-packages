# Fortran from R packages

This repository collects experimental modern Fortran translations and ports of
computational code in 0 package directories derived from R packages. Each
subdirectory is an independent
Fortran Package Manager (fpm) project with its own documentation, tests,
provenance record, and license.

These projects are unofficial and are not endorsed by the original package
authors, CRAN, or the R Foundation. They are experimental and incompletely
validated; independently verify numerical results before relying on them.

## Packages and licenses

| Package | What it does | Upstream version · License |
| --- | --- | --- |
| [`backtest`](backtest/) | Explores portfolio-based conjectures about financial instruments. | `0.3-4` · `GPL-2.0-or-later` |
| [`bayesgarch`](bayesgarch/) | Performs Bayesian estimation of GARCH models with Student-t innovations. | `2.1.10` · `GPL-2.0-or-later` |
| [`bcc1997`](bcc1997/) | Prices European options with stochastic volatility, rates, and jumps. | `0.1.1` · `GPL-2.0-or-later` |
| [`betategarch`](betategarch/) | Estimates Beta-t-EGARCH volatility models. | `3.4` · `GPL-2.0-only` |
| [`bidask`](bidask/) | Estimates bid-ask spreads efficiently from OHLC prices. | `2.1.5` · `MIT` |
| [`blmodel`](blmodel/) | Computes Black-Litterman posterior distributions. | `1.0.2` · `GPL-3.0-only` |
| [`deoptimr`](deoptimr/) | Performs global optimization using differential evolution. | `1.2-0` · `GPL-2.0-or-later` |
| [`derivmkts`](derivmkts/) | Provides derivative pricing and financial-market calculations. | `0.2.5.1` · `MIT` |
| [`epo`](epo/) | Performs enhanced portfolio optimization with correlation shrinkage. | `0.1.0.9000` · `MIT` |
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
| [`garchsk`](garchsk/) | Estimates GARCH models with conditional skewness and kurtosis. | `0.1.0` · `GPL-2.0-or-later` |
| [`garchx`](garchx/) | Fits GARCH models with exogenous covariates. | `1.7` · `GPL-2.0-or-later` |
| [`gogarch`](gogarch/) | Fits generalized orthogonal GARCH models. | `0.7-6` · `GPL-2.0-or-later` |
| [`greeks`](greeks/) | Computes option sensitivities, implied volatility, and Monte Carlo Greeks. | `1.5.6` · `MIT` |
| [`hdshop`](hdshop/) | Constructs high-dimensional shrinkage optimal portfolios. | `0.1.7` · `GPL-3.0-only` |
| [`highfrequency`](highfrequency/) | Analyzes high-frequency trade and quote data. | `1.0.2` · `GPL-2.0-or-later` |
| [`lgarch`](lgarch/) | Simulates and estimates log-GARCH models. | `0.7` · `GPL-2.0-only` |
| [`longmemo`](longmemo/) | Provides statistical methods for long-memory processes. | `1.1-4` · `GPL-2.0-or-later` |
| [`lsmontecarlo`](lsmontecarlo/) | Prices American options using least-squares Monte Carlo. | `1.0` · `GPL-3.0-only` |
| [`markowitzr`](markowitzr/) | Performs statistical inference for Markowitz portfolios. | `1.0.2.0002` · `LGPL-3.0-or-later` |
| [`msgarch`](msgarch/) | Fits Markov-switching GARCH models. | `2.51` · `GPL-2.0-or-later` |
| [`ob-analytics`](ob-analytics/) | Analyzes limit order books and liquidity. | `0.1.2` · `GPL-2.0-or-later` |
| [`opthedging`](opthedging/) | Values and optimally hedges call and put options. | `1.0` · `GPL-2.0-or-later` |
| [`optionpricing`](optionpricing/) | Prices options with efficient simulation algorithms. | `0.1.2` · `GPL-2.0-only OR GPL-3.0-only` |
| [`pa`](pa/) | Performs equity-portfolio performance attribution. | `1.2-4` · `GPL-2.0-only` |
| [`parma`](parma/) | Provides portfolio allocation and risk-management applications. | `1.7` · `GPL-3.0-or-later` |
| [`pbo`](pbo/) | Estimates the probability of backtest overfitting. | `1.3.5` · `MIT` |
| [`peerperformance`](peerperformance/) | Performs luck-corrected peer-performance analysis. | `2.4.0` · `GPL-2.0-or-later` |
| [`performanceanalytics`](performanceanalytics/) | Provides econometric tools for performance and risk analysis. | `2.1.0` · `GPL-2.0-or-later` |
| [`portfoliooptim`](portfoliooptim/) | Performs small- and large-sample portfolio optimization. | `1.1.1` · `GPL-3.0-only` |
| [`rmgarch`](rmgarch/) | Fits and analyzes multivariate GARCH models. | `1.4-2` · `GPL-3.0-only` |
| [`rnd`](rnd/) | Extracts option-implied risk-neutral densities. | `1.2` · `GPL-2.0-or-later` |
| [`robustbase`](robustbase/) | Provides fundamental robust statistical methods. | `0.99-7` · `GPL-2.0-or-later` |
| [`rquantlib`](rquantlib/) | Implements option, bond, and fixed-income pricing algorithms. | `0.4.28` · `GPL-2.0-or-later` |
| [`rtl`](rtl/) | Provides trading, risk, and analytics tools for commodities. | `1.3.9` · `MIT` |
| [`rugarch`](rugarch/) | Fits and analyzes univariate GARCH models. | `1.5-6` · `GPL-3.0-only` |
| [`sde`](sde/) | Simulates and performs inference for stochastic differential equations. | `2.0.21` · `GPL-2.0-or-later` |
| [`sharper`](sharper/) | Evaluates the statistical significance of Sharpe ratios. | `1.4.0` · `LGPL-3.0-or-later` |
| [`stochfactor`](stochfactor/) | Models univariate and factor stochastic volatility. | `stochvol 3.2.9 + factorstochvol 1.1.2` · `GPL-2.0-or-later` |
| [`strand`](strand/) | Provides a framework for investment-strategy simulation. | `0.2.3` · `GPL-3.0-only` |
| [`svdnf`](svdnf/) | Performs discrete nonlinear filtering for stochastic-volatility models. | `0.1.11` · `GPL-3.0-only` |
| [`timsac`](timsac/) | Provides time-series analysis, prediction, and control methods. | `1.3.8-6` · `GPL-2.0-or-later` |
| [`tsdyn`](tsdyn/) | Fits nonlinear time-series models with regime switching. | `11.0.5.2` · `GPL-2.0-or-later` |
| [`tseries`](tseries/) | Provides time-series analysis and computational-finance methods. | `0.10-62` · `GPL-2.0-only OR GPL-3.0-only` |
| [`tserieschaos`](tserieschaos/) | Analyzes nonlinear and chaotic time series. | `0.1-13.1` · `GPL-2.0-only` |
| [`tvm`](tvm/) | Performs time-value-of-money and interest-rate calculations. | `0.5.2` · `MIT` |
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
