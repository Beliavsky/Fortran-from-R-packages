# Fortran from R packages

This repository collects experimental modern Fortran translations and ports of
computational code in 58 package directories derived from R packages. Each
subdirectory is an independent
Fortran Package Manager (fpm) project with its own documentation, tests,
provenance record, and license.

These projects are unofficial and are not endorsed by the original package
authors, CRAN, or the R Foundation. They are experimental and incompletely
validated; independently verify numerical results before relying on them.

## Packages and licenses

| Directory | Original R package | License for the Fortran project |
| --- | --- | --- |
| `backtest-modern-fortran/` | backtest 0.3-4 | GPL-2.0-or-later |
| `bayesgarch-modern-fortran/` | bayesGARCH 2.1.10 | GPL-2.0-or-later |
| `bcc1997-fortran/` | BCC1997 0.1.1 | GPL-2.0-or-later |
| `betategarch-modern-fortran/` | betategarch 3.4 | GPL-2.0-only |
| `bidask-fortran/` | bidask 2.1.5 | MIT |
| `blmodel-fortran/` | BLModel 1.0.2 | GPL-3.0-only |
| `deoptimr-modern-fortran/` | DEoptimR 1.2-0 | GPL-2.0-or-later |
| `derivmkts-fortran/` | derivmkts 0.2.5.1 | MIT |
| `epo-fortran/` | epo 0.1.0.9000 | MIT |
| `fbasics-modern-fortran/` | fBasics 4052.98 | GPL-2.0-or-later |
| `fbonds-modern-fortran/` | fBonds 3042.78 | GPL-2.0-or-later |
| `fcopulae-modern-fortran/` | fCopulae 4052.86 | GPL-2.0-or-later |
| `fextremes-modern-fortran/` | fExtremes 4032.84 | GPL-2.0-or-later |
| `fGarch-modern-fortran/` | fGarch 4052.93 | GPL-2.0-or-later |
| `fhmm-fortran/` | fHMM 1.4.3 | GPL-3.0-only |
| `financialmath-fortran/` | FinancialMath 0.1.1 | GPL-2.0-only |
| `fmultivar-modern-fortran/` | fMultivar 4031.84 | GPL-2.0-or-later |
| `fnonlinear-modern-fortran/` | fNonlinear 4052.83 | GPL-2.0-or-later |
| `fportfolio-modern-fortran/` | fPortfolio 4023.84 | GPL-2.0-or-later |
| `garchsk-fortran/` | GARCHSK 0.1.0 | GPL-2.0-or-later |
| `garchx-modern-fortran/` | garchx 1.7 | GPL-2.0-or-later |
| `gogarch-modern-fortran/` | gogarch 0.7-6 | GPL-2.0-or-later |
| `greeks-fortran/` | greeks 1.5.6 | MIT |
| `hdshop-fortran/` | HDShOP 0.1.7 | GPL-3.0-only |
| `highfrequency-fortran/` | highfrequency 1.0.2 | GPL-2.0-or-later |
| `lgarch-modern-fortran/` | lgarch 0.7 | GPL-2.0-only |
| `longmemo-modern-fortran/` | longmemo 1.1-4 | GPL-2.0-or-later |
| `lsmontecarlo-fortran/` | LSMonteCarlo 1.0 | GPL-3.0-only |
| `markowitzr-fortran/` | MarkowitzR 1.0.2.0002 | LGPL-3.0-or-later |
| `msgarch-modern-fortran/` | MSGARCH 2.51 | GPL-2.0-or-later |
| `ob-analytics-fortran/` | obAnalytics 0.1.2 | GPL-2.0-or-later |
| `opthedging-fortran/` | OptHedging 1.0 | GPL-2.0-or-later |
| `optionpricing-fortran/` | OptionPricing 0.1.2 | GPL-2.0-only OR GPL-3.0-only |
| `pa-modern-fortran/` | pa 1.2-4 | GPL-2.0-only |
| `parma-fortran/` | parma 1.7 | GPL-3.0-or-later |
| `pbo-fortran/` | pbo 1.3.5 | MIT |
| `peerperformance-fortran/` | PeerPerformance 2.4.0 | GPL-2.0-or-later |
| `performanceanalytics-modern-fortran/` | PerformanceAnalytics 2.1.0 | GPL-2.0-or-later |
| `portfoliooptim-fortran/` | PortfolioOptim 1.1.1 | GPL-3.0-only |
| `rnd-fortran/` | RND 1.2 | GPL-2.0-or-later |
| `rmgarch-modern-fortran/` | rmgarch 1.4-2 | GPL-3.0-only |
| `rquantlib-modern-fortran/` | RQuantLib 0.4.28 | GPL-2.0-or-later |
| `robustbase-modern-fortran/` | robustbase 0.99-7 | GPL-2.0-or-later |
| `rtl-fortran/` | RTL 1.3.9 | MIT |
| `rugarch-modern-fortran/` | rugarch 1.5-6 | GPL-3.0-only |
| `sde-fortran/` | sde 2.0.21 | GPL-2.0-or-later |
| `sharper-fortran/` | SharpeR 1.4.0 | LGPL-3.0-or-later |
| `stochfactor-modern-fortran/` | stochvol 3.2.9 and factorstochvol 1.1.2 | GPL-2.0-or-later |
| `strand-fortran/` | strand 0.2.3 | GPL-3.0-only |
| `svdnf-fortran/` | SVDNF 0.1.11 | GPL-3.0-only |
| `timsac-modern-fortran/` | timsac 1.3.8-6 | GPL-2.0-or-later |
| `tsdyn-modern-fortran/` | tsDyn 11.0.5.2 | GPL-2.0-or-later |
| `tseries-modern-fortran/` | tseries 0.10-62 | GPL-2.0-only OR GPL-3.0-only |
| `tserieschaos-modern-fortran/` | tseriesChaos 0.1-13.1 | GPL-2.0-only |
| `tvm-fortran/` | tvm 0.5.2 | MIT |
| `vrtest-fortran/` | vrtest 1.2 | GPL-2.0-only |
| `yieldcurves-fortran/` | yieldcurves 0.1.0 | MIT |

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
