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

| Directory | Original R package | License for the Fortran project |
| --- | --- | --- |
| `backtest/` | backtest 0.3-4 | GPL-2.0-or-later |
| `bayesgarch/` | bayesGARCH 2.1.10 | GPL-2.0-or-later |
| `bcc1997/` | BCC1997 0.1.1 | GPL-2.0-or-later |
| `betategarch/` | betategarch 3.4 | GPL-2.0-only |
| `bidask/` | bidask 2.1.5 | MIT |
| `blmodel/` | BLModel 1.0.2 | GPL-3.0-only |
| `deoptimr/` | DEoptimR 1.2-0 | GPL-2.0-or-later |
| `derivmkts/` | derivmkts 0.2.5.1 | MIT |
| `epo/` | epo 0.1.0.9000 | MIT |
| `fbasics/` | fBasics 4052.98 | GPL-2.0-or-later |
| `fbonds/` | fBonds 3042.78 | GPL-2.0-or-later |
| `fcopulae/` | fCopulae 4052.86 | GPL-2.0-or-later |
| `fextremes/` | fExtremes 4032.84 | GPL-2.0-or-later |
| `fGarch/` | fGarch 4052.93 | GPL-2.0-or-later |
| `fhmm/` | fHMM 1.4.3 | GPL-3.0-only |
| `financialmath/` | FinancialMath 0.1.1 | GPL-2.0-only |
| `fmultivar/` | fMultivar 4031.84 | GPL-2.0-or-later |
| `fnonlinear/` | fNonlinear 4052.83 | GPL-2.0-or-later |
| `fportfolio/` | fPortfolio 4023.84 | GPL-2.0-or-later |
| `garchsk/` | GARCHSK 0.1.0 | GPL-2.0-or-later |
| `garchx/` | garchx 1.7 | GPL-2.0-or-later |
| `gogarch/` | gogarch 0.7-6 | GPL-2.0-or-later |
| `greeks/` | greeks 1.5.6 | MIT |
| `hdshop/` | HDShOP 0.1.7 | GPL-3.0-only |
| `highfrequency/` | highfrequency 1.0.2 | GPL-2.0-or-later |
| `lgarch/` | lgarch 0.7 | GPL-2.0-only |
| `longmemo/` | longmemo 1.1-4 | GPL-2.0-or-later |
| `lsmontecarlo/` | LSMonteCarlo 1.0 | GPL-3.0-only |
| `markowitzr/` | MarkowitzR 1.0.2.0002 | LGPL-3.0-or-later |
| `msgarch/` | MSGARCH 2.51 | GPL-2.0-or-later |
| `ob-analytics/` | obAnalytics 0.1.2 | GPL-2.0-or-later |
| `opthedging/` | OptHedging 1.0 | GPL-2.0-or-later |
| `optionpricing/` | OptionPricing 0.1.2 | GPL-2.0-only OR GPL-3.0-only |
| `pa/` | pa 1.2-4 | GPL-2.0-only |
| `parma/` | parma 1.7 | GPL-3.0-or-later |
| `pbo/` | pbo 1.3.5 | MIT |
| `peerperformance/` | PeerPerformance 2.4.0 | GPL-2.0-or-later |
| `performanceanalytics/` | PerformanceAnalytics 2.1.0 | GPL-2.0-or-later |
| `portfoliooptim/` | PortfolioOptim 1.1.1 | GPL-3.0-only |
| `rnd/` | RND 1.2 | GPL-2.0-or-later |
| `rmgarch/` | rmgarch 1.4-2 | GPL-3.0-only |
| `rquantlib/` | RQuantLib 0.4.28 | GPL-2.0-or-later |
| `robustbase/` | robustbase 0.99-7 | GPL-2.0-or-later |
| `rtl/` | RTL 1.3.9 | MIT |
| `rugarch/` | rugarch 1.5-6 | GPL-3.0-only |
| `sde/` | sde 2.0.21 | GPL-2.0-or-later |
| `sharper/` | SharpeR 1.4.0 | LGPL-3.0-or-later |
| `stochfactor/` | stochvol 3.2.9 and factorstochvol 1.1.2 | GPL-2.0-or-later |
| `strand/` | strand 0.2.3 | GPL-3.0-only |
| `svdnf/` | SVDNF 0.1.11 | GPL-3.0-only |
| `timsac/` | timsac 1.3.8-6 | GPL-2.0-or-later |
| `tsdyn/` | tsDyn 11.0.5.2 | GPL-2.0-or-later |
| `tseries/` | tseries 0.10-62 | GPL-2.0-only OR GPL-3.0-only |
| `tserieschaos/` | tseriesChaos 0.1-13.1 | GPL-2.0-only |
| `tvm/` | tvm 0.5.2 | MIT |
| `vrtest/` | vrtest 1.2 | GPL-2.0-only |
| `yieldcurves/` | yieldcurves 0.1.0 | MIT |

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
