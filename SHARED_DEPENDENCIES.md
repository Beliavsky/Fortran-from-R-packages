# Shared translated dependencies

Translated dependencies are kept once at the repository top level when their
APIs are compatible. Consumer packages use sibling FPM path dependencies rather
than carrying independently editable vendored copies. This reduces checkout
size and, more importantly, prevents fixes and formatting changes from drifting
between copies.

Earlier consolidation passes removed 17,752 duplicated tracked files totaling
about 155.39 MiB from the working tree.

Those passes also consolidated newly added, previously untracked packages.
The exact-source scan fell from 102 duplicate groups representing about
7.67 MiB of redundant maintained Fortran to two pre-existing `fe-r`/`FER`
test and demo groups totaling about 4.9 KiB.

The table below reflects direct sibling dependencies in the current top-level
FPM manifests, refreshed on 2026-10-14, plus the separately listed integration
adapters. It excludes development-only dependencies and does not imply that
every consumer has been rebuilt or tested in this documentation refresh.
The duplicate-file figures above describe earlier audits, not a fresh scan.

The newly added `Hmisc` and `RLRsim` translations were checked on 2026-09-20.
Neither declares FPM dependencies, so neither adds a shared-consumer entry
below. Both currently use package-local numerical helpers; their addition is
not a new consolidation or validation pass.

The 2026-09-21 additions include `fda.usc`, which uses `rfortran-core` and
`rfortran-linalg`, and `kSamples`, which uses `rfortran-core` and `SuppDists`.
The added `abind`, `ecp`, `fdapace`, `funData`, `geometry`, and `pcaPP`
packages declare no external FPM dependencies. These are manifest checks,
not new build or numerical-validation results.

The latest additions reuse shared modules in `ftsa`, `matrixcalc`, `mbend`,
and `psych`; `pdfCluster` uses `geometry`, and `rARPACK` uses
`rfortran-arpack`. The `geigen` and `QZ` packages declare direct Git
dependencies on the [Beliavsky fork of fortran-lapack](https://github.com/Beliavsky/fortran-lapack),
originally developed by Federico Perini. Those Git dependencies are not sibling
paths and therefore are not listed in the table. These observations describe
manifests, not a fresh build or regression-test pass.

The October additions also reuse `matrixStats`, `posterior`, `exactRankTests`,
`nortest`, `moments`, `Hmisc`, `TTR`, `randompack`, `expm`, `lpSolve`,
`lbfgsb3`, `splines`, `kendall`, and `randtests` as sibling dependencies.
`matrixTests` additionally declares a pinned Git dependency on the same
`fortran-lapack` fork; this is distinct from the sibling paths listed below.

| Canonical package | Consumers using the shared package |
|---|---|
| `actuar` | `mbbefd` |
| `AdequacyModel` | `BGFD` |
| `alabama` | `mbbefd`, `TruncatedNormal` |
| `anMC` | `KrigInv` |
| `ape` | `MCMCglmm` |
| `bayesm` | `compositions` |
| `coda` | `MCMCpack` |
| `COMPoissonReg` | `DiscreteDists` |
| `contfrac` | `hypergeo` |
| `copula` | `ewens` |
| `corpcor` | `REN` |
| `cubature` | `hyper2` |
| `DEoptim` | `trawl` |
| `deoptimr` | `RSDC` |
| `deSolve` | `flexsurv`, `hypergeo`, `rootSolve` |
| `DiceKriging` | `GPareto`, `mlrMBO` |
| `dplyr` | `tidyr` |
| `elliptic` | `hypergeo` |
| `exactRankTests` | `matrixTests`, `trend` |
| `expint` | `actuar`, `mev`, `new.dist` |
| `expm` | `NetworkRiskMeasures` |
| `extraDistr` | `PropCIs` |
| `fastcluster` | `cluster`, `stats` |
| `fastmatrix` | `L1pack` |
| `fitdistrplus` | `mbbefd` |
| `fitHeavyTail` | `highOrderPortfolios` |
| `forcats` | `readr` |
| `forecast` | `imputeTS` |
| `fracdiff` | `forecast`, `ufRisk` |
| `GA` | `rmoo` |
| `gamlss.dist` | `gamlss` |
| `garchx` | `tvgarch` |
| `geometry` | `pdfCluster` |
| `ghyp` | `sharpeRratio`, `tsdistributions`, `tsgarch`, `tsmarch` |
| `GPArotation` | `lavaan` |
| `gRbase` | `gRain` |
| `Hmisc` | `monobin` |
| `igraph` | `gRbase` |
| `kendall` | `trend` |
| `KFAS` | `MARSS` |
| `lbfgs` | `RcppNumerical` |
| `lbfgsb3` | `binom`, `CreditRisk`, `NFCP`, `RcppNumerical`, `roptim` |
| `leaps` | `tsa` |
| `lme4` | `gamm4` |
| `lpSolve` | `adagio`, `clue`, `irr`, `limSolve`, `linprog`, `matchingMarkets` |
| `lsei` | `nspmix` |
| `matchingR` | `matchingMarkets` |
| `Matrix` | `MatrixExtra`, `piqp` |
| `MatrixExtra` | `ECOSolveR/integration/matrixextra-adapter` |
| `matrixStats` | `loo`, `matrixTests` |
| `maxLik` | `dccmidas`, `rumidas` |
| `mclust` | `otrimle` |
| `mcmc` | `MCMCpack` |
| `mgcv` | `gamm4` |
| `minqa` | `lme4`, `survey` |
| `moments` | `matrixTests` |
| `MTS` | `SteadyStateBVAR` |
| `mvtnorm` | `ks`, `matrixNormal`, `mc2d`, `mixSPE`, `tmvtnorm` |
| `nleqslv` | `mev`, `TruncatedNormal` |
| `nlme` | `gamlss`, `segmented` |
| `NMOF` | `neighbours/integration/nmof-demo` |
| `nnet` | `forecast` |
| `nnls` | `isotone` |
| `nortest` | `matrixTests` |
| `numDeriv` | `alabama`, `compound.Cox`, `flexsurv`, `gkwdist`, `lavaan`, `pbkrtest`, `survey` |
| `optimx` | `dlm` |
| `partitions` | `hyper2`, `MM` |
| `pbivnorm` | `lavaan` |
| `pdqutils` | `sadists` |
| `polynom` | `orthopolynom` |
| `posterior` | `loo` |
| `pracma` | `new.dist`, `poweRlaw` |
| `QCSIS` | `wqc` |
| `qrng` | `TruncatedNormal` |
| `quadform` | `MM` |
| `quadprog` | `BB`, `flexsurv`, `INFOSET`, `limSolve`, `NlcOptim`, `pracma`, `quadprogXT` |
| `quantreg` | `MCMCpack` |
| `randompack` | `ESG`, `varmapack` |
| `randtests` | `trend` |
| `RcppNumerical` | `PropCIs` |
| `readr` | `tidyr` |
| `relsurv` | `flexsurv`, `mstate` |
| `Rfast` | `Rfast2` |
| `rfortran-arpack` | `bigstatsr`, `rARPACK`, `RSpectra` |
| `rfortran-compat` | `CompQuadForm`, `DPQ`, `evd`, `gmm`, `matrixdist`, `nnet`, `pearsonds`, `qrng`, `spam`, `SpatialExtremes`, `stabledist`, `statmod`, `TruncatedNormal`, `truncnorm`, `tweedie` |
| `rfortran-core` | `ape`, `bayesgarch`, `bayesm`, `binom`, `BurStFin`, `changepoint`, `cmprsk`, `corpcor`, `CreditRisk`, `credule`, `dccmidas`, `diptest`, `DiscreteWeibull`, `dlm`, `entropy`, `ESG`, `exactRankTests`, `fda`, `fda.usc`, `FinTS`, `fitdistrplus`, `fportfolio`, `fracdiff`, `fTrading`, `ftsa`, `GB2`, `geepack`, `gkwdist`, `gRain`, `gRbase`, `HDInterval`, `ineq`, `infotheo`, `irr`, `isotone`, `kde1d`, `kendall`, `kSamples`, `loo`, `matrixcalc`, `matrixStats`, `matrixTests`, `MCMCglmm`, `mice`, `mitml`, `mitools`, `monobin`, `NetworkRiskMeasures`, `nortest`, `optionstrat`, `outliers`, `pbkrtest`, `performanceanalytics`, `perm`, `posterior`, `PropCIs`, `PRROC`, `psych`, `quarks`, `randomForest`, `randtests`, `ranger`, `roll`, `rrcov`, `rugarch`, `spantest`, `splines2`, `stats`, `SteadyStateBVAR`, `strucchange`, `survey`, `trend`, `tseries`, `TTR`, `vares`, `vars`, `vrtest`, `waveslim`, `wavethresh` |
| `rfortran-linalg` | `ape`, `apt`, `bayesianOU`, `BEKKs`, `binom`, `BurStFin`, `cccp`, `CEoptim`, `changepoint`, `CLA`, `cmaes`, `cmprsk`, `compositions`, `dccmidas`, `dlm`, `esback`, `etrm`, `expm`, `fastmatrix`, `fbasics`, `fbonds`, `fcopulae`, `fda`, `fda.usc`, `fmultivar`, `fnonlinear`, `forecast`, `fportfolio`, `ftsa`, `gamm4`, `garchx`, `geepack`, `gmm`, `gogarch`, `gRbase`, `irlba`, `irr`, `ks`, `lgarch`, `lmtest`, `MARSS`, `matchingMarkets`, `matrixcalc`, `matrixdist`, `mbend`, `mclust`, `MCMCglmm`, `mice`, `mitml`, `mixsqp`, `msm`, `MultiATSM`, `nmof`, `nnet`, `pa`, `pbkrtest`, `psych`, `randomForest`, `randompack`, `Rcsdp`, `Rdsdp`, `riskParityPortfolio`, `RiskPortfolios`, `Rmalschains`, `robustbase`, `roll`, `rquantlib`, `Rssa`, `SpatialExtremes`, `statmod`, `stats`, `SteadyStateBVAR`, `stochfactor`, `strucchange`, `svd`, `tsdyn`, `tvgarch`, `urca`, `varmapack`, `vars`, `wavethresh` |
| `rfortran-optional` | `rfortran-core`, `rfortran-linalg`, `stats` |
| `rngWELL` | `randtoolbox` |
| `RobStatTM` | `RPEIF` |
| `robustbase` | `compositions`, `RobStatTM` |
| `roll` | `dccmidas` |
| `roptim` | `alabama` |
| `RPEGLMEN` | `RPESE` |
| `RPEIF` | `RPESE` |
| `rrcov` | `MASS`, `RobStatTM` |
| `Rsolnp` | `DiscreteInverseWeibull`, `DiscreteWeibull` |
| `RSpectra` | `bigstatsr`, `svd` |
| `rugarch` | `dccmidas`, `portvine`, `PWEV`, `quarks` |
| `rumidas` | `dccmidas`, `PWEV` |
| `rvinecopulib` | `portvine` |
| `SACCR` | `xVA` |
| `smoots` | `ufRisk` |
| `spacefillr` | `TruncatedNormal` |
| `splines` | `gamlss`, `mgcv`, `splines2`, `survival`, `VGAM` |
| `stinepack` | `imputeTS` |
| `stringr` | `tidyr` |
| `SuppDists` | `kSamples` |
| `survey` | `GB2` |
| `survival` | `compound.Cox`, `flexsurv`, `gamlss`, `mlr`, `mstate`, `relsurv`, `survey` |
| `svd` | `Rssa` |
| `tensorA` | `compositions` |
| `tibble` | `dplyr`, `readr`, `tidyr` |
| `tidyselect` | `dplyr`, `tidyr` |
| `TMB` | `glmmTMB` |
| `Trading` | `SACCR`, `xVA` |
| `tsdistributions` | `tsgarch`, `tsmarch` |
| `tsgarch` | `tsmarch` |
| `TTR` | `fTrading` |
| `tvm` | `yrnd` |
| `tweedie` | `statmod` |
| `urca` | `forecast` |
| `vctrs` | `dplyr`, `readr`, `tibble`, `tidyr` |
| `VGAM` | `new.dist` |
| `waveslim` | `wavethresh`, `wqc` |
| `zigg` | `Rfast` |

For the earlier consolidation passes, each canonical package and affected
consumer passed its FPM test suite before the redundant tree was removed.
Canonical packages retain the applicable licenses, notices, provenance, and
upstream reference material.
Compatibility exports were added to canonical `gamlss.dist` and `actuar` so
their consumers no longer require independently maintained API variants.
The transitional `rfortran-compat` package similarly supplies one renamed
`r_compat` module for older generated translations. Renaming the compatibility
module allows packages such as `forecast` to use old-runtime consumers and the
newer split `rfortran-core` dependency in the same build without duplicate
module names. Fourteen embedded runtime files of roughly 648 KB each were
removed. `SpatialExtremes` additionally uses `rfortran-linalg` for its
Cholesky, SPD solve, inverse, and log-determinant operations.
The shared-runtime consumers passed their tests except for `spam`, whose
bundled legacy ARPACK code still requires system BLAS, LAPACK, and ARPACK
linkage on Windows. `urca` now uses `rfortran-linalg` and its pinned
pure-Fortran LAPACK dependency, removing the system BLAS/LAPACK linker
requirement that also blocked consumers such as `forecast`. `lmtest` now uses
the same shared dependency for least squares, symmetric eigenvalues, dense
solves, and SPD inversion instead of linking system BLAS/LAPACK.
The four `TruncatedNormal` targets pass individually; an all-target FPM build
can encounter a Windows executable-file race.
The `classInt`, `deldir`, `e1071`, `fields`, and `multcomp` translations
currently retain package-local snapshots of `e1071`, `polyclip`, `proxy`,
`spam`, and `mvtnorm`, respectively. Those pinned copies are not counted as
shared consumers above; replacing them with the canonical top-level
translations requires API and regression validation in a later consolidation
pass.
Unreferenced private dependency trees were also removed from `PSDistr`,
`compositions`, `NlcOptim`, `Directional`, `tsa`, `survey`, `mstate`, and
`trawl` after their suites passed without them. The reusable ECOS-MatrixExtra
adapter remains in place and now uses the canonical top-level `MatrixExtra`
package.

Use the repository's package downloader for an individual package with shared
dependencies:

```bat
python download_build_package.py PACKAGE_NAME
```

It recursively adds sibling FPM path dependencies to the sparse checkout.

## Finding remaining exact copies

Run:

```bat
python find_duplicate_fortran_sources.py --include-untracked --minimum-bytes 1000 --top 30
```

The scanner uses tracked working-tree content, optionally includes untracked
files, and excludes original, upstream, reference, build, and Git directories.
Its report is an inventory, not an
instruction to delete files mechanically. Copies embedded directly in package
source trees, intentionally pinned forks, and merely similar implementations
need API and numerical-equivalence review before consolidation.
