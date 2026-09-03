# Shared translated dependencies

Translated dependencies are kept once at the repository top level when their
APIs are compatible. Consumer packages use sibling FPM path dependencies rather
than carrying independently editable vendored copies. This reduces checkout
size and, more importantly, prevents fixes and formatting changes from drifting
between copies.

The consolidation passes removed 17,752 duplicated tracked files totaling
about 155.39 MiB from the working tree. They established these canonical dependency
paths:

The latest pass also consolidated newly added, previously untracked packages.
The exact-source scan fell from 102 duplicate groups representing about
7.67 MiB of redundant maintained Fortran to two pre-existing `fe-r`/`FER`
test and demo groups totaling about 4.9 KiB.

| Canonical package | Consumers using the shared package |
|---|---|
| `rugarch` | `portvine`, `PWEV`, `quarks` |
| `spacefillr` | `TruncatedNormal` |
| `deSolve` | `rootSolve`, `hypergeo`, `flexsurv` |
| `numDeriv` | `alabama`, `survey`, `compound.Cox`, `gkwdist`, `lavaan`, `flexsurv`, `pbkrtest` |
| `roptim` | `alabama` |
| `alabama` | `TruncatedNormal`, `mbbefd` |
| `nleqslv` | `TruncatedNormal`, `mev` |
| `splines` | `survival`, `gamlss`, `mgcv`, `VGAM` |
| `survival` | `survey`, `gamlss`, `mlr`, `compound.Cox`, `relsurv`, `flexsurv`, `mstate` |
| `minqa` | `survey`, `lme4` |
| `survey` | `GB2` |
| `lpSolve` | `adagio`, `clue`, `limSolve`, `linprog`, `matchingMarkets` |
| `pracma` | `new.dist`, `poweRlaw` |
| `rrcov` | `RobStatTM`, `MASS` |
| `maxLik` | `rumidas` |
| `rumidas` | `PWEV` |
| `Rsolnp` | `DiscreteInverseWeibull`, `DiscreteWeibull` |
| `quadprog` | `limSolve`, `quadprogXT`, `BB`, `flexsurv`, `NlcOptim`, `pracma`, `INFOSET` |
| `Rfast` | `Rfast2` |
| `DiceKriging` | `mlrMBO`, `GPareto` |
| `partitions` | `hyper2`, `MM` |
| `cubature` | `hyper2` |
| `Matrix` | `MatrixExtra`, `piqp` |
| `MatrixExtra` | `ECOSolveR/integration/matrixextra-adapter` |
| `expint` | `actuar`, `mev`, `new.dist` |
| `VGAM` | `new.dist` |
| `nlme` | `gamlss`, `segmented` |
| `mvtnorm` | `ks`, `matrixNormal`, `mc2d`, `mixSPE`, `tmvtnorm` |
| `fitdistrplus` | `mbbefd` |
| `zigg` | `Rfast` |
| `lbfgs` | `RcppNumerical` |
| `matchingR` | `matchingMarkets` |
| `polynom` | `orthopolynom` |
| `coda` | `MCMCpack` |
| `mcmc` | `MCMCpack` |
| `quantreg` | `MCMCpack` |
| `robustbase` | `RobStatTM`, `compositions` |
| `Trading` | `SACCR`, `xVA` |
| `SACCR` | `xVA` |
| `copula` | `ewens` |
| `QCSIS` | `wqc` |
| `waveslim` | `wqc` |
| `elliptic` | `hypergeo` |
| `contfrac` | `hypergeo` |
| `pbivnorm` | `lavaan` |
| `GPArotation` | `lavaan` |
| `pdqutils` | `sadists` |
| `corpcor` | `REN` |
| `rvinecopulib` | `portvine` |
| `COMPoissonReg` | `DiscreteDists` |
| `gamlss.dist` | `gamlss` |
| `actuar` | `mbbefd` |
| `lbfgsb3` | `NFCP`, `RcppNumerical`, `roptim` |
| `ghyp` | `sharpeRratio`, `tsdistributions`, `tsgarch`, `tsmarch` |
| `tsdistributions` | `tsgarch`, `tsmarch` |
| `garchx` | `tvgarch` |
| `tensorA` | `compositions` |
| `NMOF` | `neighbours/integration/nmof-demo` |
| `leaps` | `tsa` |
| `RPEIF` | `RPESE` |
| `deoptimr` | `RSDC` |
| `tvm` | `yrnd` |
| `quadform` | `MM` |
| `DEoptim` | `trawl` |
| `rngWELL` | `randtoolbox` |
| `relsurv` | `flexsurv`, `mstate` |
| `anMC` | `KrigInv` |
| `GA` | `rmoo` |
| `fitHeavyTail` | `highOrderPortfolios` |
| `fastcluster` | `cluster` |
| `bayesm` | `compositions` |
| `ape` | `MCMCglmm` |
| `gRbase` | `gRain` |
| `TMB` | `glmmTMB` |
| `RPEGLMEN` | `RPESE` |
| `fastmatrix` | `L1pack` |
| `AdequacyModel` | `BGFD` |
| `tsgarch` | `tsmarch` |
| `RSpectra` | `bigstatsr` |
| `lsei` | `nspmix` |
| `nnls` | `isotone` |
| `fracdiff` | `forecast`, `ufRisk` |
| `smoots` | `ufRisk` |
| `RobStatTM` | `RPEIF` |
| `nnet` | `forecast` |
| `qrng` | `TruncatedNormal` |
| `tweedie` | `statmod` |
| `urca` | `forecast` |
| `rfortran-compat` | `CompQuadForm`, `DPQ`, `evd`, `gmm`, `matrixdist`, `nnet`, `pearsonds`, `qrng`, `spam`, `SpatialExtremes`, `stabledist`, `statmod`, `TruncatedNormal`, `truncnorm`, `tweedie` |
| `rfortran-optional` | `rfortran-core`, `rfortran-linalg` |
| `rfortran-core` | `ape`, `bayesgarch`, `bayesm`, `changepoint`, `cmprsk`, `corpcor`, `DiscreteWeibull`, `fda`, `FinTS`, `fitdistrplus`, `fportfolio`, `fracdiff`, `GB2`, `geepack`, `gkwdist`, `gRain`, `gRbase`, `isotone`, `MCMCglmm`, `mice`, `mitml`, `mitools`, `pbkrtest`, `performanceanalytics`, `quarks`, `randomForest`, `ranger`, `rrcov`, `rugarch`, `spantest`, `strucchange`, `survey`, `tseries`, `vares`, `vars`, `vrtest`, `waveslim` |
| `rfortran-linalg` | `ape`, `apt`, `bayesianOU`, `BEKKs`, `cccp`, `CEoptim`, `changepoint`, `CLA`, `cmaes`, `cmprsk`, `compositions`, `esback`, `etrm`, `expm`, `fastmatrix`, `fbasics`, `fbonds`, `fcopulae`, `fda`, `fmultivar`, `fnonlinear`, `fportfolio`, `garchx`, `geepack`, `gmm`, `gogarch`, `gRbase`, `irlba`, `ks`, `lgarch`, `lmtest`, `matchingMarkets`, `matrixdist`, `mclust`, `MCMCglmm`, `mice`, `mitml`, `mixsqp`, `msm`, `MultiATSM`, `nmof`, `nnet`, `pa`, `pbkrtest`, `randomForest`, `Rcsdp`, `Rdsdp`, `riskParityPortfolio`, `RiskPortfolios`, `Rmalschains`, `robustbase`, `rquantlib`, `SpatialExtremes`, `statmod`, `stochfactor`, `strucchange`, `tsdyn`, `tvgarch`, `vars` |

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
