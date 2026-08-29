# Shared translated dependencies

Translated dependencies are kept once at the repository top level when their
APIs are compatible. Consumer packages use sibling FPM path dependencies rather
than carrying independently editable vendored copies. This reduces checkout
size and, more importantly, prevents fixes and formatting changes from drifting
between copies.

The consolidation passes removed 17,752 duplicated tracked files totaling
about 155.39 MiB from the working tree. They established these canonical dependency
paths:

| Canonical package | Consumers using the shared package |
|---|---|
| `rugarch` | `portvine`, `PWEV`, `quarks` |
| `spacefillr` | `TruncatedNormal-fortran-v0.1.0` |
| `deSolve` | `rootSolve`, `hypergeo`, `flexsurv` |
| `numDeriv` | `alabama`, `survey`, `compound.Cox`, `gkwdist`, `lavaan-fortran-v0.7.0`, `flexsurv` |
| `roptim` | `alabama` |
| `alabama` | `TruncatedNormal-fortran-v0.1.0`, `mbbefd` |
| `nleqslv` | `TruncatedNormal-fortran-v0.1.0`, `mev` |
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
| `mvtnorm` | `ks`, `mc2d`, `mixSPE` |
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
| `pbivnorm` | `lavaan-fortran-v0.7.0` |
| `GPArotation` | `lavaan-fortran-v0.7.0` |
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
| `RPEGLMEN` | `RPESE` |
| `fastmatrix` | `L1pack` |
| `AdequacyModel` | `BGFD` |
| `tsgarch` | `tsmarch` |
| `RSpectra` | `bigstatsr` |
| `lsei` | `nspmix` |
| `nnls` | `isotone` |
| `fracdiff` | `ufRisk` |
| `smoots` | `ufRisk` |
| `RobStatTM` | `RPEIF` |

Each canonical package and every affected consumer passed its FPM test suite
before the redundant tree was removed. Canonical packages retain the applicable
licenses, notices, provenance, and upstream reference material.
Compatibility exports were added to canonical `gamlss.dist` and `actuar` so
their consumers no longer require independently maintained API variants.
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
python find_duplicate_fortran_sources.py --minimum-bytes 1000 --top 30
```

The scanner uses tracked working-tree content and excludes `original`, `orig`,
`reference`, build, and Git directories. Its report is an inventory, not an
instruction to delete files mechanically. Copies embedded directly in package
source trees, intentionally pinned forks, and merely similar implementations
need API and numerical-equivalence review before consolidation.
