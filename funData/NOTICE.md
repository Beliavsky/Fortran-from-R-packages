# NOTICE and provenance

## Upstream

This project is a modern Fortran translation of the computational portions of the R package **funData**, version 1.3-9 (2024-02-10).

- Upstream author and maintainer: Clara Happ-Kurz
- Upstream project: https://github.com/ClaraHapp/funData
- Upstream license: GPL-2
- Upstream citation: Happ-Kurz, C. (2020), "Object-Oriented Software for Functional Data", Journal of Statistical Software 93(5), 1-38, doi:10.18637/jss.v093.i05.

The original DESCRIPTION, NAMESPACE, R sources, and CITATION used for this translation are retained under `upstream/` for provenance. The GPL-2 license text is retained as `LICENSE`.

## Translation notes

The R package uses S4 classes and imports `abind`, `fields`, `foreach`, graphics packages, and `stats`. This Fortran translation implements its mapped numerical operations directly and does not copy or vendor those dependencies. It also does not copy BLAS, LAPACK, `r.f90`, `r_mod.f90`, or any translated R-package dependency.

The target `Fortran-from-R-packages` repository was reviewed before implementation. No sibling numerical package was needed for the translated surface, so this package keeps one package-local `dp = real64` kind and has no FPM dependencies.

For one-point irregular full-domain extrapolation, upstream documentation says the observed function is extended as a constant. The upstream `.extrapolate` implementation assigns `rep(x, 3)`, which appears to use the observation-grid value rather than the function value. This translation follows the documented constant-function-value behavior.

Random simulation, sparsification, and noise routines use Fortran `random_number` plus Box-Muller normal generation. Therefore seeded streams are not expected to reproduce R's RNG sequence exactly.
