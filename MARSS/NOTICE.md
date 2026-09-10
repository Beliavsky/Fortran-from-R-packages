# NOTICE and provenance

## Upstream package

This work is a Fortran translation of computational portions of the R package **MARSS: Multivariate Autoregressive State-Space Modeling**, upstream version **3.11.10** (dated 2025-09-10).

Upstream authors:

- Elizabeth Eli Holmes
- Eric J. Ward
- Mark D. Scheuerell
- Kellie Wills

The upstream package declares the license `GPL-2` in `DESCRIPTION`. The complete GNU General Public License version 2 text is provided in `LICENSE`. This Fortran translation is a derivative work under the same license.

The source snapshot used for this translation was supplied as `MARSS-master.zip`. For traceability, the package metadata (`DESCRIPTION`), export declarations (`NAMESPACE`), upstream citation file, and R sources used in the coverage inventory are retained under `upstream/`. These copies are provenance material from MARSS itself, not vendored third-party dependencies.

## Dependency provenance

The Fortran package declares a sibling FPM path dependency on `rfortran-linalg`. That dependency is maintained separately in the `Fortran-from-R-packages` repository and is **not** copied into this package. The translation also does not copy or vendor KFAS, mvtnorm, nlme, BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, or `rfortran-compat`.

The R package's use of KFAS, mvtnorm, and nlme informed the upstream behavior inventory, but the supported Fortran numerical core uses its own MARSS filter/smoother, deterministic RNG, and finite-difference inference routines together with the shared `rfortran-linalg` API.

## Translation scope

The translation is not an R runtime emulation. R S3 classes and method dispatch, formula/list parsing, names/dimnames, plotting, printed summaries, progress bars, data-frame construction, R warnings, and R RNG streams are outside the maintained Fortran API. MARSS fixed/free/equality constraints are represented numerically through affine f+D*beta blocks rather than R list matrices or named R objects.
