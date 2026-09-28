# NOTICE and provenance

## Upstream R package

- Package: `optR`
- Version: `1.2.5`
- Date in supplied DESCRIPTION: `2016-11-29`
- Author: Prakash (PKS Prakash) `<prakash2@uwalumni.com>`
- Upstream declared license: `GPL (>= 2)`
- Upstream description: optimization toolbox for linear systems using Gaussian elimination, LU decomposition, Gauss-Seidel, conjugate-gradient-style iteration, and Cholesky methods.

The complete source archive supplied for this translation is preserved under `upstream/`, including its DESCRIPTION, NAMESPACE, R source, manuals, MD5 file, examples embedded in documentation, authorship, and package metadata.

## Fortran translation

The maintained source in `src/`, tests in `test/`, example in `example/`, FPM metadata, and translation documentation are a modern free-form Fortran translation/derivative work of the upstream computational code. They are distributed under GPL-2.0-or-later to remain compatible with the upstream `GPL (>= 2)` declaration.

No BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, or translated R-package dependency is copied or vendored. This package currently has no external numerical dependency: the algorithms implemented by optR are themselves the subject of the translation.

The translation is not a claim of full R compatibility. R-specific formulas, S3 object construction, printing, plotting, interactive behavior, and data-frame/model-frame semantics are intentionally outside the maintained Fortran API.
