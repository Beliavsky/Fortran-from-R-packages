# Attribution, license and provenance

This package translates computational code from CRAN **ESG 1.3**, published
2023-08-29, by **Jean-Charles Croix, Thierry Moudiki, Frédéric Planchet and
Wassim Youssef**. Upstream maintainer: Wassim Youssef
<Wassim.G.Youssef@gmail.com>. Original title: *A Package for Asset Projection*.
The upstream DESCRIPTION declares **GPL (>= 2)**. This translation and its tests
are distributed under **GPL-2.0-or-later**; the complete GPL version 2 text is
provided in LICENSE. Retained original files keep their original attribution.
No additional upstream copyright notice was found in these R source files.

Official source: https://cran.r-project.org/src/contrib/ESG_1.3.tar.gz
Package page: https://cran.r-project.org/package=ESG
Archive SHA-256:
`fc4d15aac614ca8246ef29c4b1ca602ac969ec730dfce0d2b688717eec2285e4`.

The original DESCRIPTION, NAMESPACE, MD5 manifest, R sources and manual pages
are preserved byte-for-byte under upstream/. The upstream data directory is
omitted: supplied market datasets are not needed for this computational
translation. Thus upstream/MD5 describes the original archive, not the subset
retained here. No dependency package source is embedded under upstream/ or src/.

Fortran translation and test additions: Fortran-from-R-packages contributors,
2026. Source-compatible formulas are explicitly distinguished from numerical
extensions and interface changes in README.md and API-COVERAGE.md.

## Reused sibling dependencies

Before implementation, the target repository tree at commit
`cbd437d6d92f5481706fd8ea5a7829aaff7b6e24` was inspected, including rfortran-core,
rfortran-linalg, MultiRNG, randompack and rngWELL. No ESG top-level translation
existed. ESG reuses rfortran-core's `r_kinds:dp`, `r_distributions:r_pnorm` and the stable
`r_transforms:r_expm1/r_log1p` helpers.
The existing `randompack_rng_type` provides explicit, copyable local streams and
Gaussian draws; its implementation is reused rather than copied or recreated.
MultiRNG's RNG uses compiler-global state. rngWELL has different APIs and
additional upstream source restrictions. Neither is used.

- rfortran-core: MIT, Vivek Rao and contributors.
- randompack: BSD-3-Clause, Kristján Jónasson, University of Iceland; retain and
  consult that sibling's NOTICE.md and THIRD-PARTY-NOTICES for RNG, Ziggurat,
  OpenLibm and other component attribution.
- randompack depends on the unchanged rfortran-linalg sibling (MIT, Vivek Rao
  and contributors). Its own pinned FPM dependency is Federico Perini's
  fortran-lapack, revision `9982df3b660918761b993eaf0ff2e507547e25b2`, which
  compiles Fortran BLAS/LAPACK sources. ESG introduces no direct LAPACK calls,
  copied implementation, or system BLAS/LAPACK link flags. The existing sibling
  currently names the upstream perazz repository; this translation does not
  change that shared dependency. For new direct-LAPACK integrations the
  repository's prescribed FPM-compatible Beliavsky fork is
  https://github.com/Beliavsky/fortran-lapack (original author Federico Perini).

Dependency sources and their notices remain in their respective sibling
packages and are intentionally excluded from an ESG-only delivery archive.
