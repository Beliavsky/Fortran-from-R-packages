# Provenance and dependency notices

This is a numerical translation of **BurStFin 1.3**, Burns Statistics Ltd.,
maintained upstream by Pat Burns <patrick@burns-stat.com>.

- Official archive: https://cran.r-project.org/src/contrib/BurStFin_1.3.tar.gz
- Archive SHA-256: `1068e489bd67a314a25c354899ad7fbc621a7177644c14ee6bd430ba885e5320`
- Package date: 2022-04-16; CRAN publication: 2022-04-18.
- Upstream website: https://www.burns-stat.com/
- Retrieved for this translation: 2026-10-12.
- The source archive's `R/`, `man/`, `DESCRIPTION`, `NAMESPACE`, `NEWS`, and
  `MD5` files are preserved byte-for-byte under `upstream/`. The original MD5
  listing describes the original archive, not the translated tree.
- `LICENSE` reproduces the Unlimited metadata and embedded public-domain
  dedications. No third-party R implementation has been copied into maintained
  Fortran sources.

## Repository dependency inspection

Before adding helpers, the top-level repository at revision
`cbd437d6d92f5481706fd8ea5a7829aaff7b6e24` was checked for compatible APIs.
`rfortran-core` supplies `r_kinds:dp`, `r_descriptive:r_covariance`,
`r_quantiles:r_quantile_type7`, and
`r_ordering:r_sort_values_in_place`. These are reused through a sibling FPM
path, not copied. The `stats` translation currently refers to an unavailable
`r_matrix_statistics` module and is not a usable dependency in this checkout.
The covariance assembly and single-predictor factor-score regression here are
BurStFin-specific compositions, not a replacement copy of that package.

`rfortran-linalg` supplies its high-level `thin_svd` and `symmetric_eigen`
procedures. Its original manifest is unchanged. It pins Federico Perini's
original `perazz/fortran-lapack` revision
`9982df3b660918761b993eaf0ff2e507547e25b2`; that transitive dependency provides
Fortran LAPACK/BLAS with explicit interfaces. Federico Perini is the author of
fortran-lapack; the upstream LAPACK authors' notices and licenses remain in
that separately fetched package. Shared rfortran packages are credited to
Vivek Rao and contributors and retain their MIT notices.

The user-requested direct-LAPACK route is the FPM-compatible Beliavsky fork,
https://github.com/Beliavsky/fortran-lapack. This package does not need direct
LAPACK calls: it uses the sufficient high-level sibling API and therefore
neither adds that fork nor overrides the older sibling's incompatible pin.
No dependency source or compiler cache belongs in the distribution.

No `link = ["lapack", "blas"]` or system BLAS/LAPACK flags are used. Reference
R may itself use system BLAS; it is an independent test oracle only, never a
runtime dependency of the distributed Fortran package.
