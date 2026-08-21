# Changelog

## 0.2.0

- Translated the specialized BNB/BZINB expectation engine from `src/expt.cpp`.
- Translated the inverse-digamma and coupled M-step optimizer from `src/opt.cpp`.
- Translated the historical-maximum EM control logic from `src/em.cpp`.
- Switched `fit_bnb` and `fit_bzinb` to the source EM algorithms by default.
- Retained v0.1 direct likelihood fits as `fit_bnb_direct` and
  `fit_bzinb_direct`.
- Replaced numerical-Hessian BZINB SEs with the upstream score outer-product
  information matrix.
- Added the exact eight-free-parameter information matrix and a native expanded
  9 x 9 covariance including constrained `p4`.
- Added likelihood trajectories to BNB/BZINB fit result types.
- Added source-compatible `idigamma`.
- Added the upstream all-zero fit shortcut.
- Added `pairwise_bzinb_full` with full parameter records, nonzero proportions,
  and optional pair subsampling.
- Added dedicated EM expectation/score/information and pairwise parity tests.

## 0.1.0

- Initial modern Fortran/FPM translation of `bzinb` 1.0.8.
- Added BP, BZIP-A, BZIP-B, BNB, and BZINB probability and simulation routines.
- Translated BZIP-A/B EM estimators.
- Added constrained direct BNB/BZINB likelihood fitting and covariance estimates.
- Added inverse digamma, latent nondropout weights, weighted Pearson correlation,
  and pairwise BZINB rho fitting.
- Preserved the complete GPL-2 upstream source under `upstream/`.
