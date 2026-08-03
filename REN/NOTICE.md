# Notices

This package is a modern Fortran translation of the computational portions of
REN 0.1.0, "Regularization Ensemble for Robust Portfolio Optimization," by
Hardik Dixit, Shijia Wang, Bonsoo Koo, Cash Looi, and Hong Wang.

Upstream REN declares `AGPL (>= 3)` in its `DESCRIPTION`. The upstream snapshot
is retained in `upstream/REN-master` for provenance. Plotting, R data-frame
plumbing, and parallel-backend setup are intentionally excluded.

The package vendors the supplied completed corpcor translation, licensed
GPL-3.0-or-later; see `vendor/corpcor/LICENSE` and its notices.

The supplied glmnet translation is GPL-2.0-only and is therefore not linked or
redistributed in this AGPL-3.0-or-later package. The Gaussian LASSO, ridge, and
elastic-net subset needed by REN is implemented locally with coordinate descent
and deterministic cross-validation.

The combined REN/corpcor work is distributed under AGPL-3.0-or-later. No upstream authors
are claimed to endorse this translation.
