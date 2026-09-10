# Notices and provenance

## Upstream BivKLD

This directory is a Fortran translation of the computational code in R package
**BivKLD 0.1.0** (package date 2026-08-22; CRAN publication date 2026-09-02).
The supplied upstream `DESCRIPTION` names these authors:

- Saurav Singla (`aut`, `cre`)
- Rafflesia Chackochan (`aut`)

Upstream describes the kernel estimator as following Chackochan, Sankaran and
Unnikrishnan Nair (2026), *Bivariate Kullback-Leibler divergence*,
Communications in Statistics - Theory and Methods 55(1), 292-312,
DOI 10.1080/03610926.2025.2496687.

The upstream package declares `License: GPL-3`. The complete GNU GPL version 3
text is retained as `LICENSE`. Original computational R sources, package
metadata, citation information, and news are retained under `upstream/` for
attribution and traceability. Those files remain attributable to their
upstream authors.

## Dependency review

Upstream R imports the `ks` package for Gaussian KDE plus `Hscv`/`Hns`
bandwidth selection. Before translation, the `Fortran-from-R-packages`
repository was checked and found to contain a top-level `ks` translation.
That port exposes KDE and `Hns`, but its API mapping records multivariate
`Hscv` as deferred. It also declares `GPL-2.0-only` in its FPM metadata.
Consequently this package does not link, vendor, copy, or embed that sibling.
The limited bivariate Gaussian KDE, normal-scale bandwidth rule, and SCV
criterion needed by BivKLD are independently implemented here.

The local SCV implementation preserves a full positive-definite two-dimensional
bandwidth search and the direct unbinned Gaussian SCV pair-functional
criterion. It deliberately uses a streamlined pilot bandwidth and does not
claim exact parity with every pilot, pre-scaling, binning, or optimizer option
of `ks::Hscv`.

No BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, `rfortran-compat`, or translated
R-package dependency source is included.

## Translation changes

- R matrices/data frames and lists are represented by Fortran arrays and the
  `biv_sample` derived type.
- R exceptions are represented by status codes in kernel routines and IEEE
  quiet NaN in exact functions. Positive infinity is retained for discrete
  support mismatch.
- S3 printing, attributes, group-name handling, roxygen documentation, and
  other R-specific presentation/interface behavior are not translated.
- All maintained Fortran source uses one `dp = real64` kind and free-form
  source.
