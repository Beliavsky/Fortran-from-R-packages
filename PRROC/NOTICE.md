# License, attribution and provenance

## Original R package

- Package: PRROC, version 1.4
- Title: Precision-Recall and ROC Curves for Weighted and Unweighted Data
- Authors: Jan Grau and Jens Keilwagen
- Maintainer named upstream: Jan Grau, grau@informatik.uni-halle.de
- Upstream license declaration: GPL-3 (GNU GPL version 3 only)
- Upstream date: 2025-03-15; CRAN publication: 2025-03-18
- Official release archive: https://cran.r-project.org/src/contrib/PRROC_1.4.tar.gz
- Archive SHA-256: `9dbb4a7a9becb5a8f281a05fe263c4c96b187e459f1ed96b27381e838db493b6`
- Package page: https://cran.r-project.org/package=PRROC
- Source obtained and translation validated: 2026-10-11

`upstream/R/PRROC.R`, `upstream/DESCRIPTION`, `upstream/NAMESPACE`, and
`upstream/CITATION` are unmodified copies from that release. They preserve
its source and authorship notices, provide auditable source locations for
coverage metadata, and drive actual-R regression fixture regeneration.
The complete GNU GPL version 3 text is retained in `LICENSE`.

`src/prroc_api.f90` translates the numerical behavior of the two computational
exports and their helpers into pure modern Fortran. Derived-type results,
validation/status handling, stable class scaling/logarithms, bounded curve
allocation and corrected zero-weight/random-baseline edge behavior are new
translation changes. Tests, example and documentation accompany this work.
The translated/added code is distributed under GPL-3.0-only. This package
provides no additional license exception or broader permission for retained
upstream material.

## Scientific attribution

- Davis, J. and Goadrich, M. (2006). The relationship between precision-recall
  and ROC curves. Proceedings of the 23rd International Conference on
  Machine Learning, pp. 233–240. DOI: https://doi.org/10.1145/1143844.1143874
- Keilwagen, J., Grosse, I. and Grau, J. (2014). Area under Precision-Recall
  Curves for Weighted and Unweighted Data. PLOS ONE 9(3).
  DOI: https://doi.org/10.1371/journal.pone.0092209
- Grau, J., Grosse, I. and Keilwagen, J. (2015). PRROC: computing and
  visualizing precision-recall and receiver operating characteristic curves
  in R. Bioinformatics 31(15), 2595–2597.
  DOI: https://doi.org/10.1093/bioinformatics/btv153

The package's original machine-readable R citation instructions, including
scientific authors not listed as package authors, remain in
`upstream/CITATION`.

## Dependency reuse audit

The target Fortran-from-R-packages repository was inspected before
implementation. No existing top-level PRROC translation was found.
Compatible facilities were found in the existing MIT-licensed
`rfortran-core`: the `dp` kind, stable ordering, `r_log1p`, `r_log1pexp`, and
`r_logistic`. `fpm.toml` references this package through the sibling path
`../rfortran-core`. Its authorship is Vivek Rao and contributors, as stated
by its manifest; its separate license and attribution remain with that
package. No dependency source was copied into this directory. No BLAS,
LAPACK, ARPACK, `r.f90`, or `r_mod.f90` is embedded or linked externally.

The R package's `rlang` import serves `.data` in the excluded ggplot interface;
it is not a numerical dependency for the translated operations or source-
loaded R reference generation. R `graphics`, `grDevices`, and suggested
plotting/testing packages are likewise not translated or vendored.
