# Upstream provenance

Source translated: supplied `mlr-master.zip`.

Upstream metadata in the supplied tree:

- Package: mlr
- Version: 2.19.3
- Title: Machine Learning in R
- License: BSD_2_clause + file LICENSE
- URL: https://mlr.mlr-org.com and the mlr-org/mlr GitHub repository
- NeedsCompilation: yes

The only native computational C source in mlr itself is `src/smote.c`; a copy
is retained under `upstream/smote.c` for provenance. The remaining package is
R code, much of which is framework plumbing or wrappers around third-party
machine-learning packages.

The supplied `survival-fortran-v0.1.0` dependency is included verbatim under
`vendor/` with its own notices and provenance. It is used for Cox,
Kaplan-Meier and concordance calculations.
