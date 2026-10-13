# Upstream provenance

- Package: HDInterval
- Version: 0.2.4
- CRAN publication: 2022-11-17 18:10:21 UTC
- Source archive: https://cran.r-project.org/src/contrib/HDInterval_0.2.4.tar.gz
- Source archive SHA-256: `bb07f0edd660a02ed18e578c2798eb8c2db0e181a5e0c3e23db182d13e9494f6`
- Additional source view: https://github.com/cran/HDInterval/tree/34206d7f2374b8ac6d532d7c5f62e5d2d1f07562
- Upstream license designation: GPL-3, preserved here as GPL-3.0-only.

The numerical source files, `DESCRIPTION` and `NAMESPACE` in `upstream/` are
unmodified upstream files. Relevant package documentation and upstream tests are
also retained for review; only upstream figure images are omitted. The retained
upstream MD5 inventory describes the complete original archive, not the contents
of this selected-source directory. No original R source is required at runtime.

`reference/generate_reference.R` sources all six preserved R files and records
actual upstream results, rather than reimplementing them in R. Reference values
were generated on 2026-10-11 using R 4.5.0. FPM tests also check closed-form normal,
logistic, exponential and Beta(2,1) interval identities and logistic inverse-CDF
values independently of the upstream reference outputs.

The implementation uses `rfortran-core` for shared kind/sorting operations and
new, package-specific scalar searches when the existing `stats` FPM dependency
failed its compatibility build. See `NOTICE.md` for details. Upstream R's
`optimize` and `uniroot` source has not been vendored. Numerical method differences
and limits are explicit in the README and coverage metadata.
