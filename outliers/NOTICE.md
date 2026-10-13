# Notices and provenance

This work translates the computational R code of **outliers 0.15**, authored and
maintained by **Lukasz Komsta <lukasz.komsta@umlub.pl>**, released 2022-03-26.
The upstream scores manual credits **Alan Richter** for corrections. All applicable
upstream author notices, documentation, and source bytes are preserved under
`upstream/`. Upstream declares **GPL (>= 2)**; this translation and derived tables
are distributed under **GPL-2.0-or-later**, with the complete GNU GPL version 2 in
LICENSE. No warranty is provided. The translation is not represented as an upstream
release or endorsement. New Fortran, tests and generator scripts are modifications
created in 2026 for the Fortran-from-R-packages project.

Official source: https://cran.r-project.org/src/contrib/outliers_0.15.tar.gz

Archive SHA-256: `cc31d7f2faefd2c3a27f8ce78c7e67d3b321dcd6690292fad2468125e5e635fb`.
The upstream archive MD5 manifest verified all 27 listed files. DESCRIPTION,
NAMESPACE, MD5, all 14 R files and all manuals are preserved byte-for-byte.

Destination repository: https://github.com/Beliavsky/Fortran-from-R-packages

## Dependency audit

The destination repository was checked before implementation. Its focused shared
`rfortran-core` package already provides compatible normal, t, F, chi-square,
quantile, MAD, sorting and descriptive-statistics APIs. It is reused exclusively
through `rfortran-core = { path = "../rfortran-core" }`. The tested shared repository
revision is `cbd437d6d92f5481706fd8ea5a7829aaff7b6e24`.
The dependency is MIT licensed, credited to Vivek Rao and contributors, and retains
its own notices in its sibling directory. No dependency source is included here.
No rfortran-compat, BLAS, LAPACK, ARPACK, translated stats package, copied r.f90 or
system numerical libraries are needed. Package code imports the shared dp kind.

## Statistical/table attribution

The following references are retained from the upstream manuals. They identify
upstream statistical provenance; the numeric tables were transcribed from the
preserved outliers 0.15 source, with all original anomalous values intact.

* Dixon, W.J. (1950), Analysis of extreme values, Ann. Math. Stat. 21(4), 488–506.
* Dixon, W.J. (1951), Ratios involving extreme values, Ann. Math. Stat. 22(1), 68–78.
* Rorabacher, D.B. (1991), Statistical Treatment for Rejection of Deviant Values:
  Critical Values of Dixon Q Parameter and Related Subrange Ratios at the 95
  percent Confidence Level. The complete upstream citation is preserved in
  upstream/man/qdixon.Rd without editorial corrections.
* Grubbs, F.E. (1950), Sample Criteria for testing outlying observations,
  Ann. Math. Stat. 21(1), 27–58.
* Pearson, E.S. and Sekar, C.C. (1936), The efficiency of statistical tools and a
  criterion for the rejection of outlying observations, Biometrika 28(3), 308–320.
* David, H.A., Hartley, H.O. and Pearson, E.S. (1954), The distribution of the ratio,
  in a single normal sample, of range to standard deviation, Biometrika 41(3), 482–493.
* Shiffler/Schiffler attribution and original citation spelling are retained in
  upstream/man/scores.Rd.

The Fortran interpolation uses equivalent direct polynomial evaluation instead
of copying R's lm/poly implementation. Grubbs type 11 uses algebraic t inversion
inside the source's clipping and valid bracket rather than copying R's root solver.
