# Provenance, attribution and changes

This is a modified, free-form Fortran translation of computational code from
**loo 2.10.1**, not an official upstream release. Translation date: 2026-10-11.

- Released source: https://cran.r-project.org/src/contrib/loo_2.10.1.tar.gz
- Source archive SHA-256: `a4c7759ec0923f775fd7b5e4cf8341e7f7c3520d6b68d2d5b3af0aafc465b432`
- Release publication: 2026-07-24; package date: 2026-07-23.
- Upstream: https://github.com/stan-dev/loo and https://mc-stan.org/loo/
- Destination: https://github.com/Beliavsky/Fortran-from-R-packages
- Inspected target-repository snapshot: `cbd437d6d92f5481706fd8ea5a7829aaff7b6e24`
- Original source and DESCRIPTION/NAMESPACE/NEWS are retained verbatim in
  `upstream/`; `upstream/SHA256SUMS` records their hashes. No dependency source
  is included. The R files are reference material, never a runtime dependency.

## Authors and licenses

Upstream authors: Aki Vehtari, Jonah Gabry, Måns Magnusson, Yuling Yao,
Paul-Christian Bürkner, Topi Paananen, and Andrew Gelman. Contributors: Ben
Goodrich, Juho Piironen, Bruno Nicenboim, Leevi Lindgren, Visruth Srimath Kandali,
and Florence Bockting. Original authorship and maintainer details remain in
upstream/DESCRIPTION. The algorithms and retained original documentation remain
attributed to their respective upstream authors.

Released DESCRIPTION licenses the package GPL (>= 3). The translated Fortran,
new example and tests are distributed under GPL-3.0-or-later; full GPL v3 text is
in LICENSE. Upstream's repository LICENSE additionally identifies documentation
as CC BY 4.0 and includes both complete texts; that unmodified combined notice is
retained in LICENSE.upstream. Its provenance is upstream tree
`87024b152520c34f9e2373f6655234ee601b1124`, file LICENSE:
https://github.com/stan-dev/loo/blob/87024b152520c34f9e2373f6655234ee601b1124/LICENSE
This development-tree license is retained for notices only; it does not replace
the released numerical-source pin. Adapted documentation has been reorganized
and written for native Fortran arrays, statuses and explicit callbacks.

## Changes from released R

R objects, names, dynamic method dispatch, plotting, printing, parallel process
management and RNG are replaced by explicit numeric arrays, derived result types,
status values and pure computation wherever possible. Four full model-refitting
or subsampling exports remain untranslated, as listed in API_COVERAGE.md.
Callback APIs materialize the complete likelihood matrix; they are not streaming
memory-saving replicas of R's observation-by-observation parallel implementation.
Caller-supplied permutations and Bayesian-bootstrap weights make stochastic
operations reproducible without promising R RNG-sequence compatibility.
Stacking retains the original concave log-mixture score but uses a closed-simplex
vertex-exchange solver with exact one-dimensional line searches and an explicit
duality-gap convergence certificate rather than R stats::constrOptim.

Numerical improvements retain well-conditioned release outputs while correcting
release failure modes. Stable exponential tail differences, centered weighted
variance, stricter finite likelihood validation and log-domain MCSE are informed
by the development source at the same tree pin above (R/psis.R, R/helpers.R,
R/E_loo.R and R/loo.R). The Fortran implementation computes log-domain MCSE directly
rather than only as a fallback. Centered normalized weights are stored independently
of rounded raw-scale metadata, and tail selection branches before division by tiny
relative efficiency. These changes are tested by mathematical invariance identities,
not by treating release NaN or negative-variance bugs as reference answers.

## Existing dependency reuse

The target repository was inspected before implementation. Compatible sibling FPM
path dependencies are used:

- `posterior` 1.7.0 (BSD-3-Clause): `gpdfit`, `qgeneralized_pareto`, `pareto_khat`
  and `ess_mean`. The existing port was reused unchanged. Loo-specific ceiling
  tail lengths are retained locally rather than substituting posterior's floor.
- `matrixStats` 1.5.0 (Artistic-2.0): stable `log_sum_exp` and `col_vars`.
- `rfortran-core` (MIT): the shared `dp` kind, ordering, type-seven quantiles,
  stable logarithm/exponential transforms and descriptive helpers.

No copies of these packages, r.f90, r_mod.f90, BLAS, LAPACK, ARPACK or compatibility
modules are embedded. No system numerical library is linked, and neither direct
LAPACK nor rfortran-linalg is needed by these array/importance computations.

Algorithms and references include Vehtari, Gelman and Gabry (2017), Practical
Bayesian model evaluation using leave-one-out cross-validation and WAIC,
doi:10.1007/s11222-016-9696-4; Vehtari, Simpson, Gelman, Yao and Gabry (2024),
Pareto Smoothed Importance Sampling; Zhang and Stephens (2009), A New and Efficient
Estimation Method for the Generalized Pareto Distribution; Ionides (2008),
Truncated Importance Sampling; and the model-stacking/scoring references retained
in the corresponding upstream R source documentation.
