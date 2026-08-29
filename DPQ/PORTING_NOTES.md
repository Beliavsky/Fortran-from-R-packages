# Porting notes

## Translation policy

Computational/statistical code was translated to standard free-form modern
Fortran. Plotting, R object-system glue, Rmpfr dispatch, formatting, and
interactive diagnostics were omitted. The supplied MIT `r_mod.f90` is reused
for standard distribution functions, special functions, and other applicable
R-compatible helpers instead of duplicating them.

## Numerical organization

The port is split into:

- `dpq_core`: DPQ tail/log macros expressed as functions; log-space arithmetic;
  special series; Chebyshev and low-level helpers.
- `dpq_gamma_discrete`: `bd0`, Stirling corrections, Gamma/Poisson/binomial/NB
  kernels, q* replications and chi-square/Gamma starting approximations.
- `dpq_normal_beta`: Normal, beta and noncentral-beta approximations.
- `dpq_hyper`: hypergeometric approximations and exact/reference recurrences.
- `dpq_nchisq`: noncentral chi-square exact-mixture calculations and named
  approximations.
- `dpq_t`: t/noncentral-t calculations and named approximations.
- `dpq_wiener`: Wiener-germ noncentral chi-square approximations.
- `dpq_toms1006`: the separately licensed TOMS-1006-facing functionality.

## Deliberate equivalence/consolidation

DPQ is partly a laboratory for comparing historical implementations. The
Fortran port preserves public numerical functionality, but does not falsely
claim that every didactic variant follows a distinct line-for-line code path:

- `pnchisq_rc`, `pnchisq_it`, `pnchisq_v`, and `pnchisq_ss` use the stable
  centered Poisson-mixture evaluator as their common high-accuracy reference
  kernel. Named closed-form approximations remain separate.
- Several historical noncentral-t comparison entry points (`pnt3150`,
  `pntP94`, `pntChShP94`, `pntVW13`, `pntGST23_*`) share the robust Guenther/R
  series where their R-specific research scaffolding did not translate cleanly.
  `pntJW39`, `pntLrg`, JKB density, and `dtWV` retain distinct formulas.
- `lbeta_asy`, `lbeta_m`, and `lbeta_mm` expose the corresponding numerical
  endpoints but use the accurate log-beta helper where this is preferable to
  preserving a deliberately lower-accuracy comparison path.
- `bpser` uses `r_mod::pbeta` for the regularized incomplete beta result rather
  than reproducing R Mathlib's entire TOMS-708 internal call graph.
- `pnbeta_as310` uses a centered Poisson mixture of beta CDFs. The original
  AS-310-derived C source is retained under `upstream/src/310-pnbeta.c`.
- `ebd0` returns the accurate scalar deviance value. R's compensated `(yh,yl)`
  split is not represented as a public pair because the downstream Fortran
  density kernels do not require that representation.

These choices preserve useful DPQ numerical semantics while keeping the port
maintainable and explicit about where algorithm-comparison internals differ.

## TOMS 1006

`lgamma_p11` retains the Pugh coefficient formula from the GPL-3 upstream
component. For the generalized incomplete-Gamma integral, the Fortran API
computes the same mathematical quantity using `r_mod::pgamma` when `mu > 0`
and an exact finite recurrence for the supported negative-`mu`, integer-`p`
case. It is therefore API/mathematical parity, not a line-by-line port of the
upstream continued-fraction/Romberg implementation. The original GPL-3 source
is retained verbatim for provenance.

## Hypergeometric detail

`pdhyper` was ported as its actual upstream recurrence: it returns
`phyper/dhyper`, not a CDF. `phyper_r2` combines this ratio with tail swapping,
matching the purpose of R's `phyperR2` implementation.

## Accuracy and edge cases

The noncentral chi-square CDF/density use centered Poisson summation to avoid
underflow from starting a noncentral mixture at term zero. Quantiles use
bracketed bisection against the high-accuracy CDF. The noncentral-t quantile is
similarly bracketed against the Guenther/R series.

The project uses IEEE NaN/infinities for invalid/boundary cases where natural
Fortran equivalents exist. R warnings/messages are generally represented by
return values/status arguments rather than console output.

## `r_mod.f90`

`upstream/r_mod-original.f90` is the exact supplied helper file. `src/r_mod.F90`
contains formatting-only continuation/line wrapping needed to stay within
standard free-form source line limits. No helper algorithm was intentionally
changed for DPQ.
