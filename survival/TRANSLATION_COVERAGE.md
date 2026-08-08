# Translation coverage

Upstream: `survival 3.8-9`.

## Directly translated / computationally represented

| Upstream family | Fortran implementation |
|---|---|
| `survfitKM`, `src/survfitkm.c` | `kaplan_meier`, `kaplan_meier_counting` |
| `survfitAJ`, `src/survfitaj.c` | `aalen_johansen` transition-probability recursion |
| `coxph.fit`, `src/coxfit6.c` | `coxph_fit`, `coxph_fit_counting` |
| `basehaz`, Cox residual kernels | `cox_baseline`, martingale/Schoenfeld residuals |
| `survreg.fit`, built-in distributions | `survreg_fit`, `survreg_loglik` |
| `survdiff.fit`, `src/survdiff2.c` | `survdiff` |
| concordance kernels | `concordance_right` |
| `src/finegray.c` | `finegray_expand` |
| `survSplit` | `surv_split` |
| survival pseudo-values | `pseudo_survival` |
| `pspline` basis/penalty | `pspline_basis` + supplied `splines-fortran` |

## Architectural differences

The upstream Cox code is optimized around sorted arrays, Cholesky kernels, and
R memory layouts. The Fortran implementation evaluates risk sets explicitly.
This has the same Breslow/Efron partial-likelihood equations but trades speed
for a small, transparent standalone implementation, especially for start/stop
data. It is currently O(n^2) in common Cox operations rather than using all of
upstream's specialized sorted scans.

`survreg_fit` uses the same built-in density families and censoring likelihood
for ordinary right censoring, but obtains score/information by finite
differences rather than porting every derivative branch in `survregc1.c`.
Interval- and left-censored `survreg` likelihoods are not yet exposed.

The Aalen-Johansen routine exposes the transition-product recursion directly
instead of reproducing R's `Surv`/state bookkeeping, robust influence arrays,
and model-object assembly.

## Computational areas intentionally not yet exposed

These remain in `original/survival-master/` and are not mislabeled as already
translated:

- exact conditional Cox (`coxexact`, `agexact`, `clogit` exact method);
- penalized Cox/frailty fitting and all frailty control algorithms;
- Aalen additive regression (`aareg`) and tapering;
- Turnbull interval-censored survival fitting;
- robust cluster influence matrices for `survfit` and Cox fits;
- full multi-state Cox prediction and transition-specific variance machinery;
- ratetable/expected-survival and person-years kernels;
- specialized weighted concordance variants and counting-process concordance;
- `cox.zph` covariance/test assembly beyond exported Schoenfeld residuals;
- interval/left-censored and Student-t `survreg` fitting;
- sparse `Matrix`-class paths, formula parsing, model frames/matrices,
  prediction-object assembly, S3 methods, printing, and plotting.

These omissions are documented because `survival` is much larger than an
optimizer package: its R interface layer and specialized native kernels span
hundreds of routines/tests. The present release is a tested core numerical
translation, not a claim that every R-visible method has a Fortran twin.

## Floating-point time equality

The C sources often compare sorted event times with `==`. To remain warning
clean under `gfortran -Werror`, the Fortran translation uses a scale-aware
8-epsilon comparison. This only affects distinctions at near machine precision.
