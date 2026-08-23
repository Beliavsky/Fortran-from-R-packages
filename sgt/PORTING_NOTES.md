# Porting notes

## Upstream source

This translation is based on `sgt` 2.0 by Carter Davis. The supplied upstream
archive is retained unchanged under `upstream/sgt-master.zip`.

## R formula interface

`sgt.mle` accepts arbitrary R formulas for the observation and all five SGT
parameters. Reimplementing an R language evaluator is outside the numerical
scope of this port. `sgt_mle_model` instead accepts a Fortran callback that
returns, for each observation, the transformed `x` and the corresponding
`mu`, `sigma`, `lambda`, `p`, and `q`. This retains the numerical generality of
the regression likelihood while using a natural Fortran interface.

`sgt_mle_constant` provides the common constant-parameter case directly.

## Optimization

The R package calls `optimx` and by default tries Nelder-Mead and BFGS. The
Fortran port contains small native implementations of those two algorithms and
selects the fit with the lower negative log likelihood.

For the constant-parameter convenience routine, positive parameters are
optimized on log scales and `lambda` on an `atanh`/`tanh` scale. This prevents
invalid trial parameters while targeting the same likelihood. Final gradients
and Hessians are recomputed in the original parameter scale.

## Numerical derivatives

The R package uses `numDeriv` for the final gradient and Hessian. This port uses
independently implemented centered finite differences for those calculations.
The supplied dependency snapshots are kept for reference but are not active
build dependencies.

## Dependency licensing

Upstream `sgt` declares `GPL (>= 3)`. The supplied `optimx-fortran` archive is
marked GPL-2.0-only. Its bundled original `optimx` DESCRIPTION also says
`License: GPL-2`. GPL-2.0-only code cannot conservatively be combined into a
GPL-3.0-or-later binary.

The supplied `numDeriv-fortran` archive labels itself GPL-2.0-or-later, while
its bundled original DESCRIPTION says `License: GPL-2`. Because that metadata
is ambiguous, it is also retained as reference-only rather than linked.

The active `sgt-fortran` library is therefore self-contained GPL-3.0-or-later.

## Numerical details

The beta and gamma CDFs and inverse CDFs required by the SGT formulas are
implemented internally using continued fractions/series and bracketed
inversion.

The quantile implementation explicitly handles the finite-SGT split
probability `(1-lambda)/2`. At that probability the incomplete-beta quantile
is zero and the distribution quantile is exactly the mode before optional mean
centering.

The R package contains a warning-condition typo in `dsgt`, `psgt`, and `qsgt`:
the variance warning checks `goodmean` a second time instead of `goodvar`. The
Fortran API returns NaN for invalid mean/variance adjustment conditions and has
no equivalent warning side effect.
