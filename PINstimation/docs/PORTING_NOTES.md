# Porting notes

## Data orientation

`trade_counts%buys(i)` and `trade_counts%sells(i)` hold the two daily count
series. MPIN posterior columns are ordered as

```text
no information, good layer 1, bad layer 1, ..., good layer J, bad layer J
```

AdjPIN cluster order follows the upstream package:

1. no information, no shock
2. no information, shock
3. good information, no shock
4. good information, shock
5. bad information, no shock
6. bad information, shock

## Constraints

Probability parameters are optimized on a logit scale. Positive rates use a
softplus transformation. MPIN layer probabilities use a softmax with an
explicit no-information baseline, guaranteeing `sum(alpha) < 1`. MPIN rates
are parameterized as positive cumulative increments, giving increasing `mu`
without a discontinuous sorting operation during optimization.

AdjPIN restrictions are represented by `adjpin_restrictions` and tie buy/sell
parameters before likelihood evaluation.

## Optimizers

The project contains self-contained Nelder-Mead and numerical-gradient BFGS
implementations. Nelder-Mead is the default because it mirrors the main
upstream estimation path and is robust to flat mixture-likelihood regions.

## ECM adaptations

The MPIN ECM computes exact mixture posteriors and exact probability updates.
The shared baseline/informed Poisson rates do not have a convenient simultaneous
closed form, so the conditional rate M-step is solved numerically.

The AdjPIN ECM likewise uses exact six-state posterior probabilities and a
numerical complete-data maximization. This is computationally faithful but does
not reproduce every specialized rational-equation shortcut in the R source.

## PIN initialization wrappers

The upstream EA, GWJ, and YZ routines contain extensive data-dependent grids and
boundary-handling logic. The Fortran wrappers preserve their roles as distinct
initializer families but use compact moment/quantile starts before the common
stable MLE engine. They are therefore adapted interfaces, not line-for-line
ports.

## VPIN input boundary

The R function accepts raw timestamp/price/volume records and performs
calendar-aware bar completion. The Fortran `compute_vpin()` interface begins at
the numerical time-bar stage: price change, volume, and duration arrays are
supplied by the caller. It then implements bulk classification, proportional
bar splitting, fixed-volume buckets, rolling VPIN, and optional IVPIN.

`classify_trades()` accepts timestamps as numeric seconds. Its `timelag` is also
in seconds, whereas the R-facing argument is documented in smaller time units.

## Numerical constants

Factorials are evaluated with `log_gamma`, mixture likelihoods use log-sum-exp,
and zero mixture weights map to negative infinity. Random Poisson generation
uses inversion for small rates and transformed rejection for large rates.
