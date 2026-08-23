# Porting notes

## Source basis

This port follows `ChernoffDist` 0.1.0 by Haitian Xie. The upstream package says
its implementation is based on Groeneboom and Wellner (2001), *Computing
Chernoff's Distribution*, Journal of Computational and Graphical Statistics
10(2), 388-400.

## Airy dependency

The R implementation imports `gsl::airy_zero_Ai` and `gsl::airy_Ai_deriv`.
Only fixed Airy data are required: 20 zeros in one series and 100 zeros plus
derivatives in another. The Fortran translation embeds 100 high-accuracy
binary64 values of `a_k` and `Ai'(a_k)`. This removes an otherwise unnecessary
runtime GSL dependency while preserving the numerical method.

## Effective coefficient set

The R source defines 30 `a` coefficients and 30 `b` entries, but sets `m1=20`
and loops only over `k=1:20`. Several later `b` literals are also visibly
malformed/extreme, but are unreachable. The Fortran implementation therefore
stores only the first 20 coefficients actually used by the algorithm. No active
upstream computation is changed by doing so.

## Integration

R `integrate()` calls are replaced by recursive adaptive 15-point
Gauss-Kronrod integration. The two upstream improper integrals are truncated at
12 after checking their rapidly decaying tails; at the tolerances used here the
omitted tails are negligible relative to binary64 accuracy for the distribution
range of interest.

The CDF integral is split at `x=1`, where the density's internal representation
changes, even though the two representations agree continuously there.

## Quantile

Upstream uses `uniroot()` with a narrow normal-based bracket. The Fortran
routine uses a safeguarded Newton/bisection solver and the translated density as
the derivative. This is more robust in extreme tails while targeting the same
CDF.

For convenience, exact probability endpoints return infinities and the Fortran
API supports `lower_tail` and `log_p`; the R package exposes only the basic
probability argument.

## Removed side effect

`R/dChern.R` ends with a top-level `dChern(0)` call. This causes an unnecessary
density evaluation when the source is loaded and is not part of the exported
API. It is intentionally omitted.
