# Porting notes

## Scope

This is a complete translation of the computational code in `polyaAeppli`
2.0.2. Plotting occurs only in documentation examples and is not part of the
package source API.

## Probability recurrence

The translated PMF uses the same Johnson-Kotz-Kemp / Nuel recurrence as the R
package, but stores it directly in a zero-based Fortran array indexed by the
random variable value.

Lower-tail log probabilities use stable log-addition. Upper-tail log
probabilities use the dedicated forward-tail recurrence translated from
`logTailPA`, followed by the reverse accumulation translated from `hArray`.
This avoids forming `1-CDF` in extreme upper tails.

## Quantiles

The R implementation obtains an initial upper bound from `qgamma`. The Fortran
translation removes that base-R dependency: it starts from the exact
Polya-Aeppli mean/variance and doubles the bracket until the requested
log-probability is crossed, then performs integer bisection. The resulting
quantile definition remains the smallest integer satisfying the requested
lower- or upper-tail probability criterion.

Probability endpoints are represented by ordinary integer sentinels:
`huge(integer)` denotes the infinite upper quantile.

## Random generation

The upstream algorithm is preserved: draw a Poisson number of clusters and sum
that many shifted geometric random variables. The Fortran implementation uses
a Knuth Poisson generator for small lambda and PTRS transformed rejection for
large lambda, plus inverse-transform geometric generation.

Fortran's intrinsic RNG supplies the uniform stream, so a seed does not
reproduce R's RNG bit-for-bit.

## Vector recycling

The vector entry points reproduce R's `rep_len` parameter recycling. The caller
provides an output array of length equal to the maximum input-array length.
