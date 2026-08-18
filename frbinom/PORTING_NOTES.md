# Porting notes

## Scope

This is a complete translation of the computational API of `frbinom` 1.0.0.
The only graphics in the upstream repository occur in README examples and are
not part of the package API.

## Recurrence

The R implementation constructs waiting-time probabilities for the generalized
Bernoulli process and then fills a lower-triangular matrix of repeated renewal
convolutions. The Fortran `count_distribution` routine translates that
recurrence directly.

The same PMF/CDF table is shared by density, CDF, quantile, and random
generation. This avoids the repeated table reconstruction done by separate R
calls.

## Family I

The upstream admissibility condition is retained:

`0 <= c <= min(1-p, (-2p + 2^(2H-2) + sqrt(4p-p*2^(2H)+2^(4H-4)))/2)`.

Although the documentation describes `p` and `H` as open-interval parameters,
the executable R code accepts their endpoints; the Fortran validation follows
the executable code.

For `size=1`, upstream `dfrbinom` explicitly returns an ordinary Bernoulli
distribution with parameter `prob`, independent of `c`, `h`, and `start`.
The Fortran PMF table preserves that special convention.

When `c=0`, the recurrence reduces to the ordinary binomial distribution; this
identity is tested.

## Family II

The executable parameter domain is retained:

- `0.5 <= h <= 1`;
- `0 < c < 2^(2h-2)`;
- `0 < la < c`.

With `start=.false.`, the initial success probability is
`la*size^(2h-2)`. With `start=.true.`, the upstream process starts from the
ordinary family-II waiting law whose first probability is `c`; consequently
`la` has no effect on that branch for `size>1`. The Fortran routines accept
`la` consistently but preserve this model behavior.

For `size=1`, the upstream special case is Bernoulli(`la`) when `start=FALSE`
and Bernoulli(`c`) when `start=TRUE`.

## Quantiles

R searches the precomputed cumulative table with `which`. The Fortran version
performs integer bisection on the same table. Invalid probabilities return the
integer sentinel `-huge(integer)`; probabilities 0 and 1 return support
endpoints directly.

The Fortran API additionally supports `lower_tail=.false.` consistently for
both families.

## RNG

The upstream algorithms are inverse-CDF sampling using `runif`. The Fortran
translation uses the same inversion strategy with the intrinsic random-number
generator. Seeded streams therefore do not reproduce R's RNG bit-for-bit.
