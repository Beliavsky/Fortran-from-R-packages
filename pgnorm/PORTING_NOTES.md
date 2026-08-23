# Porting notes

## Sigma/scale inconsistency in upstream

The R source and Rd documentation disagree. The source standardizes as `(y-mean)/sigma`, uses default `sigma=1`, and the RNG multiplies a natural-scale variate by `sigma`. Thus the source implements a scale parameter. However, `dpgnorm` omits the corresponding `1/sigma` Jacobian. The Rd files instead describe `sigma` as the actual standard deviation and give a different rescaling involving the natural standard deviation `sigma_p`.

This port uses a coherent **scale** parameterization matching the source CDF/RNG convention and fixes the density Jacobian. `pgnorm_sd(p,scale)` returns the resulting actual standard deviation.

## Monty Python and Ziggurat lookup tables

The upstream Monty Python and Ziggurat implementations require packaged R datasets containing tuned constants for selected `p<1` values. The numerical distribution does not depend on those algorithms. In this port, the public routines with those names use the exact Nardon-Pianca/gamma construction, while `zigsetup` itself is translated and available. The p-generalized polar and rejecting-polar generators are implemented independently.

## Nardon-Pianca sampler

The source's Nardon-Pianca beta/uniform construction is distributionally equivalent to

`|X|^p / p ~ Gamma(1/p, 1)`

with an independent random sign. The Fortran implementation uses this direct identity, which is simpler and numerically robust.
