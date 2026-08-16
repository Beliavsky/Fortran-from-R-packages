# Porting notes

## Dependencies

The supplied translations are used directly:

- `VGAM-fortran`: `lerch_phi`, regularized gamma functions, gamma quantiles, normal CDF;
- `expint-fortran`: upper incomplete gamma in the Gamma-Lomax CDF;
- `pracma-fortran`: the real Lambert W -1 branch used by Lindley and Muth quantiles.

## Numerical changes

- Maxwell quantiles use a gamma quantile identity instead of `optim()`.
- Gamma-Lomax quantiles use the equivalent gamma-transform/Lomax identity.
- Slashed generalized Rayleigh CDFs use the exact incomplete-gamma antiderivative instead of repeated numerical integration.
- Continuous root-based quantiles use bracketed bisection.
- Discrete quantiles use integer bracketing/bisection instead of linear loops.
- The inverse-Gaussian CDF evaluates its exponentially weighted normal-tail term in log space.

## Corrected upstream edge cases

Several R quantile routines implement `lower.tail=FALSE` after the recycling loop and therefore only replace the final element of a vector result. The Fortran scalar API applies tail semantics consistently. `qgld` also ignores its upstream `lower.tail` argument; this is corrected.

The upstream Kumaraswamy and standard-Omega CDFs return zero rather than one at and above the upper endpoint; the Fortran CDFs use the mathematically correct endpoint values. Discrete quantiles return their actual lower support at probability zero rather than the upstream loop artifact `-1`.

`ppldd` documents real `beta`, but the formula divides by beta; beta=0 is therefore treated as invalid. Negative beta values are evaluated where finite, but quantile inversion requires the resulting CDF to be monotone and bracketable.

For the slashed generalized Rayleigh family, the upstream code rejects beta <= 2 even though the density itself exists for beta > 0; this translation preserves the actual package restriction.

RNG streams use Fortran's intrinsic generator and do not reproduce R's `runif()` stream bit-for-bit.
