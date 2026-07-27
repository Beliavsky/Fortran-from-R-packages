# Porting notes

## License and attribution

`greeks` 1.5.6 declares `MIT + file LICENSE`. The root `LICENSE` contains the
full MIT text, every translated Fortran unit has an MIT SPDX identifier and
upstream attribution, and the original source is retained unmodified.

## Typed results

R returns named numeric vectors or matrices. Fortran returns `greeks_result`,
`mc_greeks_result`, and `implied_vol_result`. The upstream name `lambda` is
represented by `elasticity` because the descriptive name is clearer in typed
code.

## Callbacks

R accepts arbitrary payoff and jump-generation functions. The Fortran version
uses explicit procedure interfaces:

- `payoff_callback`
- `payoff_derivative_callback`
- `jump_sampler_callback`
- `price_callback`

This provides compile-time checking and avoids dynamic R evaluation.

## Numerical precision

The upstream Rcpp path helpers accept several arguments as C++ `float`. The
Fortran translation uses `real(dp)` throughout, including time increments,
volatility, rates, paths, and quadrature integrals.

## Brownian path construction

The upstream `make_BM.cpp` selects `Range(s, s + paths)`, whose inclusive upper
endpoint requests one more increment than a path column contains. The Fortran
implementation consumes exactly one increment per path and time step.

## American tree

The upstream routine stores all node values discounted to time zero and then
applies a European-tree bias correction. The Fortran code uses the equivalent
standard backward-discount CRR recursion and applies the same bias correction:

```text
corrected American = tree American + exact European - tree European
```

For numerical stability, the default spot step used for American delta is
scaled to the spot price. Gamma keeps the upstream `spot / 50` step.

## Implied volatility

The upstream European routine uses unguarded Halley updates and the general
routine uses unguarded Newton updates. The Fortran routines retain those fast
updates when they remain inside a valid bracket and otherwise fall back to
bisection. This prevents negative volatility and non-finite iteration states.

American inversion starts its lower bracket at `0.003`, consistent with the
upstream general routine's low-volatility workaround for lattice models.

## Monte Carlo

The formulas follow the R implementations and Hudde and Rueschendorf (2023).
The Fortran code additionally reports pathwise standard errors.

The specialized arithmetic-Asian Black-Scholes routine reproduces the upstream
regression-control-variate estimator. It regresses each arithmetic contribution
on the difference between the geometric path contribution and its exact
expectation, then returns the fitted intercept.

The random-number generator is self-contained and deterministic. It does not
reproduce `dqrng` bit for bit, so seeded R and Fortran estimates will differ
path by path while targeting the same expectations.

## Jump diffusion

The upstream jump model multiplies the diffusion path by cumulative exponentiated
jumps without a compensating drift adjustment. That convention is preserved.
The default jump distribution is Student t with three degrees of freedom.

## Scalar versus vector parameters

Several R functions allow one model parameter to be a vector. The Fortran API
is intentionally scalar per invocation. A vectorized study should loop over the
parameter array; this avoids ambiguous output ranks and keeps all interfaces
strongly typed.
