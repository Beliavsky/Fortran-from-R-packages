# stabledist-fortran

Modern Fortran/FPM translation of the computational core of the R package
`stabledist` 0.7-2.

The port implements stable-distribution density, CDF, quantiles, random
variates, Nolan S0/S1/S2 parameterizations, stable modes, and the package's
tail helpers. Plotting, R object-system code, environment/options handling,
and test/UI scaffolding are not part of the Fortran API.

## Main API

```fortran
use stabledist

f = dstable(x, alpha, beta, gamma=1.0_dp, delta=0.0_dp, pm=0)
p = pstable(x, alpha, beta, gamma=1.0_dp, delta=0.0_dp, pm=0)
q = qstable(p, alpha, beta, gamma=1.0_dp, delta=0.0_dp, pm=0)
x = rstable(n, alpha, beta, gamma=1.0_dp, delta=0.0_dp, pm=0)
m = stable_mode(alpha, beta)
```

`pm=0`, `pm=1`, and `pm=2` correspond to the upstream Nolan S0, S1, and S2
parameterizations. S2 is mode-centered by construction.

Also exported are `c_stable_tail`, `dpareto`, `ppareto`, `stable_omega`, and
`rstable_varying`. Vector `dstable`, `pstable`, and `qstable` overloads support
R-style recycling of vector `gamma` and `delta` arguments.

## Algorithms

* Nolan (1997) integral representations for general stable density/CDF.
* Chambers-Mallows-Stuck simulation, with an explicit stable `alpha=1` limit
  rather than evaluating `tan(pi/2)` numerically.
* Upstream S0/S1/S2 location/scale transformations.
* Exact Gaussian (`alpha=2`) and Cauchy (`alpha=1,beta=0`) reductions.
* Exact Levy reduction for `alpha=1/2, |beta|=1`; this identity is already
  present in the upstream package's `inst/xtraR/Levy.R` and tests.
* Upstream Pareto tail constants/approximations.

The supplied MIT-licensed `r_mod.f90` is used for normal/Cauchy helpers, RNG,
gamma/log-gamma functions, and numerical integration. The build copy differs
from the supplied file only by free-form continuation/whitespace formatting.

Two small numerical helpers were added because the supplied module has no
R-compatible counterparts: a bracketed scalar root solver and bounded scalar
maximizer. The package-specific integration wrapper still calls `r_mod`'s
`integrate`; it only divides difficult finite intervals into panels.

## Build

```text
fpm build
fpm test
```

The project uses standard free-form Fortran 2018. `r_mod.F90` is preprocessed
because the supplied module contains its existing optional CPP RNG bridge.

## Validation

The included tests cover:

* exact Normal, Cauchy, and Levy identities;
* independent Nolan S0/S1 reference density/CDF values;
* finite-support one-sided cases;
* quantile/CDF inversion;
* stable-mode and S2 mode-centering;
* log-density/log-probability and upper-tail APIs;
* scalar and recycled vector parameters;
* Monte Carlo checks of S0/S1 random generation, including skew `alpha=1`.

As the upstream manual itself notes, calculations with `alpha` extremely close
to 1 or 0 are numerically much more challenging and should not be assumed to
have uniform high precision.

## License

Code translated from `stabledist` is distributed under GPL-2.0-or-later,
following the upstream `License: GPL (>= 2)` field. The supplied `r_mod.F90`
remains MIT-licensed. See `LICENSE`, `LICENSES/`, and `NOTICE.md`.
