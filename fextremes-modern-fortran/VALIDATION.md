# Validation record

## Environment

- GNU Fortran 14.2.0
- Debian Linux execution environment
- Fortran 2018
- No R installation was available
- No `fpm` executable was available
- No external numerical libraries are required

## Debug configuration

```text
-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror
-fcheck=all -fbacktrace
```

Command:

```sh
./build.sh debug
```

Output:

```text
Profile-likelihood and stability tests passed.
Exploratory and tail-index tests passed.
Distribution, fitting, and risk tests passed.
Preprocessing and extremal-index tests passed.
GPL-2.0-or-later source license checks passed.
debug build, tests, and applications passed.
```

## Optimized configuration

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace
```

Command:

```sh
./build.sh release
```

Output:

```text
Profile-likelihood and stability tests passed.
Exploratory and tail-index tests passed.
Distribution, fitting, and risk tests passed.
Preprocessing and extremal-index tests passed.
GPL-2.0-or-later source license checks passed.
release build, tests, and applications passed.
```

## Numerical coverage

The regression suites exercise:

- GEV and GPD distribution inversion and known moments
- Scalar and vector simulation
- GEV PWM and MLE, Gumbel PWM and MLE, and GPD PWM and MLE
- Observed and expected covariance paths
- GPD VaR/ES, tail survival, profile intervals, and threshold stability
- GEV delta-method and profile-likelihood return levels
- Block maxima, threshold selection, point-process extraction, and declustering
- Max-Frechet extremal-index simulation and block, cluster, run, and
  Ferro-Segers estimators
- Empirical survival and Pareto QQ curves
- Mean-excess and mean-residual-life calculations
- Record, subsample-record, maximum/sum, SLLN, and LIL paths
- Exceedance height/distance autocorrelations
- Pickands, Hill, and DEH tail-index calculations
- Normal mean excess, EMA, and RiskMetrics volatility
- Sample lower/upper VaR and CVaR
- Demo, CSV GPD fit, and the Danish claims tail-analysis example

## Equivalence limits

R was unavailable, so exact R-output comparison was not performed. Closed-form
identities, independently known moments, simulation recovery, monotonicity,
finite covariance checks, profile-interval containment, and exact deterministic
array cases were used instead.

Bounded Nelder-Mead replaces R's `optim`, and direct grid profiling replaces R's
spline-smoothed presentation. Exact optimizer iterations, random streams, and
floating-point endpoints are not claimed.
