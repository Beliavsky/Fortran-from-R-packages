# Translation notes

## Upstream

- R package: `nbconv`
- Upstream version: 1.0.1
- Author/copyright holder: Gregory Bedwell
- Upstream license: GPL (>= 3)

The original `DESCRIPTION`, `NAMESPACE`, and `README.md` are retained under
`upstream/` for provenance.

## Coverage

| R routine | Fortran translation |
| --- | --- |
| `nb_sum_exact` | `nb_sum_exact` |
| `nb_sum_moments` | `nb_sum_moments` |
| `nb_sum_saddlepoint` | `nb_sum_saddlepoint` |
| `dnbconv` | `dnbconv`, `dnbconv_mu`, `dnbconv_p` |
| `pnbconv` | `pnbconv`, `pnbconv_mu`, `pnbconv_p` |
| `qnbconv` | `qnbconv`, `qnbconv_mu`, `qnbconv_p` |
| `rnbconv` | `rnbconv`, `rnbconv_mu`, `rnbconv_p` |
| `nbconv_params` | `nbconv_params`, `nbconv_params_mu`, `nbconv_params_p` |

No plotting code exists in the upstream package.

## Numerical implementation details

### Furman series

The formulas in the R implementation are preserved. The Fortran translation
keeps the recurrence and final mixture sum in log space, including a
log-add-exp accumulation for the final PMF. This avoids avoidable underflow in
the original `total <- total + exp(probs)` loop.

The upstream `n.terms` default of 1000 and K-series mass tolerance default of
`1e-3` are preserved. If the truncated K distribution is outside tolerance,
`nb_sum_exact` raises `error stop`, corresponding to the R error asking for
more terms. The optional `enforce_tolerance=.false.` can be used by callers
who explicitly want an unconverged truncation, and optional `k_mass` exposes
the retained K mass.

Components with `p=1` (equivalently `mu=0`) are mathematically degenerate at
zero. They are removed before applying Furman's nondegenerate formulas. This
also makes the all-degenerate convolution return a point mass at zero rather
than entering logarithms of zero.

### Saddlepoint approximation

The same cumulant-generating function and Daniels saddlepoint density formula
are used. R's `uniroot` call is replaced by monotone bisection on the exact
CGF domain

```text
t < min_i log(1 + phi_i / mu_i).
```

The upper bracket approaches that singular endpoint adaptively from below.
The Fortran solve is intentionally tighter than the upstream
`sqrt(.Machine$double.eps)` stopping tolerance; this changes only root-solver
roundoff, not the approximation being computed.

As upstream, `normalize=.true.` normalizes over the supplied `counts` vector.
Therefore callers wanting a normalized saddlepoint approximation should pass
a sufficiently wide contiguous support.

### Moment approximation

The translated moment-matched size parameter is exactly

```text
phi = (sum(mu))^2 / sum(mu^2 / phi_i),
```

followed by the standard negative-binomial PMF with that size and total mean.

### RNG

R's `rnbinom` dependency is replaced by the exact gamma-Poisson mixture for a
negative-binomial random variable. Gamma deviates use Marsaglia-Tsang; Poisson
deviates use direct inversion for small means and transformed rejection for
large means. A Park-Miller stream is included to make tests and examples
reproducible without compiler-specific `random_number` sequences.

### R-only dependencies removed

The upstream imports `parallel`, `matrixStats`, and `stats`. They are not
required by the Fortran translation:

- `matrixStats::logSumExp` -> native log-space helpers
- `stats::dnbinom/rnbinom` -> native negative-binomial PMF/RNG
- `stats::uniroot` -> native saddlepoint bisection
- `parallel::mclapply` -> omitted; Fortran callers can parallelize independent
  count/sample loops with their preferred OpenMP/coarray/application layer

The `n.cores` argument is consequently not part of the Fortran API.

## Quantile convention

The Fortran quantile returns the smallest nonnegative integer `x` with
`F(x) >= p`, the conventional discrete quantile definition. For ordinary
non-boundary probabilities this agrees with the upstream implementation; it
also gives a defined answer when `p < F(0)`, an edge case in which the R
`max(which(cdf <= p) - 1)` expression can yield an empty result.

## Validation

The release tests use:

- fixed PMF/CDF values independently obtained by direct convolution of the two
  constituent negative-binomial distributions;
- fixed moment-matched negative-binomial values;
- fixed high-accuracy saddlepoint values from an independent root solve;
- exact cumulant-derived summary values;
- PMF/CDF/quantile consistency checks;
- degenerate `p=1` cases;
- Monte Carlo mean/variance checks for the native RNG.

The release was compiled with GNU Fortran 14.2 using:

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all
```
