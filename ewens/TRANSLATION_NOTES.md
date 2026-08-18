# Translation notes

## Upstream

- R package: `ewens`
- Upstream version: 0.1.0
- Author: Chris Hanretty
- License: MIT
- Source supplied by the user as `ewens-master.zip`

The supplied `copula-fortran` translation is vendored under
`vendor/copula-fortran` and used as an actual FPM dependency.

## API mapping

| R routine | Fortran routine | Notes |
|---|---|---|
| `dewens` | `dewens`, `dewens_log` | Integer class-label vector |
| internal frequency calculation | `class_size_spectrum` | Explicit reusable helper |
| - | `dewens_counts`, `dewens_counts_log` | Direct `m_j` frequency-spectrum API |
| `dewens_k` | `dewens_k`, `dewens_k_log` | Uses copula Stirling numbers when int64-safe |
| `ewens_k_exact` | `ewens_k_exact` | Included although not exported by upstream NAMESPACE |
| `ewens_mle` | `ewens_mle`, `ewens_mle_nk` | Bisection with adaptive upper bracket |
| `rewens` | `rewens` | Native Fortran CRP sampler |
| `gcrp` | `gcrp` | Pitman-Yor/generalized CRP sampler |
| `rgem` | `rgem` | Beta stick breaking using copula RNG primitives |
| R RNG state | `ewens_seed` | Reproducible dependency RNG seed |

## Numerical details

### Ewens PMF

The log PMF is evaluated directly using `log_gamma()` and sums of logarithms,
avoiding factorial/rising-factorial overflow.

### Number of classes

Upstream uses `copula::Stirling1(n,k)`. The attached Fortran dependency's
`stirling_first()` returns an exact signed-64-bit integer and is safe for all
entries through `n=20`. For `n>20`, the translation uses the same unsigned
Stirling recurrence in log space,

`c(n,k) = c(n-1,k-1) + (n-1)c(n-1,k)`,

so `dewens_k` remains usable after exact integer values exceed `int64`.

### MLE

For interior cases `1 < K < n`, the score equation is solved by bisection
with an adaptive upper bracket. This removes the upstream hard-coded
`uniroot(..., interval=c(1e-6,n))` limitation.

Two boundary cases are treated mathematically:

- `K=1`: the MLE is `theta=0`.
- `K=n`: the likelihood increases toward its supremum as `theta -> infinity`,
  so the routine returns positive IEEE infinity.

For `n=1`, the parameter is unidentifiable and the routine returns NaN.

### theta = 0

Upstream `dewens()` evaluates an `-Inf - (-Inf)` expression at `theta=0` and
can return NaN. The Fortran translation uses the limiting distribution:
probability one for a single class of size `n`, and zero otherwise.

### GEM boundaries

The upstream R check allows `alpha=1`, but `rbeta(1-alpha,...)` then has a
zero shape parameter. The Fortran routine follows the well-defined
Pitman-Yor range `0 <= alpha < 1`. At `theta=-alpha`, the first stick is the
proper degenerate value one.

## R-only infrastructure omitted

Roxygen documentation, package/S3 infrastructure, `.Call` registration, and R
RNG state management are not translated. The computational C algorithms are
implemented directly in Fortran.

## Validation

The release tests cover:

- an independently evaluated Ewens frequency-profile PMF;
- `P(K=1 | n=20, theta=1)=1/20`;
- normalization and exact mean of the `K` distribution at `n=50`, exercising
  the log-Stirling fallback;
- an independent MLE reference (`n=10`, `K=4`);
- MLE boundary cases;
- equality of `rewens` and `gcrp(alpha=0)` under a common seed;
- Monte Carlo agreement with the exact expected number of classes;
- GEM support and first-stick expectation.
