# Porting notes

## Probability kernels

Upstream `mnormt` already contains mature Fortran algorithms for bivariate,
trivariate and general multivariate normal/t probabilities. These algorithms
were retained and converted mechanically from fixed-form source to free-form
Fortran. The modern wrappers use explicit module APIs and no R C/Fortran
registration layer.

Some of the retained Genz compatibility kernels still use standard-but-
obsolescent constructs such as `COMMON`, `ENTRY`, and labeled `DO`. This is
intentional: changing those stateful integration internals would add risk with
no numerical benefit. All newly written project code is modern free-form
Fortran 2018 and compiles cleanly with `-Wall -Wextra -Werror -fcheck=all`.

## Random generation

The original `rtmng.f` calls R's `pnorm`, `qnorm`, and RNG entry points. It is
not used. `rmtruncnorm` is a native implementation of the same componentwise
conditional-normal Gibbs algorithm. Its default starting point is the exact
truncated mean obtained from `mom_mtruncnorm`, matching upstream behavior.

## Moment arrays

R can return arrays whose rank is determined at runtime. Standard Fortran
cannot expose an arbitrary-rank allocatable result as conveniently, so
`recintab` returns the same column-major moment table flattened to rank one.
The associated shape is `kappa + 1`; `raw_moment_at` converts a zero-based
multi-index to the corresponding value. `mom_mtruncnorm` provides a higher-
level result containing normalized raw moments, means, covariance, third and
fourth cumulants, marginal standardized cumulants, and Mardia measures.

## Limits

The upstream adaptive Genz algorithms support at most 20 dimensions. The port
preserves that numerical limit. For multivariate t probabilities, upstream
rounds non-integer degrees of freedom to the nearest integer; the port does the
same for dimensions greater than one while retaining real-valued df for the
univariate t density/CDF and for `dmt`/`rmt`.

## Validation

Tests include independent SciPy/quadrature reference values for correlated
normal probabilities and correlated truncated-normal first/second moments,
plus normal/t specialization, specialized-vs-general probability consistency,
CDF normalization, random generation, and cumulant identities.


## FPM implicit-interface settings (v0.1.1)

Current FPM defaults disable implicit typing and implicit external interfaces. The retained `biv_nt.f90`, `tvpack.f90`, and `sadmvnt.f90` kernels intentionally preserve legacy Genz code that relies on those language features. Therefore `fpm.toml` explicitly sets:

```toml
[fortran]
implicit-typing = true
implicit-external = true
source-form = "free"
```

Without these entries, recent FPM/GFortran builds can fail with `-Werror=implicit-interface` at calls such as `MVPHI` in `biv_nt.f90`.
