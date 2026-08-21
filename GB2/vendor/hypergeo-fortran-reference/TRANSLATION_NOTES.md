# Translation notes

## Scope

`hypergeo` 1.2-14 is a pure-R package.  The numerical algorithms have been translated
to modern free-format Fortran.  R vector attributes, call inspection, warnings,
option handling, and other language-specific presentation/glue were not reproduced.
There is no plotting code in the numerical translation.

All Fortran files, including the vendored deSolve solver backend, use free source form
and the `.f90` suffix.

## Main R-to-Fortran mapping

| R computational surface | Fortran surface |
| --- | --- |
| `complex_gamma`, `complex_factorial`, `lanczos` | `complex_gamma`, `complex_factorial`, `lanczos` |
| `genhypergeo`, `genhypergeo_series` | same underscore-safe names |
| `genhypergeo_contfrac_single`, `genhypergeo_contfrac` | same names; uses `contfrac-fortran` |
| `hypergeo` | generic `hypergeo` scalar/array interface |
| `hypergeo_powerseries`, `hypergeo_general` | same names |
| `f15.1.1`, `f15.3.1`, `f15.3.3`--`f15.3.14` | `f15_1_1`, `f15_3_1`, `f15_3_3`--`f15_3_14` |
| `i15.3.*`, `j15.3.*` | `i15_3_*`, `j15_3_*` |
| `hypergeo_cover1/2/3` | same names |
| Wolfram `w07.23.06.*` helpers | underscore forms `w07_23_06_*` |
| `hypergeo_gosper` | `hypergeo_gosper` |
| Buhring routines | `lpham`, `buhring_eqn11/12`, `buhring_eqn5_*`, `hypergeo_buhring` |
| `shanks`, `genhypergeo_shanks`, `hypergeo_shanks` | same names |
| residue helpers | same names; use `elliptic-fortran` |
| `hypergeo_press`, `f15.5.1` | `hypergeo_press`, `f15_5_1`; use `desolve` |
| `semicircle`, `semidash`, `straight`, `straightdash` | same names |
| `to_real`, `to_complex` | same names |

R-only `.process_args()` is replaced by explicit scalar and array Fortran interfaces.
`hypergeo_func()` and `complex_ode()` are represented internally by the typed deSolve
callback used by `hypergeo_ode_continue()` and `f15_5_1()`.

## Numerical implementation details

### Continued fractions

The generalized-hypergeometric continued-fraction formula is translated directly:

`alpha_k = z * prod(U+k) / prod(k + c(1,L))`.

In R, `prod(k + c(1,L))` is `(k+1) * prod(k+L)`.  The test suite includes the
upstream Maple regression value
`1.0007289707983569879 + 0.0086250714217251837 i` to guard this indexing detail.
The modified Lentz implementation is supplied by the attached `contfrac-fortran`
translation.

### Euler integral

The upstream `f15.3.1()` delegates numerical integration to `elliptic`.  Direct
endpoint evaluation is problematic when the Euler integrand has an integrable
algebraic singularity at 0 or 1.  The Fortran translation therefore evaluates this
specific integral with tanh-sinh quadrature and forms `t` and `1-t` independently near
the endpoints to avoid cancellation.  The optional complex detour point is retained.

The separate Cauchy/residue methods continue to use the translated `elliptic` package.

### ODE continuation

The Press/A&S 15.5.1 continuation is implemented through the translated free-format
`deSolve` package.  The top-level Fortran `hypergeo()` preserves the main R fallback
sequence (Gosper near critical points, analytic/power-series methods, continued
fraction, Euler integral) and additionally tries ODE continuation as a final fallback.
The ODE methods are also directly callable.

### Complex digamma

Several integer-difference formulas need digamma values.  The Fortran implementation
uses reflection, recurrence, and an asymptotic expansion, allowing complex arguments
without an external special-function library.

### Shanks routine

The upstream `genhypergeo_shanks()` computes the Shanks-accelerated sequence internally
but returns the ordinary final partial sum.  That observable behavior is intentionally
preserved.

### External interfaces omitted

No external CAS/shell interface is needed.  R-specific names/attributes, `match.call`,
options such as `showHGcalls`, and debug return-list structures are omitted or replaced
by explicit Fortran arguments/results.

## Dependency integration

The release vendors:

- `dependencies/contfrac-fortran`
- `dependencies/elliptic-fortran`
- `dependencies/desolve`

The root `fpm.toml` declares all three as path dependencies.  The deSolve dependency
retains its own manifest settings for the converted legacy solver backend; its source
is nevertheless entirely free format.

## Validation

The package tests are compiled with:

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all
```

The classic deSolve backend is compiled as free-format Fortran 2018 with runtime
checking; it naturally emits obsolescence diagnostics for retained COMMON/labeled-DO
constructs when such diagnostics are enabled.

Packaged tests cover:

- upstream six-point Maple regression values;
- upstream integer-difference and terminating-polynomial cases;
- all seven elementary A&S identities used by the R test suite;
- A&S continuation transformations;
- generalized power series and continued fractions;
- complex gamma;
- Gosper, Buhring, Euler-integral, and ODE continuation paths;
- Cauchy/residue evaluation and complex/real packing;
- complex A/B/C/z cases from the upstream SAGE comparison.

In addition to the packaged tests, all 100 complex-parameter rows in upstream
`tests/aab.R` were evaluated against the stored SAGE values.  The maximum absolute
difference in the release validation was approximately `6.9e-13`.
