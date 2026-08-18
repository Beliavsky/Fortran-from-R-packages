# Translation notes

## Upstream

- Package: `bivgeom`
- Version: 1.0
- Author: Alessandro Barbiero
- Model: Roy (1993), bivariate geometric distribution
- Upstream license declaration: `GPL`

The original source tree is retained under `upstream/`.

## API mapping

| R routine | Fortran routine |
| --- | --- |
| `dbivgeomRoy` | `dbivgeom_roy`, `dbivgeomroy` |
| `FbivgeomRoy` | `fbivgeom_roy`, `fbivgeomroy` |
| `SbivgeomRoy` | `sbivgeom_roy`, `sbivgeomroy` |
| `FyxbivgeomRoy` | `fyxbivgeom_roy`, `fyxbivgeomroy` |
| `EyxbivgeomRoy` | `eyxbivgeom_roy`, `eyxbivgeomroy` |
| `lambda1Roy` | `lambda1_roy`, `lambda1roy` |
| `lambda2Roy` | `lambda2_roy`, `lambda2roy` |
| `corbivgeomRoy` | `corbivgeom_roy`, `corbivgeomroy` |
| `RelbivgeomRoy` | `relbivgeom_roy`, `relbivgeomroy` |
| `rbivgeomRoy` | `rbivgeom_roy`, `rbivgeomroy` |
| `S.n` | `empirical_survival_roy`, `s_n` |
| `loglikgeomRoy` | `negative_loglik_roy`, `loglikgeomroy` |
| `minuslogRoy` | `minuslogroy` |
| `estbivgeomRoy(...,"ML")` | `fit_bivgeom_ml`, `estbivgeom_roy` |
| `estbivgeomRoy(...,"LS")` | `estimate_ls_roy`, `estbivgeom_roy` |
| `MMP`, `MM1`--`MM4` | corresponding `estimate_*_roy` routines |

## Numerical and interface differences

### Parameter domain

The upstream feasibility checks require `0 < theta1 < 1`, `theta2 > 0`, and
`0 < theta3 <= 1`, but accidentally omit `theta2 < 1`. Since the second
marginal is geometric with survival `theta2^y`, the mathematically valid
model requires `0 < theta2 < 1`. `feasible_roy` enforces this condition.

The Roy feasibility constraint

    theta3 >= (theta1 + theta2 - 1) / (theta1*theta2)

is retained, together with `theta3 > 0`.

### Random generation

The R routine samples `X` geometrically and then calls `uniroot()` on a
continuous extension of the conditional CDF before applying `ceiling()`.
The Fortran routine samples the same geometric marginal and performs exact
integer conditional-quantile inversion by bracket expansion plus binary
search. This removes root-rounding ambiguity while targeting the same
discrete conditional law.

### Empirical survival and the `copula` import

The R `S.n()` helper evaluates the empirical bivariate survival function via
an empirical CDF helper (`F.n`) supplied through the R dependency stack. In
Fortran, the same quantity is computed directly as

    count(X_i >= x and Y_i >= y) / n.

Therefore the translation is self-contained and does not require `copula`.

### Maximum likelihood

The R package delegates bounded optimization to `bbmle::mle2()` with
L-BFGS-B. The Fortran translation uses a native three-parameter bounded
Nelder-Mead optimizer with explicit feasibility projection/penalties. The
likelihood itself is translated directly and accumulated in log space.

### Least squares

The R `lm()` call is replaced by direct construction and solution of the
3-by-3 normal equations for

    log(S_n(x,y)) = beta1*x + beta2*y + beta3*x*y,

with `theta = exp(beta)`.

### Correlation and reliability

The finite truncation rule from the R source is preserved by default: the
series is truncated at twice the `1-alpha` geometric marginal quantiles,
with `alpha=1e-5`. An optional `alpha` argument is exposed in Fortran.

## Validation

The packaged tests check fixed upstream formula values, CDF/PMF identities,
normalization, conditional moments, direct double-sum correlation and
reliability, simulation moments, all estimator families, and ML likelihood
improvement.
