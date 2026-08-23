# Porting notes

## Source basis

This port is based on the attached `ecpdist` 0.2.1 source package. The complete
supplied archive is retained at `upstream/ecpdist-master.zip`.

The R package has `NeedsCompilation: no`; all substantive numerical code is in
R and there are no external numerical package dependencies beyond base R's
`stats::integrate()`.

## R-to-Fortran API mapping

| R routine | Fortran routine | Status |
|---|---|---|
| `decp` | `decp` | translated |
| `pecp` | `pecp` | translated |
| `qecp` | `qecp` | translated, option bugs corrected |
| `recp` | `recp` | translated |
| `secp` | `secp` | translated |
| `hecp` | `hecp` | translated |
| `ecp_kmoment` | `ecp_kmoment` | translated |
| `ecp_kmoment_cond` | `ecp_kmoment_cond` | translated |
| `ecp_mrl` | `ecp_mrl` | translated |
| `ecp_shape` | `ecp_shape` | translated |
| `ecp_plot` | - | plotting intentionally omitted |

## Corrected upstream `qecp` behavior

The upstream implementation computes the quantile before it executes

```r
if (!lower_tail) p <- 1 - p
```

so `lower_tail=FALSE` has no effect. It also executes `qf <- log(qf)` when
`log_p=TRUE`, whereas the conventional `q*` API meaning of `log.p` is that the
input probability is supplied on the log scale. The Fortran `qecp` implements
the intended standard semantics:

- `lower_tail=.false.` interprets the supplied probability as an upper-tail
  probability;
- `log_p=.true.` interprets the input as `log(p)`.

The underlying lower-tail quantile formula is unchanged.

## Stable probability evaluation

The upstream survival expression is

```text
[1 - exp(-phi*A)] / [1 - exp(-phi)]
```

with

```text
A = exp(lambda * [1 - exp(x^gamma)]).
```

The Fortran code evaluates the ratio on the log scale using stable `expm1`
logic. This avoids overflow for large negative `phi` and cancellation when the
survival probability is close to one.

The hazard is evaluated as `exp(log(pdf) - log(survival))`.

## Moment integration

The R package changes variables and integrates over `y in [0,1]`. Its
integrands contain powers of `log(1 - log(y)/lambda)` and are singular as
`y -> 0` even though the integral is finite.

The Fortran port uses the equivalent quantile identities

```text
E[X^k] = integral_0^1 Q(u)^k du
E[X^k | X>x] = integral_F(x)^1 Q(u)^k du / S(x)
```

and adaptive 15-point Gauss-Kronrod quadrature. The quadrature nodes never
sample the singular endpoint `u=1`. `ecp_integral_result` returns the estimate,
an absolute error estimate, and a status code.

## Invalid arguments

R signals errors. Elemental Fortran functions instead return IEEE NaN for
invalid numerical arguments so they remain usable in array expressions.
`recp` and the numerical integration routines expose integer status codes.

## Validation

The clean source tree was tested with:

```sh
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all
```

The tests include independent high-precision reference values, CDF derivative
checks, CDF/quantile inversion, negative-`phi` cases, the corrected quantile
options, numerical moments/conditional moments, shape measures, and RNG/PIT
checks.
