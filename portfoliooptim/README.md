# PortfolioOptim for modern Fortran

This project translates the numerical core of the R package `PortfolioOptim`
1.1.1 into a self-contained modern Fortran library managed by FPM.

The original package was written by Andrzej Palczewski, with Aleksandra
Dabrowska listed as a contributor. The Fortran project remains licensed under
GNU GPL version 3 only.

## Scope

The library implements:

- weighted empirical VaR, CVaR, and mean absolute deviation;
- CVaR, deviation-CVaR, lower semi-absolute deviation, and MAD portfolio risk;
- Benders decomposition for scenario portfolio optimization;
- a native two-phase simplex solver replacing `Rsymphony`;
- least-distance selection from the set of risk-optimal portfolios;
- the Zhao-Li path-following projection equations used by the R package;
- a robust simplex-plus-quadratic-projection method used by default;
- lower and upper portfolio bounds;
- a target expected-return constraint; and
- arbitrary additional linear inequality constraints.

The project does not require R, Rsymphony, Rglpk, BLAS, or LAPACK.

## Build

```text
fpm build
fpm test
fpm run portfoliooptim_demo
fpm run --example benders_risk_measures
fpm run --example benchmark_projection
```

A direct GNU Fortran validation script is also included:

```text
sh scripts/validate.sh
```

On Windows:

```bat
scripts\validate.bat
```

## Basic usage

```fortran
use portfoliooptim, only : dp, portfolio_result, bdportfolio_optim, risk_mad

real(dp) :: returns(5, 2), probabilities(5)
real(dp) :: a(2, 2), b(2), lower(2), upper(2)
type(portfolio_result) :: result

! Fill returns by scenario and asset.
probabilities = 0.20_dp

a(1, :) = 1.0_dp
a(2, :) = -1.0_dp
b = [1.0_dp, -1.0_dp]  ! sum(weights) = 1
lower = 0.0_dp
upper = 1.0_dp

result = bdportfolio_optim(returns, probabilities, 0.015_dp, risk_mad, &
  0.95_dp, a, b, lower, upper)
```

`portfolio_result` contains:

- `return_mean`: weighted asset means;
- `theta`: portfolio weights;
- `mu`: realized portfolio mean;
- `var`, `cvar`, and `mad`;
- `risk`: the selected optimization risk;
- `new_portfolio_return`: any adjusted target;
- `iterations`, `converged`, and `message`.

## Risk constants

```fortran
risk_cvar
risk_dcvar
risk_lsad
risk_mad
```

The helper `risk_code(name)` converts the corresponding text names.

## Linear constraints

The public Fortran interfaces interpret additional constraints as documented by
the R package:

```text
Aconstr * theta <= bconstr
```

A budget equality can be represented by two inequalities:

```text
 sum(theta) <= 1
-sum(theta) <= -1
```

The projection routine has an optional `upstream_constraint_sign=.true.` switch
for reproducing the sign convention actually passed to the internal R
path-following routine.

## Projection method

`portfolio_optim_projection` first finds the minimum-risk LP value and then
solves the strictly convex least-distance problem. The default implementation
uses a dense ADMM quadratic projection after the native simplex solve. This is
more stable across compilers than applying the original path-following equations
to degenerate optimal faces.

The translated Zhao-Li routine remains available as `zi_projection`. It can be
selected in the high-level interface with `use_zi=.true.`.

The benchmark is shifted consistently when lower bounds are used. Set
`upstream_benchmark_shift=.true.` to reproduce the original unshifted benchmark
behavior.

## Data layout

Returns are stored as:

```text
returns(number_of_scenarios, number_of_assets)
```

Probabilities are supplied separately rather than appended as the final matrix
column. Nonnegative probabilities are normalized internally.

## Numerical notes

The Benders implementation preserves the package's centered-return cuts and the
nonnegative auxiliary-variable convention used by the upstream LP. It avoids
constructing the full scenario LP and is therefore the preferred interface for
large scenario sets.

The projection interface constructs the full scenario LP and is intended for
moderate sample sizes, matching the package description.

See `PORTING_NOTES.md`, `COVERAGE.md`, and `VALIDATION.md` for detailed mapping,
corrections, and test information.
