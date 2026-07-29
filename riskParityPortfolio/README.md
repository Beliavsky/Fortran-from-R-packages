# riskParityPortfolio-fortran

A modern Fortran 2018/FPM port of the computational core of the R package
`riskParityPortfolio` 0.2.2.9000.

The library designs risk-budgeting and risk-parity portfolios. It includes fast
convex solvers for the vanilla long-only problem and a successive convex
approximation (SCA) solver for additional bounds, linear constraints, expected
return terms, and variance penalties.

## Implemented features

- Diagonal-covariance analytical risk-parity weights.
- Spinu cyclical coordinate descent.
- Roncalli cyclical coordinate descent.
- Choi and Chen cyclical coordinate descent.
- Newton-Nesterov solver for the logarithmic formulation.
- Roncalli active-risk-parity coordinate descent.
- All eight risk-concentration formulations from the R package.
- Analytical objective gradients and risk-vector Jacobians.
- Budget-only, box, linear equality, and linear inequality constraints.
- SCA optimization with an internal active-set quadratic-programming solver.
- Exact projection onto the budget line with box constraints.
- General linear-constraint projection and feasibility checks.
- Expected-return rewards and portfolio-variance penalties.

The port does not require R, Rcpp, Eigen, `quadprog`, `nloptr`, or `alabama`.
LAPACK and BLAS are used for dense linear systems.

## Build with FPM

```text
fpm build
fpm test
fpm run
fpm run --example constrained_portfolio
fpm run --example active_risk_parity
```

The system must provide linkable LAPACK and BLAS libraries. The included
`fpm.toml` requests `lapack` and `blas`.

## Minimal example

```fortran
program example
   use risk_parity_portfolio_mod
   implicit none

   real(dp) :: sigma(3, 3), budgets(3)
   type(risk_parity_result) :: result

   sigma = reshape([ &
      0.04_dp, 0.01_dp, 0.00_dp, &
      0.01_dp, 0.09_dp, 0.02_dp, &
      0.00_dp, 0.02_dp, 0.16_dp], [3, 3])
   budgets = [0.50_dp, 0.30_dp, 0.20_dp]

   call risk_parity_portfolio(sigma, result, b=budgets)
   if (result%status /= RPP_OK) error stop "solver failed"

   print *, result%weights
   print *, result%relative_risk_contribution
end program example
```

For constrained portfolios, pass `lower`, `upper`, `cmat`, `cvec`, `dmat`, and
`dvec`. The high-level routine always adds the budget equation `sum(w)=1`.
User-supplied equality constraints therefore should not repeat that equation.
Inequalities use the convention `dmat*w <= dvec`.

## Source organization

- `src/rpp_core.f90`: portfolio risk quantities and basic objectives.
- `src/rpp_formulations.f90`: eight concentration objectives and derivatives.
- `src/rpp_solvers.f90`: fast vanilla and active-risk solvers.
- `src/rpp_qp.f90`: projections and active-set QP solver.
- `src/rpp_api.f90`: high-level result type and SCA interface.
- `src/risk_parity_portfolio.f90`: umbrella module.
- `app/`, `example/`, and `test/`: runnable programs and tests.
- `original/`: original metadata and computational R/C++ source for provenance.

See `API.md`, `PORTING.md`, and `TESTING.md` for details.

## License

This port preserves the original package's GNU General Public License version
3 designation. It is distributed under `GPL-3.0-only`. See `LICENSE` and
`NOTICE.md`.
