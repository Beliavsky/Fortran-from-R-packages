# Porting notes

## Source package

- Package: `riskParityPortfolio`
- Source version: 0.2.2.9000
- Source date: 2021-05-31
- Original authors: Ze Vinicius and Daniel P. Palomar
- Original license: GPL-3
- Port license: GPL-3.0-only

The original metadata and computational R/C++ files are retained under
`original/` for provenance.

## Feature mapping

| R/C++ area | Fortran implementation |
| --- | --- |
| `riskParityPortfolio()` | `risk_parity_portfolio()` |
| `riskParityPortfolioSCA()` | `risk_parity_sca()` |
| diagonal solution | `diagonal_risk_parity()` |
| Spinu CCD | `risk_parity_ccd_spinu()` |
| Roncalli CCD | `risk_parity_ccd_roncalli()` |
| Choi CCD | `risk_parity_ccd_choi()` |
| Newton-Nesterov | `risk_parity_newton()` |
| active risk parity | `active_risk_parity_ccd()` |
| eight `R_*`, `g_*`, `A_*` families | `risk_objective`, `risk_vector`, `risk_jacobian`, `risk_gradient` |
| `projectBudgetLineAndBox()` | `project_budget_box()` |
| equality/inequality QP helpers | `solve_equality_qp`, `solve_qp_active_set` |
| feasibility/projection helpers | `is_feasible`, `project_linear_constraints` |

## Deliberate implementation changes

1. RcppEigen operations were replaced by native Fortran arrays and BLAS/LAPACK.
2. `quadprog` and the package's accelerated dual QP iteration were replaced by
   a self-contained primal active-set QP method.
3. Deprecated `alabama` and `slsqp` front ends were not ported. Their role is
   covered by the SCA implementation, which the R package itself recommends.
4. R scalar recycling is replaced by explicit length-`n` vectors for bounds and
   budgets.
5. R list outputs are represented by `type(risk_parity_result)`.
6. Plotting and data-frame reshaping are omitted.

## Corrected derivative inconsistency

For `rc-over-sd vs b-times-sd`, the closed-form gradient in the R source does
not agree with the derivative of its own objective for a nonuniform budget
vector. The Fortran port evaluates the analytical identity

```text
gradient R = 2 * transpose(A) * g
```

using the formulation's analytical Jacobian `A` and residual vector `g`.
Finite-difference tests verify the result. No objective definition was changed.

## Numerical details

- The Spinu and Roncalli positive quadratic roots use an algebraically
  equivalent cancellation-resistant form.
- The active-risk-parity update retains the original expression exactly because
  its finite-tradeoff behavior is part of the source algorithm.
- Constraint matrices follow the source convention `C*w=c`, `D*w<=d`.
- The high-level API appends the budget equality and box inequalities.
- Dense systems are appropriate for the package's intended portfolio sizes;
  sparse Matrix-class dispatch was not reproduced.
