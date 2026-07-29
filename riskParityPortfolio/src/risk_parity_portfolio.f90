! SPDX-License-Identifier: GPL-3.0-only
! Convenience umbrella module.
module risk_parity_portfolio_mod
   use rpp_kinds, only: dp, rpp_huge
   use rpp_api
   use rpp_qp, only: project_budget_box, project_equality, project_feasible, project_linear_constraints, is_feasible, &
                     solve_equality_qp, solve_qp_active_set
   implicit none
   public
end module risk_parity_portfolio_mod
