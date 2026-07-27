! SPDX-License-Identifier: GPL-3.0-only
! Based on PortfolioOptim 1.1.1 by Andrzej Palczewski and Aleksandra Dabrowska.
module portfoliooptim
  use portfoliooptim_kinds, only : dp
  use portfoliooptim_types, only : risk_result, lp_result, projection_result, &
    portfolio_result, risk_cvar, risk_dcvar, risk_lsad, risk_mad
  use portfoliooptim_linalg, only : diagonal_matrix
  use portfoliooptim_simplex, only : solve_lp
  use portfoliooptim_risk, only : risk_post, risk_measure, risk_code
  use portfoliooptim_benders, only : bdportfolio_optim
  use portfoliooptim_projection, only : f_func, zi_projection, quadratic_lp_projection, portfolio_optim_projection
  implicit none
  public

  interface PortfolioOptimProjection
    module procedure portfolio_optim_projection
  end interface PortfolioOptimProjection
end module portfoliooptim
