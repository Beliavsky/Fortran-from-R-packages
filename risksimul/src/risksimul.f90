! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from riskSimul 0.1.2 by Wolfgang Hormann and Ismail Basoglu.
module risksimul
   use ghyp_kinds, only : dp, i8
   use risksimul_types
   use risksimul_math, only : student_quantile, gamma_cdf, gamma_quantile, &
      orthogonal_completion, optimal_allocation_heuristic
   use risksimul_portfolio
   use risksimul_simulation
   implicit none
   public
end module risksimul
