! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
! Public umbrella module for parma-fortran.
module parma
   use parma_kinds
   use parma_types
   use parma_rng
   use parma_helpers
   use parma_linalg
   use parma_risk
   use parma_timeseries
   use parma_constraints
   use parma_utility
   use parma_cmaes
   use parma_qp
   use parma_lp
   use parma_milp
   use parma_socp
   use parma_portfolio
   use parma_methods
   implicit none
   public
end module parma
