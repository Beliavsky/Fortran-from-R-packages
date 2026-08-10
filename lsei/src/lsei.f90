! SPDX-License-Identifier: GPL-2.0-or-later
module lsei
   use lsei_kinds, only : dp
   use lsei_types
   use lsei_nnls, only : nnls_solve, pnnls_solve, ldp_solve
   use lsei_solver, only : lsi_solve, lsei_solve, qp_solve, pnnqp_solve, hfti_solve
   use lsei_utils, only : indx, mat_maxs
   implicit none
   public
end module lsei
