! SPDX-License-Identifier: GPL-2.0-or-later
module nnls
   use nnls_kinds, only : dp
   use nnls_solver, only : nnls_result, nnls_fit, nnnpls_fit, nnls_solve, &
      NNLS_SUCCESS, NNLS_BAD_DIMENSIONS, NNLS_ITERATION_LIMIT
   implicit none
   public
end module nnls
