! SPDX-License-Identifier: GPL-2.0-or-later
module coneproj
   use coneproj_kinds, only : dp
   use coneproj_types
   use coneproj_core, only : cone_a, cone_b, qprog
   use coneproj_shape
   use coneproj_regression, only : constreg_fit, shapereg_fit
   use coneproj_stats, only : seed_rng, regularized_beta, student_t_cdf, student_t_quantile
   use coneproj_linalg, only : matrix_rank, column_basis
   implicit none
   public

contains

   subroutine qr_decomp(x, result)
      real(dp), intent(in) :: x(:,:)
      type(qr_result), intent(out) :: result
      call column_basis(x, result%q, result%rank, 1.0e-8_dp)
   end subroutine qr_decomp

end module coneproj
