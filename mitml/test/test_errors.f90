! SPDX-License-Identifier: GPL-2.0-or-later
! Upstream mitml 0.4-5 (2023-03-08), authored by Simon Grund,
! Alexander Robitzsch, and Oliver Luedtke; upstream license GPL (>= 2).
! Modern free-form Fortran translation for Fortran-from-R-packages.
program test_errors
   use mitml, only : MITML_ERR_ARGUMENT, MITML_ERR_DIMENSION, cluster_means, d2_test, dp, mi_test_result, pool_estimates, &
      pooled_estimates
   implicit none
   real(dp) :: means(2)
   real(dp) :: qhat(1, 1)
   type(mi_test_result) :: test_result
   type(pooled_estimates) :: pool_result
   integer :: status

   qhat = 1.0_dp
   call pool_estimates(qhat, pool_result, uhat=reshape([1.0_dp], [1, 1, 1]))
   if (pool_result%status /= MITML_ERR_ARGUMENT) error stop "FAIL: pooling should require two imputations with uhat"
   call d2_test([1.0_dp], 1, test_result)
   if (test_result%status /= MITML_ERR_ARGUMENT) error stop "FAIL: D2 should require two imputations"
   call cluster_means([1.0_dp, 2.0_dp], [1], means, status)
   if (status /= MITML_ERR_DIMENSION) error stop "FAIL: cluster dimension error"

end program test_errors
