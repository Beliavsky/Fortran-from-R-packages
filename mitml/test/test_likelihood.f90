! SPDX-License-Identifier: GPL-2.0-or-later
! Upstream mitml 0.4-5 (2023-03-08), authored by Simon Grund,
! Alexander Robitzsch, and Oliver Luedtke; upstream license GPL (>= 2).
! Modern free-form Fortran translation for Fortran-from-R-packages.
program test_likelihood
   use mitml, only : dp, gaussian_lm_loglik, gaussian_lmm_loglik, MITML_ERR_ARGUMENT, MITML_ERR_DIMENSION, MITML_OK
   implicit none

   real(dp) :: beta(2), tau(1, 1), value
   real(dp) :: x(4, 2), y(4), z(4, 1)
   integer :: cluster(4), status

   y = [1.1_dp, 0.7_dp, 2.0_dp, 1.6_dp]
   x(:, 1) = 1.0_dp
   x(:, 2) = [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp]
   beta = [1.0_dp, 0.5_dp]
   z(:, 1) = 1.0_dp
   cluster = [10, 10, 20, 20]
   tau(1, 1) = 0.4_dp

   call gaussian_lm_loglik(y, x, beta, 0.8_dp, value, status)
   call require(status == MITML_OK, 'gaussian_lm_loglik status')
   call require(abs(value - (-4.266967030190271_dp)) < 1.0e-12_dp, 'gaussian_lm_loglik value')

   call gaussian_lmm_loglik(y, x, z, cluster, beta, tau, 0.8_dp, value, status)
   call require(status == MITML_OK, 'gaussian_lmm_loglik status')
   call require(abs(value - (-1.018735077931526_dp)) < 1.0e-12_dp, 'gaussian_lmm_loglik value')

   call gaussian_lm_loglik(y, x, beta, 0.0_dp, value, status)
   call require(status == MITML_ERR_ARGUMENT, 'gaussian_lm_loglik variance guard')

   call gaussian_lmm_loglik(y, x(1:3, :), z, cluster, beta, tau, 0.8_dp, value, status)
   call require(status == MITML_ERR_DIMENSION, 'gaussian_lmm_loglik dimension guard')

   print '(a)', 'test_likelihood: PASS'

contains

   subroutine require(condition, label)
      logical, intent(in) :: condition !! True when the tested invariant is satisfied.
      character(len=*), intent(in) :: label !! Short diagnostic label printed if the invariant fails.

      if (.not. condition) then
         print '(a,1x,a)', 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine require

end program test_likelihood
