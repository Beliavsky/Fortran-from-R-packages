! SPDX-License-Identifier: GPL-3.0-only
! Derived from R package ftsa 6.7 by Rob Hyndman and Han Lin Shang.
program test_covariance
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ftsa, only : dp, ftsa_facf, ftsa_long_run_covariance
   implicit none
   real(dp) :: rows(6, 2), cols(2, 6), bandwidth
   real(dp), allocatable :: rho(:), cov(:, :)
   integer :: i, status

   do i = 1, 6
      rows(i, 1) = real(i, dp)
      rows(i, 2) = 2.0_dp * real(i, dp) + 1.0_dp
   end do
   cols = transpose(rows)
   call ftsa_facf(rows, [0, 1, 2], rho, status)
   call check(status == 0 .and. size(rho) == 3, 'facf status')
   call check(abs(rho(1)) < 1.0e-15_dp .and. rho(2) > 0.0_dp, 'facf values')
   call ftsa_long_run_covariance(cols, cov, bandwidth, status=status)
   call check(status == 0 .and. all(shape(cov) == [2, 2]), 'long-run shape')
   call check(maxval(abs(cov - transpose(cov))) < 1.0e-10_dp, 'long-run symmetry')
   call check(ieee_is_finite(bandwidth) .and. bandwidth > 0.0_dp, 'bandwidth')

contains
   subroutine check(condition, message)
      logical, intent(in) :: condition !! Assertion condition expected to be true.
      character(len=*), intent(in) :: message !! Diagnostic emitted if the assertion fails.
      if (.not. condition) error stop message
   end subroutine check
end program test_covariance
