! SPDX-License-Identifier: GPL-3.0-only
program test_unified
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use garchito, only : dp, garchito_control, garchito_result, &
      garchito_success, garchito_max_iterations, unified_est
   implicit none

   integer, parameter :: n = 120
   real(dp) :: rv(n), returns(n), h
   integer :: i
   type(garchito_control) :: control
   type(garchito_result) :: fit

   returns = [(0.012_dp*sin(0.31_dp*real(i, dp)) + &
      0.006_dp*cos(0.11_dp*real(i, dp)), i = 1, n)]
   h = 2.0e-5_dp / (1.0_dp - 0.22_dp - 0.60_dp)
   rv(1) = h * (1.0_dp + 0.12_dp*sin(0.7_dp))
   do i = 2, n
      h = 2.0e-5_dp + 0.60_dp*h + 0.22_dp*returns(i - 1)**2
      rv(i) = h * (1.0_dp + 0.12_dp*sin(0.7_dp*real(i, dp)))
   end do

   control%max_iterations = 2200
   control%max_evaluations = 50000
   control%tolerance = 1.0e-8_dp
   call unified_est(rv, returns, fit, control)

   call check(fit%convergence == garchito_success .or. &
      fit%convergence == garchito_max_iterations, 'optimizer status')
   call check(size(fit%coefficients) == 3, 'coefficient count')
   call check(size(fit%sigma) == n, 'sigma length')
   call check(all(fit%sigma > 0.0_dp), 'positive sigma')
   call check(fit%pred > 0.0_dp, 'positive forecast')
   call check(ieee_is_finite(fit%objective), 'finite objective')
   call check(fit%coefficients(2) + fit%coefficients(3) < 1.0_dp, &
      'stationarity constraint')
   call check(abs(fit%coefficients(1) - 2.0e-5_dp) < 1.0e-5_dp, &
      'omega recovery')
   call check(abs(fit%coefficients(2) - 0.22_dp) < 0.10_dp, &
      'beta recovery')
   call check(abs(fit%coefficients(3) - 0.60_dp) < 0.10_dp, &
      'gamma recovery')

   print '(a)', 'test_unified: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,a)') 'FAIL: ', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_unified
