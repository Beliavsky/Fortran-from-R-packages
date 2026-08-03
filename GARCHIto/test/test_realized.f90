! SPDX-License-Identifier: GPL-3.0-only
program test_realized
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use garchito, only : dp, garchito_control, garchito_result, &
      garchito_success, garchito_max_iterations, garchito_invalid_input, &
      realized_est
   implicit none

   integer, parameter :: n = 100
   real(dp) :: rv(n), jv(n), h, med_jv
   integer :: i
   type(garchito_control) :: control
   type(garchito_result) :: fit

   control%max_iterations = 2400
   control%max_evaluations = 60000
   control%tolerance = 1.0e-8_dp

   h = 1.5e-5_dp / (1.0_dp - 0.25_dp - 0.58_dp)
   rv(1) = h * 1.05_dp
   do i = 2, n
      h = 1.5e-5_dp + 0.58_dp*h + 0.25_dp*rv(i - 1)
      rv(i) = h * (1.0_dp + 0.10_dp*cos(0.47_dp*real(i, dp)))
   end do
   call realized_est(rv, fit, control=control)
   call valid_fit(fit, 3, n, 3)

   jv = [(1.0e-5_dp*(1.0_dp + 0.4_dp*sin(0.19_dp*real(i, dp)))**2, &
      i = 1, n)]
   med_jv = 1.0e-5_dp
   h = (1.0e-5_dp + 0.30_dp*med_jv) / (1.0_dp - 0.20_dp - 0.65_dp)
   rv(1) = h * (1.0_dp + 0.08_dp*cos(0.4_dp))
   do i = 2, n
      h = 1.0e-5_dp + 0.65_dp*h + 0.20_dp*rv(i - 1) + &
          0.30_dp*jv(i - 1)
      rv(i) = h * (1.0_dp + 0.08_dp*cos(0.4_dp*real(i, dp)))
   end do
   call realized_est(rv, fit, jv, control)
   call valid_fit(fit, 4, n, 4)

   rv(4) = -1.0_dp
   call realized_est(rv, fit, jv, control)
   call check(fit%convergence == garchito_invalid_input, &
      'negative rv rejected')

   print '(a)', 'test_realized: PASS'

contains

   subroutine valid_fit(result, ncoef, nobs, gamma_index)
      type(garchito_result), intent(in) :: result
      integer, intent(in) :: ncoef, nobs, gamma_index
      call check(result%convergence == garchito_success .or. &
         result%convergence == garchito_max_iterations, 'optimizer status')
      call check(size(result%coefficients) == ncoef, 'coefficient count')
      call check(size(result%coefficient_names) == ncoef, 'coefficient names')
      call check(size(result%sigma) == nobs, 'sigma length')
      call check(all(result%sigma > 0.0_dp), 'positive sigma')
      call check(result%pred > 0.0_dp, 'positive forecast')
      call check(ieee_is_finite(result%objective), 'finite objective')
      call check(result%coefficients(2) + result%coefficients(gamma_index) < 1.0_dp, &
         'stationarity constraint')
   end subroutine valid_fit

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,a)') 'FAIL: ', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_realized
