! SPDX-License-Identifier: GPL-3.0-only
program test_option
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use garchito, only : dp, garchito_control, garchito_result, &
      garchito_success, garchito_max_iterations, realized_est_option
   implicit none

   integer, parameter :: n = 60
   real(dp) :: rv(n), jv(n), nv(n), h
   integer :: i
   type(garchito_control) :: control
   type(garchito_result) :: fit

   control%max_iterations = 3000
   control%max_evaluations = 60000
   control%tolerance = 1.0e-7_dp

   jv = [(1.0e-5_dp*(1.0_dp + 0.4_dp*sin(0.19_dp*real(i, dp)))**2, &
      i = 1, n)]
   h = (1.0e-5_dp + 0.30_dp*sum(jv)/real(n, dp)) / &
      (1.0_dp - 0.20_dp - 0.65_dp)
   rv(1) = h * (1.0_dp + 0.12_dp*sin(0.7_dp))
   do i = 2, n
      h = 1.0e-5_dp + 0.65_dp*h + 0.20_dp*rv(i - 1) + &
          0.30_dp*jv(i - 1)
      rv(i) = h * (1.0_dp + 0.12_dp*sin(0.7_dp*real(i, dp)))
   end do
   nv = 4.0e-6_dp + 1.15_dp*rv + &
      [(3.0e-6_dp*sin(0.6_dp*real(i, dp)), i = 1, n)]

   call realized_est_option(rv, nv, fit, homogeneous=.true., control=control)
   call valid_fit(fit, 6, 3)
   call realized_est_option(rv, nv, fit, homogeneous=.false., control=control)
   call valid_fit(fit, 7, 3)
   call realized_est_option(rv, nv, fit, jv, .true., control)
   call valid_fit(fit, 7, 4)
   call realized_est_option(rv, nv, fit, jv, .false., control)
   call valid_fit(fit, 8, 4)

   print '(a)', 'test_option: PASS'

contains

   subroutine valid_fit(result, ncoef, gamma_index)
      type(garchito_result), intent(in) :: result
      integer, intent(in) :: ncoef, gamma_index
      call check(result%convergence == garchito_success .or. &
         result%convergence == garchito_max_iterations, 'optimizer status')
      call check(size(result%coefficients) == ncoef, 'coefficient count')
      call check(size(result%sigma) == n, 'sigma length')
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

end program test_option
