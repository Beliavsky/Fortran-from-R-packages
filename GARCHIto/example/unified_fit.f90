! SPDX-License-Identifier: GPL-3.0-only
program unified_fit
   use garchito, only : dp, garchito_control, garchito_result, unified_est
   implicit none

   integer, parameter :: n = 100
   real(dp) :: rv(n), returns(n), h
   integer :: i
   type(garchito_control) :: control
   type(garchito_result) :: fit

   returns = [(0.01_dp*sin(0.25_dp*real(i, dp)), i = 1, n)]
   h = 2.0e-5_dp / (1.0_dp - 0.25_dp - 0.60_dp)
   rv(1) = h
   do i = 2, n
      h = 2.0e-5_dp + 0.60_dp*h + 0.25_dp*returns(i - 1)**2
      rv(i) = h * (1.0_dp + 0.08_dp*cos(0.43_dp*real(i, dp)))
   end do

   control%max_iterations = 2200
   call unified_est(rv, returns, fit, control)
   call print_result(fit)

contains

   subroutine print_result(result)
      type(garchito_result), intent(in) :: result
      integer :: j
      print '(a,a)', 'status: ', trim(result%message)
      do j = 1, size(result%coefficients)
         print '(a16,2x,es16.8)', trim(result%coefficient_names(j)), &
            result%coefficients(j)
      end do
      print '(a,es16.8)', 'one-step forecast: ', result%pred
   end subroutine print_result

end program unified_fit
