! SPDX-License-Identifier: GPL-3.0-only
program realized_fit
   use garchito, only : dp, garchito_control, garchito_result, realized_est
   implicit none

   integer, parameter :: n = 100
   real(dp) :: rv(n), jv(n), h
   integer :: i
   type(garchito_control) :: control
   type(garchito_result) :: fit

   jv = [(8.0e-6_dp*(1.0_dp + 0.35_dp*sin(0.17_dp*real(i, dp)))**2, &
      i = 1, n)]
   h = (1.2e-5_dp + 0.25_dp*sum(jv)/real(n, dp)) / &
      (1.0_dp - 0.18_dp - 0.67_dp)
   rv(1) = h
   do i = 2, n
      h = 1.2e-5_dp + 0.67_dp*h + 0.18_dp*rv(i - 1) + &
          0.25_dp*jv(i - 1)
      rv(i) = h * (1.0_dp + 0.08_dp*sin(0.51_dp*real(i, dp)))
   end do

   control%max_iterations = 2500
   call realized_est(rv, fit, jv, control)
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

end program realized_fit
