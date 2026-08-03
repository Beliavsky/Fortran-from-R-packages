! SPDX-License-Identifier: GPL-3.0-only
program option_fit
   use garchito, only : dp, garchito_control, garchito_result, &
      realized_est_option
   implicit none

   integer, parameter :: n = 80
   real(dp) :: rv(n), jv(n), nv(n), h
   integer :: i
   type(garchito_control) :: control
   type(garchito_result) :: fit

   jv = [(1.0e-5_dp*(1.0_dp + 0.3_dp*sin(0.21_dp*real(i, dp)))**2, &
      i = 1, n)]
   h = (1.0e-5_dp + 0.30_dp*sum(jv)/real(n, dp)) / &
      (1.0_dp - 0.20_dp - 0.65_dp)
   rv(1) = h
   do i = 2, n
      h = 1.0e-5_dp + 0.65_dp*h + 0.20_dp*rv(i - 1) + &
          0.30_dp*jv(i - 1)
      rv(i) = h * (1.0_dp + 0.10_dp*sin(0.37_dp*real(i, dp)))
   end do
   nv = 5.0e-6_dp + 1.10_dp*rv + &
      [(2.0e-6_dp*cos(0.53_dp*real(i, dp)), i = 1, n)]

   control%max_iterations = 3000
   call realized_est_option(rv, nv, fit, jv=jv, homogeneous=.false., &
      control=control)
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

end program option_fit
