! SPDX-License-Identifier: GPL-3.0-only
program test_validation
   use garchito, only : dp, garchito_control, garchito_result, &
      garchito_invalid_input, unified_est, realized_est_option
   implicit none

   real(dp) :: rv3(3), rv4(4), ret4(4), nv3(3), jv2(2)
   type(garchito_control) :: control
   type(garchito_result) :: fit

   rv3 = [1.0_dp, 1.1_dp, 0.9_dp]
   rv4 = [1.0_dp, 1.1_dp, 0.9_dp, 1.05_dp]
   ret4 = [0.1_dp, -0.1_dp, 0.05_dp, -0.02_dp]
   nv3 = [1.0_dp, 1.0_dp, 1.0_dp]
   jv2 = [0.1_dp, 0.2_dp]

   call unified_est(rv3, ret4, fit)
   call check(fit%convergence == garchito_invalid_input, 'length mismatch')

   call realized_est_option(rv3, nv3, fit, jv2)
   call check(fit%convergence == garchito_invalid_input, 'jv mismatch')

   control%max_iterations = 1
   control%max_evaluations = 1000
   call unified_est(rv4, ret4, fit, control)
   call check(allocated(fit%coefficients), 'best estimate on iteration limit')
   call check(allocated(fit%sigma), 'sigma on iteration limit')
   call check(fit%iterations <= control%max_iterations, 'iteration count bounded')

   print '(a)', 'test_validation: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,a)') 'FAIL: ', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_validation
