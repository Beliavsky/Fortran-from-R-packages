! SPDX-License-Identifier: GPL-3.0-only
program demo_garchito
   use garchito, only : dp, garchito_control, garchito_result, &
      unified_est, realized_est
   implicit none

   integer, parameter :: n = 80
   real(dp) :: rv(n), returns(n), jv(n), h
   integer :: i
   type(garchito_control) :: control
   type(garchito_result) :: unified, realized

   returns = [(0.012_dp*sin(0.28_dp*real(i, dp)), i = 1, n)]
   jv = [(8.0e-6_dp*(1.0_dp + 0.25_dp*cos(0.19_dp*real(i, dp)))**2, &
      i = 1, n)]
   h = 1.5e-5_dp / (1.0_dp - 0.25_dp - 0.58_dp)
   rv(1) = h
   do i = 2, n
      h = 1.5e-5_dp + 0.58_dp*h + 0.25_dp*returns(i - 1)**2
      rv(i) = h * (1.0_dp + 0.06_dp*sin(0.41_dp*real(i, dp)))
   end do

   control%max_iterations = 2200
   call unified_est(rv, returns, unified, control)
   call realized_est(rv, realized, jv, control)

   print '(a)', 'GARCHIto-fortran demo'
   print '(a,es14.6)', 'Unified forecast: ', unified%pred
   print '(a,es14.6)', 'Realized forecast: ', realized%pred
   print '(a,i0)', 'Unified evaluations: ', unified%evaluations
   print '(a,i0)', 'Realized evaluations: ', realized%evaluations
end program demo_garchito
