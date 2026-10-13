! SPDX-License-Identifier: BSD-3-Clause
program diagnose_draws
   use posterior, only: dp, rhat, ess_bulk, ess_tail, mcse_mean, pareto_diags, pareto_diagnostic
   implicit none
   real(dp) :: x(100, 4), excess(100)
   type(pareto_diagnostic) :: tail
   integer :: i, j
   do j = 1, 4
      do i = 1, 100
         x(i, j) = sin(real(i * (2 * j + 1), dp)) + cos(real(i * 7 + j, dp))
      end do
   end do
   excess = [(log(101.0_dp / real(i, dp)), i = 1, 100)]
   tail = pareto_diags(excess, tail = 'right', r_eff = 1.0_dp)
   print '(a,f10.5)', 'Rank-normalized split Rhat: ', rhat(x)
   print '(a,f10.3)', 'Bulk ESS: ', ess_bulk(x)
   print '(a,f10.3)', 'Tail ESS: ', ess_tail(x)
   print '(a,f10.5)', 'Mean MCSE: ', mcse_mean(x)
   print '(a,f10.5)', 'Right-tail Pareto shape: ', tail%khat
end program diagnose_draws
