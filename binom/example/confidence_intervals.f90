! SPDX-License-Identifier: GPL-2.0-or-later
program confidence_intervals
   use binom
   implicit none
   type(binom_interval) :: intervals(11)
   type(binom_design_result) :: design
   type(binom_optim_result) :: fit
   integer :: i
   intervals = binom_confint_all(3.0_dp, 10.0_dp)
   do i = 1, size(intervals)
      if (intervals(i)%status /= binom_ok) error stop 'Interval calculation failed'
      print '(a13,3f12.7)', binom_methods(i), intervals(i)%mean, intervals(i)%lower, intervals(i)%upper
   end do
   print '(a,f10.6)', 'Wilson coverage at p=.3, n=10: ', binom_coverage(0.3_dp, 10, 'wilson')
   print '(a,f10.6)', 'Wilson expected width: ', binom_length(0.3_dp, 10, 'wilson')
   design = cloglog_sample_size(p_alt=0.7_dp, p=0.4_dp, recompute_power=.true.)
   if (design%status /= binom_ok) error stop 'Sample-size calculation failed'
   print '(a,f9.0,a,f10.6)', 'Cloglog sample size: ', design%n, ' attained power: ', design%power
   call binom_optim(5, fit, method='wilson', k=1)
   if (fit%convergence /= 0) error stop 'Coverage optimization failed'
   print '(a,f12.8,a,es12.4)', 'Optimized endpoint confidence: ', fit%par(1), ' error integral: ', fit%value
end program confidence_intervals
