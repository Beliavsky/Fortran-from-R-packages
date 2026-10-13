! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program permutation_examples
   use, intrinsic :: iso_fortran_env, only : int64
   use perm_api
   implicit none
   type(perm_result) :: result
   type(perm_control_type) :: control
   type(perm_rng) :: rng
   real(dp) :: scores(6)
   integer :: status
   scores = [1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, 8.0_dp]
   call perm_ts(scores(:3), scores(4:), result, method='exact.network')
   if (result%status /= perm_ok) error stop 'Exact two-sample example failed'
   print '(a,f10.6)', 'Exact two-sided p: ', result%p_value
   print '(a,f10.6)', 'Difference of means: ', result%estimate
   call perm_control(control, status, nmc=19999, seed=20261011_int64, confidence=0.99_dp)
   if (status /= perm_ok) error stop 'Control construction failed'
   call perm_ks(scores, [1, 1, 2, 2, 3, 3], result, method='exact.mc', control=control)
   if (result%status /= perm_ok) error stop 'K-sample example failed'
   print '(a,f10.6,a,2f10.6)', 'K-sample MC p: ', result%p_value, '  99% interval: ', result%confidence_interval
   call perm_trend(scores, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], result, method='pclt')
   if (result%status /= perm_ok) error stop 'Trend example failed'
   print '(a,f10.6,a,f10.6)', 'Trend asymptotic p: ', result%p_value, '  correlation: ', result%estimate
   control%set_seed = .false.
   rng%state = 271828_int64
   call twosample_exact_mc(scores, [1, 1, 1, 0, 0, 0], result, control, 'less', rng)
   if (result%status /= perm_ok) error stop 'Continued RNG example failed'
   print '(a,f10.6)', 'Explicit-stream one-sided MC p: ', result%p_value
end program permutation_examples
