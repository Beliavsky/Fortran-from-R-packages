! SPDX-License-Identifier: GPL-2.0-or-later
program dip_diagnostics
   use diptest_api, only : dp, dip_compute, dip_result, dip_ok
   implicit none
   type(dip_result) :: result
   real(dp), parameter :: x(12) = [198.0_dp, 198.0_dp, 1961.0_dp, 2898.0_dp, 3184.0_dp, 3687.0_dp, &
      4336.0_dp, 4987.0_dp, 5661.0_dp, 6530.0_dp, 7476.0_dp, 8555.0_dp]
   call dip_compute(x, result)
   if (result%status /= dip_ok) error stop 'diagnostic calculation failed'
   print '(a,f12.8)', 'Dip: ', result%statistic
   print '(a,2i6)', 'Modal indices: ', result%low, result%high
   print '(a,2f12.3)', 'Modal interval: ', result%lower, result%upper
   print '(a,*(i0,1x))', 'Modal GCM, descending: ', result%gcm
   print '(a,*(i0,1x))', 'Modal LCM, ascending: ', result%lcm
   print '(a,*(i0,1x))', 'Full GCM, descending: ', result%full_gcm
   print '(a,*(i0,1x))', 'Full LCM, ascending: ', result%full_lcm
end program dip_diagnostics
