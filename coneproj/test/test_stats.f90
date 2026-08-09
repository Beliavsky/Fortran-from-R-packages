program test_stats
   use coneproj
   implicit none
   real(dp) :: p,q
   p=regularized_beta(0.5_dp,2.0_dp,2.0_dp)
   if (abs(p-0.5_dp)>1.0e-10_dp) error stop 'beta cdf'
   p=student_t_cdf(0.0_dp,10.0_dp)
   if (abs(p-0.5_dp)>1.0e-12_dp) error stop 't cdf'
   q=student_t_quantile(0.975_dp,10.0_dp)
   if (abs(q-2.2281388519649385_dp)>1.0e-7_dp) then
      print *, q
      error stop 't quantile'
   end if
   print *, 'test_stats: PASS'
end program test_stats
