program test_rng
   use skewunit
   implicit none
   integer, parameter :: n = 60000
   real(dp) :: x(n), meanv, varv

   call seed_skewunit_rng(123456)
   call rskewunit_vec(n,x,lambda=0.0_dp,family1=family_triang,family2=family_none)
   meanv = sum(x)/real(n,dp)
   varv = sum((x-meanv)**2)/real(n-1,dp)
   if (abs(meanv-0.5_dp) > 0.006_dp) error stop 'triangular RNG mean'
   if (abs(varv-1.0_dp/24.0_dp) > 0.003_dp) error stop 'triangular RNG variance'

   call seed_skewunit_rng(13579)
   call rskewunit_vec(n,x,lambda=-0.4_dp,delta=1.4_dp, &
      family1=family_triang,family2=family_jsb)
   meanv = sum(x)/real(n,dp)
   if (meanv < 0.38_dp .or. meanv > 0.49_dp) error stop 'skew RNG direction'

   print '(a)', 'test_rng: PASS'
end program test_rng
