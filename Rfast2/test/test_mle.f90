program test_mle
   use rfast2
   implicit none
   type(mle_result) :: r
   real(dp) :: x(6)

   x = [1.0_dp,1.1_dp,1.4_dp,2.0_dp,3.0_dp,5.0_dp]
   r = powerlaw_mle(x)
   if (r%status /= 0 .or. size(r%param) /= 2) error stop 1
   if (abs(r%param(2)-1.0_dp) > 1.0e-12_dp) error stop 2
   if (r%param(1) <= 1.0_dp) error stop 3
   r = cauchy0_mle([-2.0_dp,-1.0_dp,-0.4_dp,0.3_dp,1.2_dp,2.1_dp])
   if (r%status /= 0 .or. r%param(1) <= 0.0_dp) error stop 4
   r = sp_mle([0.1_dp,0.2_dp,0.4_dp,0.7_dp])
   if (abs(r%param(1)+4.0_dp/sum(log([0.1_dp,0.2_dp,0.4_dp,0.7_dp]))) > 1.0e-12_dp) error stop 5
   print '(a)', 'test_mle: PASS'
end program test_mle
