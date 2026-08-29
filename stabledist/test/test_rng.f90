program test_rng
   use r_compat, only: dp, set_seed_int
   use stabledist
   implicit none
   real(dp), allocatable :: x(:)
   real(dp) :: emp,ref
   integer, parameter :: n=50000

   call set_seed_int(12345)
   x=rstable(n,1.5_dp,0.4_dp,gamma=2.0_dp,delta=0.3_dp,pm=0)
   emp=real(count(x<=0.7_dp),dp)/real(n,dp)
   ref=pstable(0.7_dp,1.5_dp,0.4_dp,gamma=2.0_dp,delta=0.3_dp,pm=0)
   if(abs(emp-ref)>0.012_dp)then
      print *,'rng S0 alpha1.5',emp,ref
      error stop 1
   end if

   call set_seed_int(54321)
   x=rstable(n,1.0_dp,0.6_dp,gamma=2.3_dp,delta=-0.2_dp,pm=0)
   emp=real(count(x<=0.7_dp),dp)/real(n,dp)
   ref=pstable(0.7_dp,1.0_dp,0.6_dp,gamma=2.3_dp,delta=-0.2_dp,pm=0)
   if(abs(emp-ref)>0.012_dp)then
      print *,'rng S0 alpha1 skew',emp,ref
      error stop 1
   end if

   call set_seed_int(24680)
   x=rstable(n,1.0_dp,0.6_dp,gamma=2.3_dp,delta=-0.2_dp,pm=1)
   emp=real(count(x<=0.7_dp),dp)/real(n,dp)
   ref=pstable(0.7_dp,1.0_dp,0.6_dp,gamma=2.3_dp,delta=-0.2_dp,pm=1)
   if(abs(emp-ref)>0.012_dp)then
      print *,'rng S1 alpha1 skew',emp,ref
      error stop 1
   end if

   print *, 'test_rng: PASS'
end program test_rng
