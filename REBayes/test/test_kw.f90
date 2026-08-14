program test_kw
   use rebayes_kinds, only : dp
   use rebayes_kw
   implicit none
   real(dp) :: a(3,2),d(2),w(3)
   type(kw_result)::r
   type(kw_control)::c
   a=reshape([0.8_dp,0.2_dp,0.5_dp, 0.2_dp,0.8_dp,0.5_dp],[3,2])
   d=1.0_dp;w=[0.4_dp,0.4_dp,0.2_dp]
   c%tol=1.0e-10_dp;c%max_iter=10000
   call kw_fit(a,d,w,r,c)
   if(maxval(abs(r%f-[0.5_dp,0.5_dp]))>1.0e-7_dp)error stop "kw weights"
   if(r%kkt_gap>1.0e-7_dp)error stop "kw kkt"
   print *,"test_kw: PASS"
end program test_kw
