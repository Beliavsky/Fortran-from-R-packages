program test_api_tails
   use r_compat, only: dp, set_seed_int
   use stabledist
   implicit none
   real(dp) :: p,pu,d,ld,lp,lpu,x
   real(dp), allocatable :: dv(:),pv(:),rv(:)
   real(dp), parameter :: xv(3)=[-1.0_dp,0.0_dp,1.0_dp]
   real(dp), parameter :: gv(2)=[1.0_dp,2.0_dp]
   real(dp), parameter :: de(2)=[0.0_dp,0.5_dp]

   x=1.3_dp
   p=pstable(x,1.4_dp,-0.3_dp)
   pu=pstable(x,1.4_dp,-0.3_dp,lower_tail=.false.)
   if(abs(p+pu-1.0_dp)>3e-10_dp)error stop 'lower+upper != 1'
   lp=pstable(x,1.4_dp,-0.3_dp,log_p=.true.)
   lpu=pstable(x,1.4_dp,-0.3_dp,lower_tail=.false.,log_p=.true.)
   if(abs(exp(lp)-p)>3e-10_dp .or. abs(exp(lpu)-pu)>3e-10_dp)error stop 'log probability mismatch'
   d=dstable(x,1.4_dp,-0.3_dp)
   ld=dstable(x,1.4_dp,-0.3_dp,log_=.true.)
   if(abs(exp(ld)-d)>3e-10_dp)error stop 'log density mismatch'

   dv=dstable(xv,1.5_dp,0.2_dp,gamma=gv,delta=de)
   pv=pstable(xv,1.5_dp,0.2_dp,gamma=gv,delta=de)
   if(size(dv)/=3.or.size(pv)/=3)error stop 'vector API size'
   if(any(dv<=0.0_dp).or.any(pv<0.0_dp).or.any(pv>1.0_dp))error stop 'vector API range'

   call set_seed_int(111)
   rv=rstable_varying(8,1.5_dp,0.2_dp,gv,de)
   if(size(rv)/=8)error stop 'rstable_varying size'
   if(any(rv/=rv))error stop 'rstable_varying NaN'

   print *, 'test_api_tails: PASS'
end program test_api_tails
