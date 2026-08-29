program test_interpolation
use tweedie
implicit none
real(dp),parameter::pv(10)=[1.15_dp,1.25_dp,1.35_dp,1.45_dp,1.75_dp, &
   2.5_dp,3.5_dp,4.5_dp,6.0_dp,8.5_dp]
real(dp),parameter::xh(10)=[0.1_dp,0.3_dp,0.5_dp,0.8_dp,0.9_dp, &
   0.9_dp,0.9_dp,0.9_dp,0.5_dp,0.3_dp]
real(dp)::xix,phi,interp_d,invert_d,rel
integer::i,failures
failures=0
do i=1,size(pv)
   xix=0.5_dp*xh(i)
   phi=xix/(1.0_dp-xix)
   interp_d=dtweedie(1.0_dp,1.0_dp,phi,pv(i))
   invert_d=dtweedie_inversion(1.0_dp,1.0_dp,phi,pv(i))
   rel=abs(interp_d-invert_d)/max(abs(invert_d),tiny(1.0_dp))
   if(rel>1.0e-6_dp)then
      print '(a,f6.2,a,es12.4)','grid p=',pv(i),' relative error=',rel
      failures=failures+1
   end if
end do
if(failures/=0)error stop 'test_interpolation failed'
print '(a)','test_interpolation: PASS'
end program test_interpolation
