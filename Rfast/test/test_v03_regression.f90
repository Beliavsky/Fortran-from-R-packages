program test_v03_regression
   use rfast
   implicit none
   integer, parameter :: n1=100, n2=120, n3=80
   real(dp) :: x1(n1,1), y1(n1), x2(n2,1), y2(n2), y3(n3), z, u
   integer :: i, failures
   integer :: yo(10)
   type(regression_result) :: nl
   type(weibull_regression_result) :: wr
   type(tobit_result) :: tr
   type(ordinal_result) :: ol

   failures=0
   do i=1,n1
      x1(i,1)=sin(real(i,dp)*0.37_dp)
      y1(i)=exp(0.3_dp+1.1_dp*x1(i,1))+0.05_dp*cos(real(i,dp)*1.23_dp)+0.08_dp
   end do
   nl=normlog_regression(y1,x1)
   call check(nl%status==0,'normlog status')
   call check(abs(nl%beta(1)-0.3657352942_dp)<2.0e-7_dp,'normlog intercept')
   call check(abs(nl%beta(2)-1.0497434518_dp)<2.0e-7_dp,'normlog slope')
   call check(abs(nl%deviance-0.1488329984_dp)<2.0e-8_dp,'normlog objective')

   do i=1,n2
      x2(i,1)=sin(real(i,dp)*0.43_dp)
      u=modulo(abs(sin(real(i,dp)*12.9898_dp)*43758.5453_dp),1.0_dp)
      u=min(1.0_dp-1.0e-8_dp,max(1.0e-8_dp,u))
      y2(i)=exp(0.2_dp+0.7_dp*x2(i,1))*(-log(u))**(1.0_dp/1.7_dp)
   end do
   wr=weibull_regression(y2,x2)
   call check(wr%status==0,'weibull regression status')
   call check(abs(wr%beta(1)-0.1081103520_dp)<2.0e-6_dp,'weibull intercept')
   call check(abs(wr%beta(2)-0.8135535763_dp)<2.0e-6_dp,'weibull slope')
   call check(abs(wr%shape-1.6018584101_dp)<2.0e-6_dp,'weibull shape')
   call check(abs(wr%loglik+104.17730631_dp)<2.0e-5_dp,'weibull loglik')

   do i=1,n3
      z=0.35_dp+0.65_dp*sin(real(i,dp)*0.73_dp)+0.25_dp*cos(real(i,dp)*1.31_dp)
      y3(i)=max(0.0_dp,z)
   end do
   tr=tobit_mle(y3)
   call check(tr%status==0,'tobit status')
   call check(abs(tr%location-0.3169655968_dp)<2.0e-7_dp,'tobit location')
   call check(abs(tr%scale-0.5607892654_dp)<2.0e-7_dp,'tobit scale')
   call check(abs(tr%loglik+69.200218474_dp)<2.0e-7_dp,'tobit loglik')

   yo=[1,1,1,2,2,2,2,3,3,3]
   ol=ordinal_mle(yo,ORDINAL_LOGIT)
   call check(ol%status==0,'ordinal status')
   call check(maxval(abs(ol%threshold-[log(0.3_dp/0.7_dp),log(0.7_dp/0.3_dp)]))<1.0e-12_dp,'ordinal thresholds')
   call check(abs(ol%loglik-(3.0_dp*log(0.3_dp)+4.0_dp*log(0.4_dp)+3.0_dp*log(0.3_dp)))<1.0e-12_dp,'ordinal loglik')

   if(failures==0)then
      print *,'test_v03_regression: PASS'
   else
      print *,'test_v03_regression: FAIL',failures
      error stop 1
   end if
contains
   subroutine check(ok,name)
      logical,intent(in)::ok
      character(*),intent(in)::name
      if(.not.ok)then
         print *,'FAIL: ',trim(name)
         failures=failures+1
      end if
   end subroutine check
end program test_v03_regression
