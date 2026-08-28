program test_missing_ordinal
   use lavaan
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   implicit none
   real(dp) :: x(3,2), mu(2), sigma(2,2), ll1, ll2, rho, ll
   integer :: info
   integer :: tab(2,2)
   real(dp), allocatable :: t1(:), t2(:)

   x=reshape([0.2_dp,-0.1_dp,0.4_dp, 0.3_dp,0.5_dp,-0.2_dp],[3,2])
   mu=[0.0_dp,0.0_dp]
   sigma=reshape([1.0_dp,0.3_dp,0.3_dp,1.0_dp],[2,2])
   ll1=mvn_loglik_complete(x,mu,sigma,info)
   call check(info==0,'complete likelihood')
   ll2=mvn_loglik_missing(x,mu,sigma,info)
   call check(abs(ll1-ll2)<1e-10_dp,'missing equals complete')
   x(2,2)=ieee_value(0.0_dp,ieee_quiet_nan)
   ll2=mvn_loglik_missing(x,mu,sigma,info)
   call check(info==0 .and. ll2<0,'missing pattern likelihood')

   tab=reshape([200,100,100,200],[2,2])
   call polychoric_table(tab,rho,t1,t2,ll)
   call check(abs(t1(1))<2e-8_dp .and. abs(t2(1))<2e-8_dp,'binary thresholds')
   call check(abs(rho-0.5_dp)<2e-5_dp,'polychoric known rho')
   print '(a)', 'test_missing_ordinal: PASS'
contains
   subroutine check(cond,label)
      logical,intent(in)::cond
      character(len=*),intent(in)::label
      if(.not.cond) then
      write(*,'(a,1x,a)') 'FAIL:',trim(label)
      error stop 1
      end if
   end subroutine check
end program test_missing_ordinal
