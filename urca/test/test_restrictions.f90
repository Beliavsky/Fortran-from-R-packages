program test_restrictions
   use urca_kinds, only : dp
   use urca_types, only : johansen_result, restriction_result
   use urca_cointegration, only : johansen_test, JO_TRACE, JO_NONE, JO_LONGRUN
   use urca_restrictions
   implicit none
   integer, parameter :: n = 240
   real(dp) :: d(n,4), x(n,3), h(3,2)
   integer :: i, iu
   type(johansen_result) :: jo
   type(restriction_result) :: r

   open(newunit=iu,file='test/reference_data.txt',status='old',action='read')
   do i=1,n
      read(iu,*) d(i,:)
   end do
   close(iu)
   x=d(:,2:4)
   jo=johansen_test(x,JO_TRACE,JO_NONE,2,JO_LONGRUN)
   call assert_zero(jo%info,'Johansen status')
   h=0.0_dp
   h(1,1)=1.0_dp
   h(2,2)=1.0_dp

   r=beta_restriction_test(jo,h,1)
   call assert_zero(r%info,'blr status')
   call assert_close(r%statistic,23.203235948811248_dp,1.0e-8_dp,'blr')

   r=alpha_restriction_test(jo,h,1)
   call assert_zero(r%info,'alr status')
   call assert_close(r%statistic,15.67438347881805_dp,1.0e-8_dp,'alr')

   r=alpha_beta_restriction_test(jo,h,h,1)
   call assert_zero(r%info,'ablr status')
   call assert_close(r%statistic,23.461224117377938_dp,1.0e-8_dp,'ablr')

   r=partly_known_beta_test(jo,h(:,1:1),2)
   call assert_zero(r%info,'bh5 status')
   call assert_close(r%statistic,56.64603846973627_dp,5.0e-8_dp,'bh5')

   r=iterated_partly_known_beta_test(jo,h,2,1)
   call assert_true(r%info>=0,'bh6 status')
   call assert_close(r%statistic,56.63732560007072_dp,5.0e-8_dp,'bh6')

   r=linear_trend_lr_test(x,2,1)
   call assert_zero(r%info,'lttest status')
   call assert_close(r%statistic,2.9542314288676184_dp,1.0e-8_dp,'lttest')

   print '(a)', 'test_restrictions: PASS'
contains
   subroutine assert_close(actual,expected,tol,name)
      real(dp), intent(in) :: actual, expected, tol
      character(len=*), intent(in) :: name
      if(abs(actual-expected)>tol)then
         write(*,'(a,2(1x,es24.16))') 'FAIL '//trim(name),actual,expected
         error stop 1
      end if
   end subroutine
   subroutine assert_zero(actual,name)
      integer, intent(in) :: actual
      character(len=*), intent(in) :: name
      if(actual/=0)then
         write(*,'(a,1x,i0)') 'FAIL '//trim(name),actual
         error stop 1
      end if
   end subroutine
   subroutine assert_true(ok,name)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: name
      if(.not.ok)then
         write(*,'(a)') 'FAIL '//trim(name)
         error stop 1
      end if
   end subroutine
end program test_restrictions
