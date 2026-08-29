program test_mackinnon
   use urca_kinds, only : dp
   use urca_mackinnon
   implicit none
   integer :: info
   real(dp) :: q, p, tab(6,8)

   q=mackinnon_quantile(0.05_dp,0,MACK_NC,MACK_TAU,info)
   call assert_zero(info,'MacKinnon NC asymptotic')
   call assert_close(q,-1.9408467713813484_dp,2.0e-10_dp,'MacKinnon NC q')
   p=mackinnon_pvalue(q,0,MACK_NC,MACK_TAU,info)
   call assert_prob(p,'MacKinnon NC inverse')

   q=mackinnon_quantile(0.05_dp,25,MACK_C,MACK_TAU,info)
   call assert_zero(info,'MacKinnon C N25')
   call assert_close(q,-2.9862201645679702_dp,2.0e-10_dp,'MacKinnon C N25 q')
   p=mackinnon_pvalue(q,25,MACK_C,MACK_TAU,info)
   call assert_prob(p,'MacKinnon C N25 inverse')

   q=mackinnon_quantile(0.05_dp,100,MACK_CT,MACK_NORM,info)
   call assert_zero(info,'MacKinnon CT norm N100')
   call assert_close(q,-20.471209982449302_dp,2.0e-8_dp,'MacKinnon CT norm q')
   p=mackinnon_pvalue(q,100,MACK_CT,MACK_NORM,info)
   call assert_prob(p,'MacKinnon CT norm inverse')

   q=mackinnon_quantile(0.05_dp,500,MACK_CTT,MACK_TAU,info)
   call assert_zero(info,'MacKinnon CTT N500')
   call assert_close(q,-3.8440857220791051_dp,2.0e-10_dp,'MacKinnon CTT q')
   p=mackinnon_pvalue(q,500,MACK_CTT,MACK_TAU,info)
   call assert_prob(p,'MacKinnon CTT inverse')

   call unitroot_table(MACK_C,MACK_TAU,tab,info)
   call assert_zero(info,'unitroot table')
   call assert_true(all(tab==tab),'unitroot table finite')

   print '(a)', 'test_mackinnon: PASS'
contains
   subroutine assert_close(actual,expected,tol,name)
      real(dp), intent(in) :: actual, expected, tol
      character(len=*), intent(in) :: name
      if(abs(actual-expected)>tol)then
         write(*,'(a,2(1x,es24.16))') 'FAIL '//trim(name),actual,expected
         error stop 1
      end if
   end subroutine
   subroutine assert_prob(actual,name)
      real(dp), intent(in) :: actual
      character(len=*), intent(in) :: name
      if(abs(actual-0.05_dp)>5.0e-5_dp)then
         write(*,'(a,1x,es24.16)') 'FAIL '//trim(name),actual
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
end program test_mackinnon
