program test_breaks
   use urca_kinds, only : dp
   use urca_types, only : johansen_result
   use urca_breaks, only : johansen_level_shift
   implicit none
   integer, parameter :: n = 240
   real(dp) :: d(n,4), x(n,3)
   integer :: i, iu
   type(johansen_result) :: r

   open(newunit=iu,file='test/reference_data.txt',status='old',action='read')
   do i=1,n
      read(iu,*) d(i,:)
   end do
   close(iu)
   x=d(:,2:4)

   r=johansen_level_shift(x,.false.,2)
   call assert_zero(r%info,'cajolst status')
   call assert_int(r%break_point,192,'cajolst break')
   call assert_close(r%lambda(1),0.33166931555173257_dp,1.0e-10_dp,'cajolst lambda1')
   call assert_close(r%lambda(2),0.21373211440875084_dp,1.0e-10_dp,'cajolst lambda2')
   call assert_close(r%lambda(3),0.032935068710980744_dp,1.0e-10_dp,'cajolst lambda3')
   call assert_close(r%teststat(1),7.712230814688656_dp,2.0e-8_dp,'cajolst trace1')
   call assert_close(r%teststat(2),53.81283193407263_dp,2.0e-8_dp,'cajolst trace2')
   call assert_close(r%teststat(3),121.9839525022844_dp,3.0e-8_dp,'cajolst trace3')

   print '(a)', 'test_breaks: PASS'
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
   subroutine assert_int(actual,expected,name)
      integer, intent(in) :: actual, expected
      character(len=*), intent(in) :: name
      if(actual/=expected)then
         write(*,'(a,2(1x,i0))') 'FAIL '//trim(name),actual,expected
         error stop 1
      end if
   end subroutine
end program test_breaks
