program test_granger
   use lmtest, only : dp, test_result, granger_test
   implicit none
   integer, parameter :: n=120
   real(dp) :: x(n), y(n), v
   integer :: i
   type(test_result) :: tr
   x=0.0_dp
   y=0.0_dp
   do i=2,n
      v=real(i,dp)
      x(i)=0.55_dp*x(i-1)+0.3_dp*sin(0.37_dp*v)+0.12_dp*cos(0.11_dp*v)
      y(i)=0.45_dp*y(i-1)+0.75_dp*x(i-1)+0.15_dp*sin(0.83_dp*v)
   end do
   tr=granger_test(x,y,order=2)
   call check(tr%statistic,107.97877053759603_dp,3.0e-9_dp)
   call check(tr%p_value,6.032963797298744e-27_dp,2.0e-37_dp)
   if(abs(tr%df1-2.0_dp)>1.0e-14_dp .or. abs(tr%df2-113.0_dp)>1.0e-14_dp) error stop 1
contains
   subroutine check(actual,expected,tol)
      real(dp),intent(in)::actual,expected,tol
      if(abs(actual-expected)>tol)then
         print *,actual,expected
         error stop 1
      end if
   end subroutine check
end program test_granger
