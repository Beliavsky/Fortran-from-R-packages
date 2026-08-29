program test_nonnested
   use lmtest, only : dp, nonnested_test_result, pair_test_result, &
      cox_test, j_test, encompassing_test, pe_test
   implicit none
   integer, parameter :: n=80
   real(dp) :: x(n,3), z(n,3), y(n), yraw(n), ylog(n), v
   integer :: i
   type(nonnested_test_result) :: nt
   type(pair_test_result) :: pt

   do i=1,n
      v=real(i,dp)
      x(i,:)=[1.0_dp,(v-40.0_dp)/20.0_dp,sin(0.3_dp*v)]
      z(i,:)=[1.0_dp,(v-40.0_dp)/20.0_dp,cos(0.21_dp*v)]
      y(i)=1.0_dp+2.0_dp*x(i,2)-0.7_dp*x(i,3)+ &
         (0.18_dp+0.002_dp*v)*sin(0.91_dp*v)+0.08_dp*cos(0.17_dp*v)
   end do

   nt=cox_test(y,x,z)
   call check(nt%estimate(1),-0.010875939491692826_dp,2.0e-11_dp)
   call check(nt%statistic(1),-0.04058185983851577_dp,2.0e-10_dp)

   nt=j_test(y,x,z)
   call check(nt%estimate(1),0.020485624254740577_dp,2.0e-10_dp)
   call check(nt%p_value(1),0.9688789176522763_dp,2.0e-10_dp)

   pt=encompassing_test(y,x,z)
   call check(pt%first%statistic,0.0015321845576856952_dp,2.0e-10_dp)
   call check(pt%first%p_value,0.9688789176523341_dp,2.0e-10_dp)
   call check(pt%second%statistic,468.4879937954972_dp,2.0e-7_dp)

   do i=1,n
      v=real(i,dp)
      yraw(i)=exp(1.0_dp+0.2_dp*x(i,2)-0.1_dp*x(i,3)+ &
         0.04_dp*sin(0.7_dp*v)+0.02_dp*cos(0.13_dp*v))
      ylog(i)=log(yraw(i))
   end do
   nt=pe_test(ylog,x,yraw,z,.true.,.false.)
   call check(nt%estimate(1),0.02383245403879549_dp,2.0e-10_dp)
   call check(nt%p_value(1),0.5488048241703208_dp,2.0e-9_dp)
contains
   subroutine check(actual,expected,tol)
      real(dp),intent(in)::actual,expected,tol
      if(abs(actual-expected)>tol)then
         print *,actual,expected
         error stop 1
      end if
   end subroutine check
end program test_nonnested
