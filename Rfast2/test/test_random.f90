program test_random
   use rfast2
   use iso_fortran_env, only : int64
   implicit none
   real(dp), allocatable :: a(:),b(:),x(:)
   integer, allocatable :: s(:)
   integer :: i

   call set_seed(12345_int64)
   a = runif(12)
   call set_seed(12345_int64)
   b = runif(12)
   if (maxval(abs(a-b)) > 0.0_dp) error stop 1
   if (minval(a) < 0.0_dp .or. maxval(a) > 1.0_dp) error stop 2
   s = sample_int(20,20,replace=.false.)
   if (minval(s) < 1 .or. maxval(s) > 20) error stop 3
   do i=1,20
      if (count(s==i) /= 1) error stop 4
   end do
   call set_seed(99_int64)
   x = rgamma_fast(50000,2.5_dp,1.2_dp)
   if (abs(sum(x)/real(size(x),dp)-2.5_dp/1.2_dp) > 0.08_dp) error stop 5
   x = rcauchy_fast(10000,0.0_dp,1.0_dp)
   if (count(x < 0.0_dp) < 4500 .or. count(x < 0.0_dp) > 5500) error stop 6
   print '(a)', 'test_random: PASS'
end program test_random
