program test_arrays
   use rfast2
   implicit none
   real(dp) :: x(4),q,tm
   real(dp), allocatable :: z(:),m(:,:)
   integer :: g(4)
   real(dp) :: a(4,2),acc(2)
   integer :: pred(4,2)

   x = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
   q = quantile_rfast2(x,0.25_dp)
   if (abs(q-1.75_dp) > 1.0e-12_dp) error stop 1
   tm = trim_mean([1.0_dp,2.0_dp,3.0_dp,4.0_dp,100.0_dp],0.2_dp)
   if (abs(tm-3.0_dp) > 1.0e-12_dp) error stop 2
   z = intersect_real([1.0_dp,2.0_dp,2.0_dp,3.0_dp],[2.0_dp,4.0_dp])
   if (size(z) /= 1 .or. abs(z(1)-2.0_dp) > 1.0e-12_dp) error stop 3
   a = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,10.0_dp,20.0_dp,30.0_dp,40.0_dp],[4,2])
   g = [1,1,2,2]
   m = col_group_sum(a,g,2)
   if (maxval(abs(m-reshape([3.0_dp,7.0_dp,30.0_dp,70.0_dp],[2,2]))) > 1.0e-12_dp) error stop 4
   pred = reshape([1,0,1,0,1,1,0,0],[4,2])
   acc = col_accuracy([1,0,1,0],pred)
   if (abs(acc(1)-1.0_dp) > 1.0e-12_dp) error stop 5
   print '(a)', 'test_arrays: PASS'
end program test_arrays
