program test_random
   use skellam, only : dp, i8, rskellam, seed_random_number
   implicit none

   integer(i8), allocatable :: sample(:)
   real(dp) :: mean_value, variance_value
   integer, parameter :: n = 50000

   call seed_random_number(1729)
   sample = rskellam(n, 8.5_dp, 10.25_dp)
   mean_value = sum(real(sample, dp))/real(n, dp)
   variance_value = sum((real(sample, dp) - mean_value)**2)/real(n - 1, dp)
   if (abs(mean_value + 1.75_dp) > 0.08_dp) then
      print *, 'random mean failed:', mean_value
      error stop 1
   end if
   if (abs(variance_value - 18.75_dp) > 0.30_dp) then
      print *, 'random variance failed:', variance_value
      error stop 1
   end if
   print '(a)', 'test_random: PASS'
end program test_random
