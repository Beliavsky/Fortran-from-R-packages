program random_example
   use skellam, only : dp, i8, rskellam, seed_random_number
   implicit none
   integer(i8), allocatable :: sample(:)
   real(dp) :: mean_value, variance_value

   call seed_random_number(12345)
   sample = rskellam(10000, 8.5_dp, 10.25_dp)
   mean_value = sum(real(sample, dp))/real(size(sample), dp)
   variance_value = sum((real(sample, dp) - mean_value)**2)/real(size(sample) - 1, dp)
   print '(a,f10.5)', 'sample mean:     ', mean_value
   print '(a,f10.5)', 'sample variance: ', variance_value
end program random_example
