program test_random
   use biasedurn
   implicit none
   integer, allocatable :: xf(:), xw(:), mf(:,:), mw(:,:)
   integer :: m(2), failures
   real(dp) :: odds(2), em, sample_mean

   failures = 0
   call biasedurn_seed(12345)
   xf = rfnchypergeo(12000, 10, 15, 8, 2.5_dp)
   em = meanfnchypergeo(10, 15, 8, 2.5_dp)
   sample_mean = sum(real(xf, dp)) / real(size(xf), dp)
   if (abs(sample_mean - em) > 0.06_dp) failures = failures + 1

   xw = rwnchypergeo(12000, 10, 15, 8, 2.5_dp)
   em = meanwnchypergeo(10, 15, 8, 2.5_dp)
   sample_mean = sum(real(xw, dp)) / real(size(xw), dp)
   if (abs(sample_mean - em) > 0.06_dp) failures = failures + 1

   m = [10, 15]
   odds = [2.5_dp, 1.0_dp]
   mf = rmfnchypergeo(6000, m, 8, odds)
   sample_mean = sum(real(mf(1,:), dp)) / real(size(mf, 2), dp)
   if (abs(sample_mean - meanfnchypergeo(10, 15, 8, 2.5_dp)) > 0.08_dp) &
      failures = failures + 1

   mw = rmwnchypergeo(6000, m, 8, odds)
   sample_mean = sum(real(mw(1,:), dp)) / real(size(mw, 2), dp)
   if (abs(sample_mean - meanwnchypergeo(10, 15, 8, 2.5_dp)) > 0.08_dp) &
      failures = failures + 1

   if (failures == 0) then
      print *, 'test_random: PASS'
   else
      print *, 'test_random: FAIL', failures
      error stop 1
   end if
end program test_random
