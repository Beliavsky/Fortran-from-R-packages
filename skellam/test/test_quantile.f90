program test_quantile
   use skellam, only : dp, i8, qskellam, pskellam
   implicit none

   real(dp), parameter :: probabilities(5) = [0.001_dp, 0.05_dp, 0.5_dp, 0.95_dp, 0.999_dp]
   integer(i8), parameter :: expected(5) = [-6_i8, -3_i8, 1_i8, 5_i8, 8_i8]
   integer(i8) :: quantiles(5), q
   real(dp) :: p, cdf_at, cdf_before
   integer :: i

   quantiles = qskellam(probabilities, 3.0_dp, 2.0_dp)
   if (any(quantiles /= expected)) then
      print *, 'quantile reference failed:', quantiles
      error stop 1
   end if
   do i = 1, size(probabilities)
      p = probabilities(i)
      q = quantiles(i)
      cdf_at = pskellam(real(q, dp), 3.0_dp, 2.0_dp)
      cdf_before = pskellam(real(q - 1_i8, dp), 3.0_dp, 2.0_dp)
      if (cdf_at < p - 1.0e-14_dp .or. cdf_before >= p) error stop 'quantile inversion failed'
   end do

   q = qskellam(log(0.05_dp), 3.0_dp, 2.0_dp, log_p=.true.)
   if (q /= -3_i8) error stop 'log probability quantile failed'
   q = qskellam(0.05_dp, 3.0_dp, 2.0_dp, lower_tail=.false.)
   if (q /= 5_i8) error stop 'upper-tail quantile failed'
   if (qskellam(0.0_dp, 5.0_dp, 0.0_dp) /= 0_i8) error stop 'bounded lower support failed'
   print '(a)', 'test_quantile: PASS'
end program test_quantile
