program demo_skellam
   use skellam, only : dp, i8, dskellam, pskellam, qskellam, rskellam, &
      seed_random_number, skellam_mle_result, fit_skellam_mle
   implicit none
   integer(i8), allocatable :: sample(:)
   type(skellam_mle_result) :: fit

   print '(a,f12.8)', 'P(X = 1), X ~ Skellam(3,2): ', dskellam(1_i8, 3.0_dp, 2.0_dp)
   print '(a,f12.8)', 'P(X <= 1):                 ', pskellam(1.0_dp, 3.0_dp, 2.0_dp)
   print '(a,i0)', '95% quantile:                ', qskellam(0.95_dp, 3.0_dp, 2.0_dp)
   call seed_random_number(20260801)
   sample = rskellam(1000, 10.0_dp, 6.0_dp)
   call fit_skellam_mle(sample, fit)
   print '(a,2f11.5)', 'fitted rates: ', fit%lambda1, fit%lambda2
end program demo_skellam
