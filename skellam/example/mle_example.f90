program mle_example
   use skellam, only : dp, i8, rskellam, seed_random_number, skellam_mle_result, fit_skellam_mle
   implicit none
   integer(i8), allocatable :: sample(:)
   type(skellam_mle_result) :: fit

   call seed_random_number(2026)
   sample = rskellam(2000, 10.0_dp, 6.0_dp)
   call fit_skellam_mle(sample, fit)
   print '(a,l1)', 'converged: ', fit%converged
   print '(a,2f12.6)', 'lambda estimates: ', fit%lambda1, fit%lambda2
   print '(a,2f12.6)', 'standard errors:  ', fit%standard_error1, fit%standard_error2
   print '(a,f14.4)', 'log likelihood:   ', fit%log_likelihood
end program mle_example
