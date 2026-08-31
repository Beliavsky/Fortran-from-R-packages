program test_profile
   use ordinal, only : dp, clm_problem, clm_fit_result, init_clm_problem, fit_clm, clm_profile_likelihood, &
                       clm_profile_confidence_interval, clm_wald_confidence_interval, link_logit, threshold_flexible
   implicit none
   integer, parameter :: n = 90
   integer :: y(n), status, i, beta_index
   integer :: profile_status(3)
   real(dp) :: x(n, 1), eta, u, lower, upper, wald_lower, wald_upper
   real(dp) :: values(3), profile_loglik(3)
   type(clm_problem) :: problem
   type(clm_fit_result) :: fit

   do i = 1, n
      x(i, 1) = -1.8_dp + 3.6_dp*real(i - 1, dp)/real(n - 1, dp)
      eta = 0.85_dp*x(i, 1)
      u = real(mod(37*i + 11, 97), dp)/97.0_dp
      if (u < 1.0_dp/(1.0_dp + exp(-(-0.7_dp - eta)))) then
         y(i) = 1
      else if (u < 1.0_dp/(1.0_dp + exp(-(0.8_dp - eta)))) then
         y(i) = 2
      else
         y(i) = 3
      end if
   end do

   call init_clm_problem(problem, y, x, 3, link=link_logit, threshold=threshold_flexible, status=status)
   if (status /= 0) error stop 'test_profile: initialization failed'
   call fit_clm(problem, fit, max_iter=250, grad_tol=1.0e-8_dp)
   if (.not. fit%converged) error stop 'test_profile: CLM fit did not converge'

   beta_index = size(fit%alpha) + 1
   values = [fit%beta(1) - 0.20_dp, fit%beta(1), fit%beta(1) + 0.20_dp]
   call clm_profile_likelihood(problem, fit, beta_index, values, profile_loglik, profile_status, max_iter=250, &
                               grad_tol=1.0e-8_dp)
   if (any(profile_status /= 0)) error stop 'test_profile: profile optimization failed'
   if (abs(profile_loglik(2) - fit%loglik) > 2.0e-5_dp) error stop 'test_profile: profile peak does not match MLE'
   if (profile_loglik(1) > profile_loglik(2) + 1.0e-7_dp) error stop 'test_profile: lower profile exceeds MLE'
   if (profile_loglik(3) > profile_loglik(2) + 1.0e-7_dp) error stop 'test_profile: upper profile exceeds MLE'

   call clm_wald_confidence_interval(fit, beta_index, 0.90_dp, wald_lower, wald_upper, status)
   if (status /= 0) error stop 'test_profile: Wald interval failed'
   if (wald_lower >= fit%beta(1) .or. wald_upper <= fit%beta(1)) error stop 'test_profile: Wald interval misses MLE'

   call clm_profile_confidence_interval(problem, fit, beta_index, 0.90_dp, lower, upper, status, max_iter=250, &
                                        grad_tol=1.0e-7_dp)
   if (status /= 0) error stop 'test_profile: profile interval failed'
   if (lower >= fit%beta(1) .or. upper <= fit%beta(1)) error stop 'test_profile: profile interval misses MLE'
   if (lower >= upper) error stop 'test_profile: profile interval is reversed'
   if (abs(lower - wald_lower) > 0.9_dp .or. abs(upper - wald_upper) > 0.9_dp) then
      error stop 'test_profile: profile and Wald intervals disagree implausibly'
   end if

   print *, 'test_profile: PASS'
end program test_profile
