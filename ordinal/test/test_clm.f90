program test_clm
   use ordinal, only : dp, clm_problem, clm_fit_result, init_clm_problem, fit_clm, clm_predict_proba, &
      link_logit, threshold_flexible
   implicit none
   integer, parameter :: n = 360, k = 3
   real(dp) :: x(n, 1), probs(n, k), u, eta, p1, p2
   integer :: y(n), i, status
   type(clm_problem) :: problem
   type(clm_fit_result) :: fit

   do i = 1, n
      x(i, 1) = -2.0_dp + 4.0_dp*real(i - 1, dp)/real(n - 1, dp)
      eta = 1.1_dp*x(i, 1)
      p1 = logistic(-0.6_dp - eta)
      p2 = logistic(0.9_dp - eta)
      u = real(mod(37*i, 997), dp)/997.0_dp
      if (u < p1) then
         y(i) = 1
      else if (u < p2) then
         y(i) = 2
      else
         y(i) = 3
      end if
   end do
   call init_clm_problem(problem, y, x, k, link=link_logit, threshold=threshold_flexible, status=status)
   if (status /= 0) error stop 'init_clm_problem failed'
   call fit_clm(problem, fit, max_iter=220, grad_tol=2.0e-6_dp)
   if (.not. fit%converged) then
      print *, 'clm status=', fit%status, ' gradient=', maxval(abs(fit%gradient))
      error stop 'fit_clm did not converge'
   end if
   if (fit%theta(1) >= fit%theta(2)) error stop 'threshold order'
   if (abs(fit%beta(1) - 1.1_dp) > 0.35_dp) then
      print *, 'beta=', fit%beta(1)
      error stop 'clm beta parity tolerance'
   end if
   if (abs(fit%theta(1) + 0.6_dp) > 0.35_dp .or. abs(fit%theta(2) - 0.9_dp) > 0.35_dp) then
      print *, 'theta=', fit%theta
      error stop 'clm threshold parity tolerance'
   end if
   call clm_predict_proba(fit%theta, fit%beta, x, link_logit, probs, status=status)
   if (status /= 0) error stop 'clm_predict_proba failed'
   do i = 1, n
      if (abs(sum(probs(i, :)) - 1.0_dp) > 2.0e-13_dp) error stop 'probabilities do not sum to one'
      if (any(probs(i, :) < 0.0_dp)) error stop 'negative probability'
   end do
   if (maxval(abs(fit%gradient)) > 2.0e-4_dp) error stop 'final gradient too large'
   print *, 'test_clm: PASS'
contains
   pure elemental real(dp) function logistic(z) result(p)
      real(dp), intent(in) :: z !! Scalar predictor transformed by the logistic CDF.
      if (z >= 0.0_dp) then
         p = 1.0_dp/(1.0_dp + exp(-z))
      else
         p = exp(z)/(1.0_dp + exp(z))
      end if
   end function logistic
end program test_clm
