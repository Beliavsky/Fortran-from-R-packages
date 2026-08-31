program clm_example
   use ordinal, only : dp, clm_problem, clm_fit_result, init_clm_problem, fit_clm, link_logit, threshold_flexible
   implicit none
   integer, parameter :: n = 80
   integer :: y(n), status, i
   real(dp) :: x(n, 1), eta, p1, p2, u
   type(clm_problem) :: problem
   type(clm_fit_result) :: fit

   do i = 1, n
      x(i, 1) = -2.0_dp + 4.0_dp*real(i - 1, dp)/real(n - 1, dp)
      eta = 0.9_dp*x(i, 1)
      p1 = logistic(-0.7_dp - eta)
      p2 = logistic(0.8_dp - eta)
      u = real(mod(43*i, 251), dp)/251.0_dp
      if (u < p1) then
         y(i) = 1
      else if (u < p2) then
         y(i) = 2
      else
         y(i) = 3
      end if
   end do
   call init_clm_problem(problem, y, x, 3, link=link_logit, threshold=threshold_flexible, status=status)
   if (status /= 0) error stop 'failed to initialize example model'
   call fit_clm(problem, fit, grad_tol=2.0e-6_dp)
   if (.not. fit%converged) error stop 'example fit failed to converge'
   print '(a,2f12.6)', 'thresholds: ', fit%theta
   print '(a,f12.6)', 'beta:       ', fit%beta(1)
   print '(a,f12.6)', 'logLik:     ', fit%loglik
contains
   pure elemental real(dp) function logistic(z) result(p)
      real(dp), intent(in) :: z !! Scalar predictor transformed by the logistic CDF.
      if (z >= 0.0_dp) then
         p = 1.0_dp/(1.0_dp + exp(-z))
      else
         p = exp(z)/(1.0_dp + exp(z))
      end if
   end function logistic
end program clm_example
