program test_clmm
   use ordinal, only : dp, clmm_problem, clmm_fit_result, init_clmm_problem, fit_clmm, &
      link_logit, threshold_flexible
   implicit none
   integer, parameter :: ng = 12, per_group = 12, n = ng*per_group, k = 3
   real(dp) :: x(n, 1), eta, u, p1, p2, randeff(ng)
   integer :: y(n), group(n), g, j, i, status
   type(clmm_problem) :: problem
   type(clmm_fit_result) :: fit

   do g = 1, ng
      randeff(g) = 0.90_dp*sin(1.7_dp*real(g, dp))
      do j = 1, per_group
         i = (g - 1)*per_group + j
         group(i) = 100 + 3*g
         x(i, 1) = -1.5_dp + 3.0_dp*real(j - 1, dp)/real(per_group - 1, dp)
         eta = 0.8_dp*x(i, 1) + randeff(g)
         p1 = logistic(-0.5_dp - eta)
         p2 = logistic(0.8_dp - eta)
         u = real(mod(71*i + 19*g, 1009), dp)/1009.0_dp
         if (u < p1) then
            y(i) = 1
         else if (u < p2) then
            y(i) = 2
         else
            y(i) = 3
         end if
      end do
   end do
   call init_clmm_problem(problem, y, x, group, k, link=link_logit, threshold=threshold_flexible, status=status)
   if (status /= 0) error stop 'init_clmm_problem failed'
   call fit_clmm(problem, fit, max_iter=180, grad_tol=2.0e-5_dp)
   if (.not. fit%converged) then
      print *, 'clmm status=', fit%status, ' gradient=', maxval(abs(fit%gradient))
      error stop 'fit_clmm did not converge'
   end if
   if (fit%tau <= 0.0_dp .or. fit%tau > 3.0_dp) then
      print *, 'tau=', fit%tau
      error stop 'unreasonable random effect scale'
   end if
   if (size(fit%ranef_mode) /= ng) error stop 'random effect mode count'
   if (fit%theta(1) >= fit%theta(2)) error stop 'mixed threshold order'
   print *, 'test_clmm: PASS'
contains
   pure elemental real(dp) function logistic(z) result(p)
      real(dp), intent(in) :: z !! Scalar predictor transformed by the logistic CDF.
      if (z >= 0.0_dp) then
         p = 1.0_dp/(1.0_dp + exp(-z))
      else
         p = exp(z)/(1.0_dp + exp(z))
      end if
   end function logistic
end program test_clmm
