program test_clm_nominal
   use ordinal, only : dp, clm_problem, clm_fit_result, init_clm_problem, fit_clm, clm_predict_proba_nominal, &
                       clm_analytic_gradient, clm_analytic_hessian, link_logit, threshold_flexible
   use ordinal_numerics, only : numerical_gradient, numerical_hessian
   implicit none
   integer, parameter :: n = 72, k = 3
   integer :: y(n), i, status
   real(dp) :: x(n, 1), nom(n, 1), latent, theta1, theta2
   real(dp) :: par(5), ga(5), gn(5), ha(5, 5), hn(5, 5)
   real(dp) :: pnew(4, k), xnew(4, 1), nnew(4, 1)
   type(clm_problem) :: problem
   type(clm_fit_result) :: fit

   do i = 1, n
      x(i, 1) = (real(mod(i - 1, 9), dp) - 4.0_dp)/4.0_dp
      nom(i, 1) = real(mod((i - 1)/6, 2), dp)
      theta1 = -0.9_dp + 0.35_dp*nom(i, 1)
      theta2 = 0.95_dp - 0.15_dp*nom(i, 1)
      latent = 0.65_dp*x(i, 1) + 1.15_dp*sin(0.83_dp*real(i, dp))
      if (latent <= theta1) then
         y(i) = 1
      else if (latent <= theta2) then
         y(i) = 2
      else
         y(i) = 3
      end if
   end do

   call init_clm_problem(problem, y, x, k, link=link_logit, threshold=threshold_flexible, nominal_x=nom, status=status)
   if (status /= 0) error stop 'nominal problem initialization failed'

   par = [-0.8_dp, 0.9_dp, 0.20_dp, -0.10_dp, 0.45_dp]
   call clm_analytic_gradient(par, problem, ga, status)
   if (status /= 0) error stop 'analytic nominal gradient failed'
   call numerical_gradient(problem, par, gn, 2.0e-6_dp)
   if (maxval(abs(ga - gn)) > 2.0e-5_dp) error stop 'analytic nominal gradient mismatch'

   call clm_analytic_hessian(par, problem, ha, status)
   if (status /= 0) error stop 'analytic nominal Hessian failed'
   call numerical_hessian(problem, par, hn, 2.0e-4_dp)
   if (maxval(abs(ha - hn)) > 2.0e-3_dp) error stop 'analytic nominal Hessian mismatch'

   call fit_clm(problem, fit, max_iter=250, grad_tol=1.0e-7_dp)
   if (.not. fit%converged) error stop 'nominal CLM did not converge'
   if (size(fit%alpha) /= 4 .or. any(shape(fit%nominal_alpha) /= [2, 2])) error stop 'nominal coefficient shape is wrong'
   if (fit%max_gradient > 2.0e-5_dp) error stop 'nominal CLM final gradient is too large'

   xnew(:, 1) = [-0.75_dp, -0.25_dp, 0.25_dp, 0.75_dp]
   nnew(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp]
   call clm_predict_proba_nominal(fit%alpha, fit%beta, xnew, nnew, k, threshold_flexible, link_logit, pnew, &
                                  status=status)
   if (status /= 0) error stop 'nominal prediction failed'
   do i = 1, size(pnew, 1)
      if (abs(sum(pnew(i, :)) - 1.0_dp) > 1.0e-12_dp) error stop 'nominal predicted probabilities do not sum to one'
      if (any(pnew(i, :) <= 0.0_dp)) error stop 'nominal predicted probability is not positive'
   end do

   print *, 'test_clm_nominal: PASS'
end program test_clm_nominal
