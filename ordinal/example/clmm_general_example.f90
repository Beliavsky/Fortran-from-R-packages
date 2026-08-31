program clmm_general_example
   use ordinal, only : dp, clmm_laplace_problem, clmm_laplace_fit_result, init_clmm_laplace_problem, &
                       fit_clmm_laplace, link_logit, threshold_flexible
   implicit none
   integer, parameter :: n_group = 8
   integer, parameter :: per_group = 5
   integer, parameter :: n = n_group*per_group
   integer :: y(n), group(n, 1), term_q(1), status, i, g
   real(dp) :: x(n, 1), re_z(n, 2), latent, u, start(6)
   type(clmm_laplace_problem) :: problem
   type(clmm_laplace_fit_result) :: fit

   term_q = [2]
   do i = 1, n
      g = (i - 1)/per_group + 1
      group(i, 1) = g
      x(i, 1) = -1.2_dp + 2.4_dp*real(mod(3*i, 17), dp)/16.0_dp
      re_z(i, :) = [1.0_dp, x(i, 1)]
      latent = 0.75_dp*x(i, 1) + 0.30_dp*sin(real(g, dp)) + 0.12_dp*cos(real(g, dp))*x(i, 1)
      u = real(mod(29*i + 7, 101), dp)/101.0_dp
      if (u < 1.0_dp/(1.0_dp + exp(-(-0.65_dp - latent)))) then
         y(i) = 1
      else if (u < 1.0_dp/(1.0_dp + exp(-(0.85_dp - latent)))) then
         y(i) = 2
      else
         y(i) = 3
      end if
   end do

   call init_clmm_laplace_problem(problem, y, x, group, re_z, term_q, 3, link=link_logit, &
                                  threshold=threshold_flexible, status=status)
   if (status /= 0) error stop 'clmm_general_example: initialization failed'

   start = [-0.6_dp, 0.8_dp, 0.5_dp, log(0.45_dp), 0.0_dp, log(0.25_dp)]
   call fit_clmm_laplace(problem, fit, start=start, max_iter=100, grad_tol=8.0e-4_dp)
   if (.not. fit%converged) error stop 'clmm_general_example: fit did not converge'

   print '(a,f12.6)', 'logLik: ', fit%loglik
   print '(a,f12.6)', 'beta:   ', fit%beta(1)
   print '(a,2f12.6)', 'random-effect variances: ', fit%re_cov(1, 1), fit%re_cov(2, 2)
end program clmm_general_example
