program test_clmm_general
   use ordinal, only : dp, clmm_problem, clmm_laplace_problem, clmm_laplace_fit_result, &
                       init_clmm_problem, init_clmm_laplace_problem, clmm_nll, clmm_laplace_nll, &
                       clmm_laplace_modes, fit_clmm_laplace, link_logit, threshold_flexible
   implicit none
   integer, parameter :: ng = 6, per_group = 6, n = ng*per_group, k = 3
   integer :: y(n), group(n), group1(n, 1), groups2(n, 2), tq1(1), tq2(2)
   integer :: i, j, g, status
   real(dp) :: x(n, 1), z1(n, 1), z2(n, 2), zcross(n, 2), eta, u, p1, p2, re(ng)
   real(dp) :: par_scalar(4), par_slope(6), par_cross(5), nll_scalar, nll_general, nll_slope, nll_cross
   real(dp), allocatable :: modes(:), condvar(:, :)
   type(clmm_problem) :: scalar_problem
   type(clmm_laplace_problem) :: general_problem, slope_problem, crossed_problem
   type(clmm_laplace_fit_result) :: slope_fit

   do g = 1, ng
      re(g) = 0.8_dp*sin(1.37_dp*real(g, dp))
      do j = 1, per_group
         i = (g - 1)*per_group + j
         group(i) = 20 + 2*g
         group1(i, 1) = group(i)
         x(i, 1) = -1.0_dp + 2.0_dp*real(j - 1, dp)/real(per_group - 1, dp)
         z1(i, 1) = 1.0_dp
         z2(i, 1) = 1.0_dp
         z2(i, 2) = x(i, 1)
         eta = 0.7_dp*x(i, 1) + re(g)
         p1 = logistic(-0.6_dp - eta)
         p2 = logistic(0.9_dp - eta)
         u = real(mod(149*i + 43*g + 17, 997), dp)/997.0_dp
         if (u < p1) then
            y(i) = 1
         else if (u < p2) then
            y(i) = 2
         else
            y(i) = 3
         end if
      end do
   end do

   tq1 = [1]
   call init_clmm_problem(scalar_problem, y, x, group, k, link=link_logit, threshold=threshold_flexible, &
                          nAGQ=1, status=status)
   if (status /= 0) error stop 'scalar Laplace initialization failed'
   call init_clmm_laplace_problem(general_problem, y, x, group1, z1, tq1, k, link=link_logit, &
                                  threshold=threshold_flexible, status=status)
   if (status /= 0) error stop 'general scalar initialization failed'
   par_scalar = [-0.6_dp, 0.9_dp, 0.7_dp, log(0.8_dp)]
   nll_scalar = clmm_nll(par_scalar, scalar_problem)
   nll_general = clmm_laplace_nll(par_scalar, general_problem)
   if (abs(nll_scalar - nll_general) > 2.0e-8_dp) error stop 'general Laplace does not match scalar Laplace'

   call init_clmm_laplace_problem(slope_problem, y, x, group1, z2, [2], k, link=link_logit, &
                                  threshold=threshold_flexible, status=status)
   if (status /= 0) error stop 'random-slope initialization failed'
   par_slope = [-0.6_dp, 0.9_dp, 0.7_dp, log(0.75_dp), 0.12_dp, log(0.30_dp)]
   nll_slope = clmm_laplace_nll(par_slope, slope_problem)
   if (.not. (nll_slope < huge(1.0_dp)/1000.0_dp)) error stop 'random-slope Laplace nll is not finite'
   allocate(modes(slope_problem%nrandom), condvar(slope_problem%nrandom, slope_problem%nrandom))
   call clmm_laplace_modes(slope_problem, par_slope(:2), par_slope(3:3), par_slope(4:), modes, condvar, status)
   if (status /= 0) error stop 'random-slope conditional modes failed'
   if (any([(condvar(i, i), i = 1, size(modes))] <= 0.0_dp)) error stop 'random-slope conditional variance is not positive'
   deallocate(modes, condvar)
   call fit_clmm_laplace(slope_problem, slope_fit, start=par_slope, max_iter=80, grad_tol=5.0e-4_dp)
   if (.not. slope_fit%converged) error stop 'general random-slope Laplace fit did not converge'
   if (slope_fit%max_gradient > 1.0e-3_dp) error stop 'general random-slope fit gradient is too large'

   tq2 = [1, 1]
   do i = 1, n
      groups2(i, 1) = group1(i, 1)
      groups2(i, 2) = 100 + mod(i - 1, 4)
      zcross(i, :) = 1.0_dp
   end do
   call init_clmm_laplace_problem(crossed_problem, y, x, groups2, zcross, tq2, k, link=link_logit, &
                                  threshold=threshold_flexible, status=status)
   if (status /= 0) error stop 'crossed-term initialization failed'
   par_cross = [-0.6_dp, 0.9_dp, 0.7_dp, log(0.7_dp), log(0.4_dp)]
   nll_cross = clmm_laplace_nll(par_cross, crossed_problem)
   if (.not. (nll_cross < huge(1.0_dp)/1000.0_dp)) error stop 'crossed-term Laplace nll is not finite'

   print *, 'test_clmm_general: PASS'
contains
   pure elemental real(dp) function logistic(z) result(p)
      real(dp), intent(in) :: z !! Scalar predictor transformed by the logistic CDF.
      if (z >= 0.0_dp) then
         p = 1.0_dp/(1.0_dp + exp(-z))
      else
         p = exp(z)/(1.0_dp + exp(z))
      end if
   end function logistic
end program test_clmm_general
