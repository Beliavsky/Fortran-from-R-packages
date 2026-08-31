program test_clm_scale
   use ordinal, only : dp, clm_problem, clm_fit_result, init_clm_problem, fit_clm, &
      link_probit, threshold_equidistant
   implicit none
   integer, parameter :: n = 240, k = 4
   real(dp) :: x(n, 1), sx(n, 1), eta, sigma, u, cuts(3), cdf(3)
   integer :: y(n), i, status
   type(clm_problem) :: problem
   type(clm_fit_result) :: fit
   cuts = [-1.2_dp, 0.0_dp, 1.2_dp]
   do i = 1, n
      x(i, 1) = -1.5_dp + 3.0_dp*real(i - 1, dp)/real(n - 1, dp)
      sx(i, 1) = merge(0.5_dp, -0.5_dp, mod(i, 2) == 0)
      eta = 0.7_dp*x(i, 1)
      sigma = exp(0.35_dp*sx(i, 1))
      cdf = 0.5_dp*erfc(-(cuts - eta)/(sigma*sqrt(2.0_dp)))
      u = real(mod(53*i, 991), dp)/991.0_dp
      if (u < cdf(1)) then
         y(i) = 1
      else if (u < cdf(2)) then
         y(i) = 2
      else if (u < cdf(3)) then
         y(i) = 3
      else
         y(i) = 4
      end if
   end do
   call init_clm_problem(problem, y, x, k, link=link_probit, threshold=threshold_equidistant, scale_x=sx, status=status)
   if (status /= 0) error stop 'scale model initialization failed'
   call fit_clm(problem, fit, max_iter=260, grad_tol=5.0e-6_dp)
   if (.not. fit%converged) then
      print *, 'status=', fit%status, 'grad=', maxval(abs(fit%gradient))
      error stop 'scale clm failed to converge'
   end if
   if (size(fit%zeta) /= 1) error stop 'missing scale coefficient'
   if (fit%theta(1) >= fit%theta(2) .or. fit%theta(2) >= fit%theta(3)) error stop 'equidistant threshold order'
   if (abs((fit%theta(2) - fit%theta(1)) - (fit%theta(3) - fit%theta(2))) > 1.0e-10_dp) error stop 'not equidistant'
   print *, 'test_clm_scale: PASS'
end program test_clm_scale
