program test_clmm_quadrature
   use ordinal, only : dp, clmm_problem, init_clmm_problem, clmm_nll, link_logit, threshold_flexible
   implicit none
   integer, parameter :: ng = 8, per_group = 10, n = ng*per_group, k = 3
   integer :: y(n), group(n), g, j, i, status
   real(dp) :: x(n, 1), eta, u, p1, p2, re(ng), par(4), nla, nagq, nghq
   type(clmm_problem) :: pla, pagq, pghq

   do g = 1, ng
      re(g) = 0.85_dp*sin(1.31_dp*real(g, dp))
      do j = 1, per_group
         i = (g - 1)*per_group + j
         group(i) = 10*g + 3
         x(i, 1) = -1.2_dp + 2.4_dp*real(j - 1, dp)/real(per_group - 1, dp)
         eta = 0.75_dp*x(i, 1) + re(g)
         p1 = logistic(-0.55_dp - eta)
         p2 = logistic(0.85_dp - eta)
         u = real(mod(173*i + 97*g + 31, 1009), dp)/1009.0_dp
         if (u < p1) then
            y(i) = 1
         else if (u < p2) then
            y(i) = 2
         else
            y(i) = 3
         end if
      end do
   end do

   call init_clmm_problem(pla, y, x, group, k, link=link_logit, threshold=threshold_flexible, nAGQ=1, status=status)
   if (status /= 0) error stop 'Laplace initialization failed'
   call init_clmm_problem(pagq, y, x, group, k, link=link_logit, threshold=threshold_flexible, nAGQ=7, status=status)
   if (status /= 0) error stop 'adaptive GHQ initialization failed'
   call init_clmm_problem(pghq, y, x, group, k, link=link_logit, threshold=threshold_flexible, nAGQ=-7, status=status)
   if (status /= 0) error stop 'nonadaptive GHQ initialization failed'

   par = [-0.55_dp, 0.85_dp, 0.75_dp, log(0.85_dp)]
   nla = clmm_nll(par, pla)
   nagq = clmm_nll(par, pagq)
   nghq = clmm_nll(par, pghq)
   if (.not. (nla < huge(1.0_dp)/1000.0_dp)) error stop 'Laplace nll is not finite'
   if (.not. (nagq < huge(1.0_dp)/1000.0_dp)) error stop 'adaptive GHQ nll is not finite'
   if (.not. (nghq < huge(1.0_dp)/1000.0_dp)) error stop 'nonadaptive GHQ nll is not finite'
   if (abs(nagq - nghq) > 3.0e-2_dp) error stop 'adaptive and nonadaptive GHQ disagree unexpectedly'
   if (abs(nla - nagq) > 0.25_dp) error stop 'Laplace approximation is unexpectedly far from seven-point AGQ'

   print *, 'test_clmm_quadrature: PASS'
contains
   pure elemental real(dp) function logistic(z) result(p)
      real(dp), intent(in) :: z !! Scalar predictor transformed by the logistic CDF.
      if (z >= 0.0_dp) then
         p = 1.0_dp/(1.0_dp + exp(-z))
      else
         p = exp(z)/(1.0_dp + exp(z))
      end if
   end function logistic
end program test_clmm_quadrature
