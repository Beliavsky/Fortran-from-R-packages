program constrained_portfolio
   use risk_parity_portfolio_mod
   implicit none
   integer, parameter :: n = 8
   real(dp) :: vol(n), corr(n, n), sigma(n, n)
   real(dp) :: group_constraint(1, n), group_limit(1), lower(n), upper(n)
   type(risk_parity_result) :: result
   integer :: i, j

   vol = [0.05_dp, 0.05_dp, 0.07_dp, 0.10_dp, 0.15_dp, 0.15_dp, 0.15_dp, 0.18_dp]
   corr = reshape([ &
      100, 80, 60,-20,-10,-20,-20,-20, &
       80,100, 40,-20,-20,-10,-20,-20, &
       60, 40,100, 50, 30, 20, 20, 30, &
      -20,-20, 50,100, 60, 60, 50, 60, &
      -10,-20, 30, 60,100, 90, 70, 70, &
      -20,-10, 20, 60, 90,100, 60, 70, &
      -20,-20, 20, 50, 70, 60,100, 70, &
      -20,-20, 30, 60, 70, 70, 70,100], [n, n]) / 100.0_dp
   do j = 1, n
      do i = 1, n
         sigma(i, j) = corr(i, j) * vol(i) * vol(j)
      end do
   end do

   ! Require at least 30 percent in assets 5 through 8:
   ! -sum(w(5:8)) <= -0.30.
   group_constraint = 0.0_dp
   group_constraint(1, 5:8) = -1.0_dp
   group_limit = -0.30_dp
   lower = 0.0_dp
   upper = 0.50_dp

   call risk_parity_portfolio(sigma, result, lower=lower, upper=upper, &
                              dmat=group_constraint, dvec=group_limit, &
                              formulation=FORM_RC_OVER_VAR_VS_B)
   if (result%status /= RPP_OK) error stop 'constrained risk parity solver failed'

   write(*, '(a,8f10.6)') 'weights: ', result%weights
   write(*, '(a,8f10.6)') 'relative risk contributions: ', &
                          result%relative_risk_contribution
   write(*, '(a,f10.6)') 'weight in assets 5:8: ', sum(result%weights(5:8))
   write(*, '(a,l1)') 'feasible: ', result%feasible
end program constrained_portfolio
