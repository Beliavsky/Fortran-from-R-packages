program demo_risk_parity
   use risk_parity_portfolio_mod
   implicit none
   integer, parameter :: n = 5
   real(dp) :: sigma(n, n), budgets(n)
   type(risk_parity_result) :: result
   integer :: i

   sigma = reshape([ &
      0.0400_dp, 0.0060_dp, 0.0040_dp, 0.0020_dp, 0.0010_dp, &
      0.0060_dp, 0.0625_dp, 0.0075_dp, 0.0030_dp, 0.0020_dp, &
      0.0040_dp, 0.0075_dp, 0.0900_dp, 0.0090_dp, 0.0045_dp, &
      0.0020_dp, 0.0030_dp, 0.0090_dp, 0.1225_dp, 0.0105_dp, &
      0.0010_dp, 0.0020_dp, 0.0045_dp, 0.0105_dp, 0.1600_dp], [n, n])
   budgets = [0.30_dp, 0.25_dp, 0.20_dp, 0.15_dp, 0.10_dp]

   call risk_parity_portfolio(sigma, result, b=budgets)
   if (result%status /= RPP_OK) error stop 'risk parity solver failed'

   write(*, '(a)') ' asset       weight       target rrc     achieved rrc'
   do i = 1, n
      write(*, '(i6,3f17.8)') i, result%weights(i), budgets(i), &
                              result%relative_risk_contribution(i)
   end do
   write(*, '(a,f17.8)') ' portfolio variance: ', result%variance
   write(*, '(a,i0)') ' iterations: ', result%iterations
end program demo_risk_parity
