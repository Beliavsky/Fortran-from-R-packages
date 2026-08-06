program basic_uncorbets
  use uncorbets, only : dp, torsion_result, effective_bets_result, torsion, effective_bets
  implicit none
  real(dp) :: sigma(3, 3), weights(3)
  type(torsion_result) :: tresult
  type(effective_bets_result) :: bets

  sigma = reshape([1.00_dp, 0.35_dp, 0.15_dp, &
                   0.35_dp, 1.50_dp, 0.25_dp, &
                   0.15_dp, 0.25_dp, 2.00_dp], [3, 3])
  weights = 1.0_dp / 3.0_dp
  tresult = torsion(sigma)
  bets = effective_bets(weights, sigma, tresult%matrix)
  write(*, '(a, 3f11.7)') 'Diversification probabilities: ', bets%probability
  write(*, '(a, f11.7)') 'Effective number of bets:      ', bets%enb
end program basic_uncorbets
