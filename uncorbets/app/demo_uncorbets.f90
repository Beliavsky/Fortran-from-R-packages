program demo_uncorbets
  use uncorbets, only : dp, torsion_result, effective_bets_result, &
      max_effective_bets_result, torsion, effective_bets, max_effective_bets
  implicit none
  real(dp) :: sigma(4, 4), weights(4)
  type(torsion_result) :: tresult
  type(effective_bets_result) :: initial
  type(max_effective_bets_result) :: optimum
  integer :: i

  sigma = reshape([0.040_dp, 0.012_dp, 0.006_dp, 0.004_dp, &
                   0.012_dp, 0.090_dp, 0.015_dp, 0.008_dp, &
                   0.006_dp, 0.015_dp, 0.160_dp, 0.012_dp, &
                   0.004_dp, 0.008_dp, 0.012_dp, 0.250_dp], [4, 4])
  weights = 0.25_dp
  tresult = torsion(sigma, model='minimum-torsion', method='exact')
  if (.not. tresult%status%ok()) error stop tresult%status%message
  initial = effective_bets(weights, sigma, tresult%matrix)
  optimum = max_effective_bets(weights, sigma, tresult%matrix)

  write(*, '(a, f10.6)') 'Equal-weight ENB: ', initial%enb
  write(*, '(a, f10.6)') 'Maximum ENB:      ', optimum%enb
  write(*, '(a)') 'Optimized weights:'
  do i = 1, size(optimum%weights)
    write(*, '(i3, 2x, f12.8)') i, optimum%weights(i)
  end do
end program demo_uncorbets
