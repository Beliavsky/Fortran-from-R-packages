program test_effective_bets
  use uncorbets, only : dp, torsion_result, effective_bets_result, torsion, effective_bets
  implicit none
  real(dp) :: sigma(4, 4), b(4)
  type(torsion_result) :: tresult
  type(effective_bets_result) :: result

  sigma = reshape([0.040_dp, 0.012_dp, 0.006_dp, 0.004_dp, &
                   0.012_dp, 0.090_dp, 0.015_dp, 0.008_dp, &
                   0.006_dp, 0.015_dp, 0.160_dp, 0.012_dp, &
                   0.004_dp, 0.008_dp, 0.012_dp, 0.250_dp], [4, 4])
  b = 0.25_dp
  tresult = torsion(sigma)
  call assert_true(tresult%status%ok(), 'torsion status')
  result = effective_bets(b, sigma, tresult%matrix)
  call assert_true(result%status%ok(), 'effective bets status')
  call assert_true(abs(sum(result%probability) - 1.0_dp) < 1.0e-10_dp, &
      'probabilities sum to one')
  call assert_true(all(result%probability > 0.0_dp), 'positive probabilities')
  call assert_true(result%enb > 1.0_dp .and. result%enb <= 4.0_dp + 1.0e-10_dp, &
      'ENB bounds')
  print '(a)', 'test_effective_bets: PASS'
contains
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAILED: ' // message
      error stop 1
    end if
  end subroutine assert_true
end program test_effective_bets
