program test_errors
  use uncorbets, only : dp, torsion_result, effective_bets_result, torsion, &
      effective_bets, sqrtm, uncorbets_invalid_input, uncorbets_not_pos_semidefinite
  implicit none
  real(dp) :: nonsquare(2, 3), indefinite(2, 2), sigma(2, 2), tmat(2, 2)
  real(dp) :: bbad(3)
  type(torsion_result) :: tresult
  type(effective_bets_result) :: eb

  nonsquare = 1.0_dp
  tresult = torsion(nonsquare)
  call assert_true(tresult%status%code == uncorbets_invalid_input, &
      'nonsquare covariance rejected')

  indefinite = reshape([1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp], [2, 2])
  tresult = sqrtm(indefinite)
  call assert_true(tresult%status%code == uncorbets_not_pos_semidefinite, &
      'indefinite square root rejected')

  sigma = reshape([1.0_dp, 0.2_dp, 0.2_dp, 1.0_dp], [2, 2])
  tmat = 0.0_dp
  tmat(1, 1) = 1.0_dp
  tmat(2, 2) = 1.0_dp
  bbad = 1.0_dp / 3.0_dp
  eb = effective_bets(bbad, sigma, tmat)
  call assert_true(eb%status%code == uncorbets_invalid_input, &
      'dimension mismatch rejected')
  print '(a)', 'test_errors: PASS'
contains
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAILED: ' // message
      error stop 1
    end if
  end subroutine assert_true
end program test_errors
