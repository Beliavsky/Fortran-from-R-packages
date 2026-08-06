program test_max_effective_bets
  use uncorbets, only : dp, torsion_result, effective_bets_result, &
      max_effective_bets_result, status_type, torsion, effective_bets, &
      effective_bets_gradient, max_effective_bets
  implicit none
  real(dp) :: sigma(4, 4), x0(4), xp(4), xm(4), numeric_gradient(4)
  real(dp), allocatable :: analytic_gradient(:)
  real(dp) :: analytic_enb, h
  integer :: i
  type(torsion_result) :: tresult
  type(effective_bets_result) :: initial
  type(max_effective_bets_result) :: optimum
  type(effective_bets_result) :: plus_result, minus_result
  type(status_type) :: gradient_status

  sigma = reshape([0.040_dp, 0.012_dp, 0.006_dp, 0.004_dp, &
                   0.012_dp, 0.090_dp, 0.015_dp, 0.008_dp, &
                   0.006_dp, 0.015_dp, 0.160_dp, 0.012_dp, &
                   0.004_dp, 0.008_dp, 0.012_dp, 0.250_dp], [4, 4])
  x0 = 0.25_dp
  tresult = torsion(sigma)
  call assert_true(tresult%status%ok(), 'torsion status')
  initial = effective_bets(x0, sigma, tresult%matrix)
  call effective_bets_gradient(x0, sigma, tresult%matrix, analytic_enb, &
      analytic_gradient, gradient_status)
  call assert_true(gradient_status%ok(), 'analytic gradient status')
  h = 1.0e-6_dp
  do i = 1, 4
    xp = x0
    xm = x0
    xp(i) = xp(i) + h
    xm(i) = xm(i) - h
    plus_result = effective_bets(xp, sigma, tresult%matrix)
    minus_result = effective_bets(xm, sigma, tresult%matrix)
    numeric_gradient(i) = (plus_result%enb - minus_result%enb) / (2.0_dp * h)
  end do
  call assert_true(abs(analytic_enb - initial%enb) < 1.0e-12_dp, &
      'gradient ENB value')
  call assert_true(maxval(abs(analytic_gradient - numeric_gradient)) < 1.0e-5_dp, &
      'analytic gradient matches finite differences')
  optimum = max_effective_bets(x0, sigma, tresult%matrix, &
      tolerance=1.0e-11_dp, maxeval=10000, maxiter=5000)
  call assert_true(optimum%status%ok(), 'optimizer status')
  call assert_true(optimum%converged, 'optimizer convergence')
  call assert_true(abs(sum(optimum%weights) - 1.0_dp) < 1.0e-12_dp, &
      'weights sum to one')
  call assert_true(all(optimum%weights >= 0.0_dp), 'nonnegative weights')
  call assert_true(optimum%enb >= initial%enb - 1.0e-10_dp, 'objective improves')
  call assert_true(optimum%enb > 3.999_dp, 'maximum ENB approaches dimension')
  call assert_true(size(optimum%gradient) == 4 .and. &
      all(shape(optimum%hessian) == [4, 4]), 'diagnostic shapes')
  print '(a)', 'test_max_effective_bets: PASS'
contains
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAILED: ' // message
      error stop 1
    end if
  end subroutine assert_true
end program test_max_effective_bets
