module test_model_context
  use sgt, only : dp
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
  implicit none
  real(dp), allocatable :: yy(:), zz(:)
contains
  subroutine regression_model(theta, index, x, mu, sigma, lambda, p, q, status)
    real(dp), intent(in) :: theta(:)
    integer, intent(in) :: index
    real(dp), intent(out) :: x, mu, sigma, lambda, p, q
    integer, intent(out) :: status
    x = yy(index) - theta(1) - theta(2) * zz(index)
    mu = 0.0_dp
    sigma = exp(theta(3))
    lambda = 0.0_dp
    p = 2.0_dp
    q = ieee_value(0.0_dp, ieee_positive_inf)
    status = 0
  end subroutine regression_model
end module test_model_context

program test_model_mle
  use sgt
  use test_model_context
  implicit none
  integer, parameter :: n = 101
  real(dp) :: theta0(3), u, e
  type(sgt_mle_result) :: fit
  integer :: i, k, failures
  failures = 0
  allocate(yy(n), zz(n))
  do i = 1, n
    zz(i) = -1.0_dp + 2.0_dp * real(i - 1, dp) / real(n - 1, dp)
    k = modulo(37 * (i - 1), n) + 1
    u = (real(k, dp) - 0.5_dp) / real(n, dp)
    e = qsgt(u, 0.0_dp, 0.55_dp, 0.0_dp, 2.0_dp, &
      ieee_value(0.0_dp, ieee_positive_inf))
    yy(i) = 1.1_dp - 0.7_dp * zz(i) + e
  end do
  theta0 = [0.8_dp, -0.4_dp, log(0.8_dp)]
  call sgt_mle_model(n, theta0, regression_model, fit, max_iter=500)
  if (fit%convcode /= 0) failures = failures + 1
  if (abs(fit%estimate(1) - 1.1_dp) > 0.03_dp) failures = failures + 1
  if (abs(fit%estimate(2) + 0.7_dp) > 0.06_dp) failures = failures + 1
  if (abs(exp(fit%estimate(3)) - 0.55_dp) > 0.025_dp) failures = failures + 1
  if (failures /= 0) then
    print *, fit%estimate, exp(fit%estimate(3))
    error stop 1
  end if
  print '(a)', 'test_model_mle: PASS'
end program test_model_mle
