! SPDX-License-Identifier: GPL-2.0-or-later
program test_bhhh_and_inference
  use maxlik, only: dp, maxlik_problem, maxlik_control, maxlik_result, initialize_problem, max_lik, &
    robust_covariance_matrix
  implicit none

  real(dp), parameter :: data(8) = [0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp, 2.5_dp, 3.0_dp, 3.5_dp, 4.0_dp]
  type(maxlik_problem) :: problem
  type(maxlik_control) :: control
  type(maxlik_result) :: result
  real(dp), allocatable :: robust(:, :)
  real(dp) :: expected_mu, expected_sigma, start(2)
  integer :: status

  expected_mu = sum(data) / real(size(data), dp)
  expected_sigma = sqrt(sum((data - expected_mu)**2) / real(size(data), dp))
  start = [0.0_dp, 0.0_dp]
  call initialize_problem(problem, 2, objective, size(data))
  problem%gradient => gradient
  problem%hessian => hessian
  problem%scores => scores
  control%iterlim = 400
  control%gradtol = 1.0e-7_dp

  call max_lik(problem, start, result, 'bhhh', control)
  call assert_true(result%converged, 'BHHH converged')
  call assert_close(result%estimate(1), expected_mu, 2.0e-5_dp, 'normal mean')
  call assert_close(exp(result%estimate(2)), expected_sigma, 5.0e-5_dp, 'normal sigma')
  call assert_true(allocated(result%covariance), 'covariance allocated')
  call assert_true(allocated(result%gradient_obs), 'observation scores allocated')
  call robust_covariance_matrix(result%hessian, result%gradient_obs, result%active, robust, status)
  call assert_true(status == 0, 'robust covariance status')
  call assert_true(robust(1, 1) >= 0.0_dp .and. robust(2, 2) >= 0.0_dp, 'robust covariance diagonal')
  print '(a)', 'test_bhhh_and_inference: PASS'

contains

  subroutine objective(x, value, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    real(dp) :: sigma, residual(size(data))
    sigma = exp(x(2))
    residual = (data - x(1)) / sigma
    value = -real(size(data), dp) * (0.5_dp * log(2.0_dp * acos(-1.0_dp)) + x(2)) &
      - 0.5_dp * sum(residual * residual)
    status = 0
  end subroutine objective

  subroutine gradient(x, g, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    integer, intent(out) :: status
    real(dp) :: sigma2, residual(size(data))
    sigma2 = exp(2.0_dp * x(2))
    residual = data - x(1)
    g(1) = sum(residual) / sigma2
    g(2) = -real(size(data), dp) + sum(residual * residual) / sigma2
    status = 0
  end subroutine gradient

  subroutine hessian(x, h, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: h(:, :)
    integer, intent(out) :: status
    real(dp) :: sigma2, residual(size(data))
    sigma2 = exp(2.0_dp * x(2))
    residual = data - x(1)
    h(1, 1) = -real(size(data), dp) / sigma2
    h(1, 2) = -2.0_dp * sum(residual) / sigma2
    h(2, 1) = h(1, 2)
    h(2, 2) = -2.0_dp * sum(residual * residual) / sigma2
    status = 0
  end subroutine hessian

  subroutine scores(x, s, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: s(:, :)
    integer, intent(out) :: status
    real(dp) :: sigma2, residual(size(data))
    sigma2 = exp(2.0_dp * x(2))
    residual = data - x(1)
    s(:, 1) = residual / sigma2
    s(:, 2) = -1.0_dp + residual * residual / sigma2
    status = 0
  end subroutine scores

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', 'FAIL: ' // label
      print '(a,i0,2a)', 'code=', result%code, ' message=', trim(result%message)
      if (allocated(result%estimate)) print '(a,2f14.6)', 'estimate:', result%estimate
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    call assert_true(abs(actual - expected) <= tolerance, label)
  end subroutine assert_close

end program test_bhhh_and_inference
