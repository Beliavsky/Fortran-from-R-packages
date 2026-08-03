! SPDX-License-Identifier: GPL-2.0-or-later
program test_stochastic
  use maxlik, only: dp, maxlik_problem, maxlik_control, maxlik_result, initialize_problem, max_lik
  implicit none

  real(dp), parameter :: targets(10, 2) = reshape([ &
    0.5_dp, 1.5_dp, 1.0_dp, 2.0_dp, 1.5_dp, 0.0_dp, 2.5_dp, 1.0_dp, 0.5_dp, 1.5_dp, &
   -1.5_dp,-2.5_dp,-2.0_dp,-1.0_dp,-2.5_dp,-1.5_dp,-2.0_dp,-3.0_dp,-1.0_dp,-2.0_dp], [10, 2])
  type(maxlik_problem) :: problem
  type(maxlik_control) :: control
  type(maxlik_result) :: result
  real(dp) :: optimum(2), start(2)

  optimum = sum(targets, dim=1) / real(size(targets, 1), dp)
  start = [-3.0_dp, 3.0_dp]
  call initialize_problem(problem, 2, objective, size(targets, 1))
  problem%gradient => gradient
  problem%hessian => hessian
  problem%scores => scores

  control%iterlim = 1000
  control%learning_rate = 0.08_dp
  control%batch_size = 4
  control%patience = 80
  control%patience_step = 1
  control%reltol = 1.0e-10_dp
  control%random_seed = 31415
  call max_lik(problem, start, result, 'adam', control)
  call assert_close(result%estimate(1), optimum(1), 8.0e-2_dp, 'Adam x1')
  call assert_close(result%estimate(2), optimum(2), 8.0e-2_dp, 'Adam x2')

  control%iterlim = 800
  control%learning_rate = 0.01_dp
  control%batch_size = size(targets, 1)
  control%sga_momentum = 0.5_dp
  control%patience = 0
  call max_lik(problem, start, result, 'sga', control)
  call assert_true(result%converged, 'SGA converged')
  call assert_close(result%estimate(1), optimum(1), 2.0e-4_dp, 'SGA x1')
  call assert_close(result%estimate(2), optimum(2), 2.0e-4_dp, 'SGA x2')
  print '(a)', 'test_stochastic: PASS'

contains

  subroutine objective(x, value, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    integer :: i
    value = 0.0_dp
    do i = 1, size(targets, 1)
      value = value - 0.5_dp * sum((x - targets(i, :))**2)
    end do
    status = 0
  end subroutine objective

  subroutine gradient(x, g, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    integer, intent(out) :: status
    g = sum(targets, dim=1) - real(size(targets, 1), dp) * x
    status = 0
  end subroutine gradient

  subroutine hessian(x, h, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: h(:, :)
    integer, intent(out) :: status
    h = 0.0_dp
    h(1, 1) = -real(size(targets, 1), dp) + 0.0_dp * x(1)
    h(2, 2) = -real(size(targets, 1), dp)
    status = 0
  end subroutine hessian

  subroutine scores(x, s, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: s(:, :)
    integer, intent(out) :: status
    integer :: i
    do i = 1, size(targets, 1)
      s(i, :) = targets(i, :) - x
    end do
    status = 0
  end subroutine scores

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', 'FAIL: ' // label
      if (allocated(result%estimate)) print '(a,2f14.6)', 'estimate:', result%estimate
      print '(a,2f14.6)', 'optimum:', optimum
      print '(a,i0,2a)', 'code=', result%code, ' message=', trim(result%message)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    call assert_true(abs(actual - expected) <= tolerance, label)
  end subroutine assert_close

end program test_stochastic
