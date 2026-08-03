! SPDX-License-Identifier: GPL-2.0-or-later
program logistic_bhhh
  use maxlik, only: dp, maxlik_problem, maxlik_control, maxlik_result, initialize_problem, max_lik
  implicit none

  real(dp), parameter :: xdata(12) = [-2.0_dp, -1.5_dp, -1.0_dp, -0.5_dp, 0.0_dp, 0.5_dp, &
    1.0_dp, 1.5_dp, 2.0_dp, 2.5_dp, 3.0_dp, 3.5_dp]
  real(dp), parameter :: ydata(12) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
    0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
  type(maxlik_problem) :: problem
  type(maxlik_control) :: control
  type(maxlik_result) :: result

  call initialize_problem(problem, 2, objective, size(xdata))
  problem%scores => scores
  control%iterlim = 500
  call max_lik(problem, [0.0_dp, 0.0_dp], result, 'bhhh', control)

  print '(a,2f12.6)', 'logit coefficients:', result%estimate
  print '(a,f12.6)', 'log likelihood:    ', result%maximum

contains

  pure real(dp) function logistic(z) result(p)
    real(dp), intent(in) :: z
    if (z >= 0.0_dp) then
      p = 1.0_dp / (1.0_dp + exp(-z))
    else
      p = exp(z) / (1.0_dp + exp(z))
    end if
  end function logistic

  subroutine objective(beta, value, status)
    real(dp), intent(in) :: beta(:)
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    real(dp) :: p
    integer :: i
    value = 0.0_dp
    do i = 1, size(xdata)
      p = min(1.0_dp - epsilon(1.0_dp), max(epsilon(1.0_dp), logistic(beta(1) + beta(2) * xdata(i))))
      value = value + ydata(i) * log(p) + (1.0_dp - ydata(i)) * log(1.0_dp - p)
    end do
    status = 0
  end subroutine objective

  subroutine scores(beta, score, status)
    real(dp), intent(in) :: beta(:)
    real(dp), intent(out) :: score(:, :)
    integer, intent(out) :: status
    real(dp) :: residual
    integer :: i
    do i = 1, size(xdata)
      residual = ydata(i) - logistic(beta(1) + beta(2) * xdata(i))
      score(i, :) = [residual, residual * xdata(i)]
    end do
    status = 0
  end subroutine scores

end program logistic_bhhh
