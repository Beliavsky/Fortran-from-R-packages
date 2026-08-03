! SPDX-License-Identifier: GPL-2.0-or-later
program adam_scores
  use maxlik, only: dp, maxlik_problem, maxlik_control, maxlik_result, initialize_problem, max_lik
  implicit none

  real(dp), parameter :: targets(8) = [0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp, 2.5_dp, 3.0_dp, 3.5_dp, 4.0_dp]
  type(maxlik_problem) :: problem
  type(maxlik_control) :: control
  type(maxlik_result) :: result

  call initialize_problem(problem, 1, objective, size(targets))
  problem%scores => scores
  control%iterlim = 1000
  control%learning_rate = 0.08_dp
  control%batch_size = 4
  control%patience = 100
  call max_lik(problem, [0.0_dp], result, 'adam', control)

  print '(a,f12.6)', 'Adam estimate: ', result%estimate(1)
  print '(a,f12.6)', 'sample mean:   ', sum(targets) / real(size(targets), dp)

contains

  subroutine objective(x, value, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    value = -0.5_dp * sum((targets - x(1))**2)
    status = 0
  end subroutine objective

  subroutine scores(x, score, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: score(:, :)
    integer, intent(out) :: status
    score(:, 1) = targets - x(1)
    status = 0
  end subroutine scores

end program adam_scores
