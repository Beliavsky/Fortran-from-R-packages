! SPDX-License-Identifier: GPL-2.0-or-later
program normal_mle
  use maxlik, only: dp, maxlik_problem, maxlik_result, initialize_problem, max_lik
  implicit none

  real(dp), parameter :: observations(10) = [1.2_dp, 0.7_dp, 1.8_dp, 2.1_dp, 1.5_dp, &
    0.9_dp, 1.1_dp, 1.7_dp, 2.0_dp, 1.4_dp]
  type(maxlik_problem) :: problem
  type(maxlik_result) :: result

  call initialize_problem(problem, 2, log_likelihood, size(observations))
  problem%gradient => likelihood_gradient
  problem%hessian => likelihood_hessian
  problem%scores => observation_scores
  call max_lik(problem, [1.0_dp, log(0.5_dp)], result, 'nr')

  print '(a,f12.6)', 'log likelihood: ', result%maximum
  print '(a,f12.6)', 'mu:             ', result%estimate(1)
  print '(a,f12.6)', 'sigma:          ', exp(result%estimate(2))
  print '(a,2f12.6)', 'standard errors:', result%std_error

contains

  subroutine log_likelihood(x, value, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    real(dp) :: residual(size(observations)), sigma
    sigma = exp(x(2))
    residual = (observations - x(1)) / sigma
    value = -real(size(observations), dp) * (0.5_dp * log(2.0_dp * acos(-1.0_dp)) + x(2)) &
      - 0.5_dp * sum(residual**2)
    status = 0
  end subroutine log_likelihood

  subroutine likelihood_gradient(x, g, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    integer, intent(out) :: status
    real(dp) :: residual(size(observations)), sigma2
    sigma2 = exp(2.0_dp * x(2))
    residual = observations - x(1)
    g(1) = sum(residual) / sigma2
    g(2) = -real(size(observations), dp) + sum(residual**2) / sigma2
    status = 0
  end subroutine likelihood_gradient

  subroutine likelihood_hessian(x, h, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: h(:, :)
    integer, intent(out) :: status
    real(dp) :: residual(size(observations)), sigma2
    sigma2 = exp(2.0_dp * x(2))
    residual = observations - x(1)
    h(1, 1) = -real(size(observations), dp) / sigma2
    h(1, 2) = -2.0_dp * sum(residual) / sigma2
    h(2, 1) = h(1, 2)
    h(2, 2) = -2.0_dp * sum(residual**2) / sigma2
    status = 0
  end subroutine likelihood_hessian

  subroutine observation_scores(x, scores, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: scores(:, :)
    integer, intent(out) :: status
    real(dp) :: residual(size(observations)), sigma2
    sigma2 = exp(2.0_dp * x(2))
    residual = observations - x(1)
    scores(:, 1) = residual / sigma2
    scores(:, 2) = -1.0_dp + residual**2 / sigma2
    status = 0
  end subroutine observation_scores

end program normal_mle
