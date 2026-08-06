! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the RPEGLMEN computational code.

module rpeglmen_likelihood
  use rpeglmen_kinds, only : dp
  use rpeglmen_math, only : safe_exp, vector_norm2, digamma_value
  implicit none
  private

  public :: exp_negative_log_likelihood
  public :: grad_exp_negative_log_likelihood
  public :: gamma_negative_log_likelihood
  public :: grad_gamma_negative_log_likelihood
  public :: gamma_joint_objective_gradient
  public :: predict_mean, mse_value, mse_gradient

contains

  real(dp) function mse_value(coefficients, a, b) result(value)
    real(dp), intent(in) :: coefficients(:)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp) :: residual(size(b))

    residual = matmul(a, coefficients) - b
    value = 0.5_dp * dot_product(residual, residual)
  end function mse_value

  function mse_gradient(coefficients, a, b) result(gradient)
    real(dp), intent(in) :: coefficients(:)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp) :: gradient(size(coefficients))

    gradient = matmul(transpose(a), matmul(a, coefficients) - b)
  end function mse_gradient

  real(dp) function exp_negative_log_likelihood(coefficients, a, b) result(value)
    real(dp), intent(in) :: coefficients(:)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp) :: eta(size(b))

    eta = matmul(a, coefficients)
    value = sum(eta + b * safe_exp(-eta))
  end function exp_negative_log_likelihood

  function grad_exp_negative_log_likelihood(coefficients, a, b, normalize) result(gradient)
    real(dp), intent(in) :: coefficients(:)
    real(dp), intent(in) :: a(:, :), b(:)
    logical, intent(in), optional :: normalize
    real(dp) :: gradient(size(coefficients))
    real(dp) :: eta(size(b)), nrm
    logical :: normalize_local

    eta = matmul(a, coefficients)
    gradient = matmul(transpose(a), 1.0_dp - b * safe_exp(-eta))

    normalize_local = .false.
    if (present(normalize)) normalize_local = normalize
    if (normalize_local) then
      nrm = vector_norm2(gradient)
      if (nrm > 0.0_dp) gradient = gradient / nrm
    end if
  end function grad_exp_negative_log_likelihood

  real(dp) function gamma_negative_log_likelihood(coefficients, a, b, shape) result(value)
    real(dp), intent(in) :: coefficients(:)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp), intent(in) :: shape
    real(dp) :: eta(size(b))

    eta = matmul(a, coefficients)
    value = sum(-shape * log(shape) + shape * eta + log_gamma(shape) &
      - (shape - 1.0_dp) * log(b) + shape * b * safe_exp(-eta)) / real(size(b), dp)
  end function gamma_negative_log_likelihood

  function grad_gamma_negative_log_likelihood(coefficients, a, b, shape, normalize) result(gradient)
    real(dp), intent(in) :: coefficients(:)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp), intent(in) :: shape
    logical, intent(in), optional :: normalize
    real(dp) :: gradient(size(coefficients))
    real(dp) :: eta(size(b)), nrm
    logical :: normalize_local

    eta = matmul(a, coefficients)
    gradient = matmul(transpose(a), shape * (1.0_dp - b * safe_exp(-eta))) / real(size(b), dp)

    normalize_local = .false.
    if (present(normalize)) normalize_local = normalize
    if (normalize_local) then
      nrm = vector_norm2(gradient)
      if (nrm > 0.0_dp) gradient = gradient / nrm
    end if
  end function grad_gamma_negative_log_likelihood

  subroutine gamma_joint_objective_gradient(theta, a, b, value, gradient)
    real(dp), intent(in) :: theta(:)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp), intent(out) :: value
    real(dp), intent(out) :: gradient(size(theta))
    integer :: p
    real(dp) :: shape, eta(size(b)), mean_shape_gradient

    p = size(theta) - 1
    shape = safe_exp(theta(p + 1))
    eta = matmul(a, theta(1:p))

    value = sum(-shape * log(shape) + shape * eta + log_gamma(shape) &
      - (shape - 1.0_dp) * log(b) + shape * b * safe_exp(-eta)) / real(size(b), dp)

    gradient(1:p) = matmul(transpose(a), shape * (1.0_dp - b * safe_exp(-eta))) / real(size(b), dp)
    mean_shape_gradient = sum(-log(shape) - 1.0_dp + eta + digamma_value(shape) &
      - log(b) + b * safe_exp(-eta)) / real(size(b), dp)
    gradient(p + 1) = shape * mean_shape_gradient
  end subroutine gamma_joint_objective_gradient

  function predict_mean(a, coefficients) result(mu)
    real(dp), intent(in) :: a(:, :), coefficients(:)
    real(dp) :: mu(size(a, 1))

    mu = safe_exp(matmul(a, coefficients))
  end function predict_mean

end module rpeglmen_likelihood
