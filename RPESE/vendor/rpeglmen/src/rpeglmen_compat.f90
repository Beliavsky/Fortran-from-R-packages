! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the RPEGLMEN computational code.

module rpeglmen_compat
  use rpeglmen_kinds, only : dp
  use rpeglmen_types, only : enet_options, fit_result
  use rpeglmen_likelihood, only : exp_negative_log_likelihood, &
    grad_exp_negative_log_likelihood
  use rpeglmen_solver, only : glmnet_exp_fixed, glm_gamma_net_fixed
  use rpeglmen_cv, only : glmnet_exp, fit_glm_gamma_net
  implicit none
  private

  public :: fit_glmGammaNet, cv_glmGammaNet, glmGammaNet, fitGlmCv
  public :: ExpNegativeLogLikelihood_cpp, GradExpNegativeLogLikelihood_cpp
  public :: ProxGradDescent_cpp

contains

  subroutine fit_glmGammaNet(a, b, result, options, exponential_dist)
    real(dp), intent(in) :: a(:, :), b(:)
    type(fit_result), intent(out) :: result
    type(enet_options), intent(in), optional :: options
    logical, intent(in), optional :: exponential_dist
    logical :: use_exponential

    use_exponential = .false.
    if (present(exponential_dist)) use_exponential = exponential_dist
    if (use_exponential) then
      call glmnet_exp(a, b, result, options)
    else
      call fit_glm_gamma_net(a, b, result, options)
    end if
  end subroutine fit_glmGammaNet

  subroutine cv_glmGammaNet(a, b, result, options)
    real(dp), intent(in) :: a(:, :), b(:)
    type(fit_result), intent(out) :: result
    type(enet_options), intent(in), optional :: options

    call fit_glm_gamma_net(a, b, result, options)
  end subroutine cv_glmGammaNet

  subroutine glmGammaNet(a, b, lambda_en, shape, result, options, start)
    real(dp), intent(in) :: a(:, :), b(:), lambda_en, shape
    type(fit_result), intent(out) :: result
    type(enet_options), intent(in), optional :: options
    real(dp), intent(in), optional :: start(:)

    call glm_gamma_net_fixed(a, b, shape, lambda_en, result, options, start)
  end subroutine glmGammaNet

  subroutine fitGlmCv(a, b, result, options)
    real(dp), intent(in) :: a(:, :), b(:)
    type(fit_result), intent(out) :: result
    type(enet_options), intent(in), optional :: options

    call glmnet_exp(a, b, result, options)
  end subroutine fitGlmCv

  real(dp) function ExpNegativeLogLikelihood_cpp(x, a, b) result(value)
    real(dp), intent(in) :: x(:), a(:, :), b(:)

    value = exp_negative_log_likelihood(x, a, b)
  end function ExpNegativeLogLikelihood_cpp

  function GradExpNegativeLogLikelihood_cpp(x, a, b, normalize_gradient) result(gradient)
    real(dp), intent(in) :: x(:), a(:, :), b(:)
    logical, intent(in), optional :: normalize_gradient
    real(dp) :: gradient(size(x))
    logical :: normalize_local

    normalize_local = .false.
    if (present(normalize_gradient)) normalize_local = normalize_gradient
    gradient = grad_exp_negative_log_likelihood(x, a, b, normalize_local)
  end function GradExpNegativeLogLikelihood_cpp

  subroutine ProxGradDescent_cpp(a, b, lambda, result, options)
    real(dp), intent(in) :: a(:, :), b(:), lambda
    type(fit_result), intent(out) :: result
    type(enet_options), intent(in), optional :: options

    call glmnet_exp_fixed(a, b, lambda, result, options)
  end subroutine ProxGradDescent_cpp

end module rpeglmen_compat
