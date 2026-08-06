! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the RPEGLMEN computational code.

module rpeglmen_solver
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rpeglmen_kinds, only : dp
  use rpeglmen_types, only : enet_options, fit_result, clear_fit_result, &
    rpe_success, rpe_invalid_input, rpe_no_convergence, rpe_line_search_failure, &
    model_exponential, model_gamma
  use rpeglmen_math, only : vector_norm2, all_finite
  use rpeglmen_penalty, only : prox_en, regularizer_en
  use rpeglmen_likelihood, only : exp_negative_log_likelihood, &
    grad_exp_negative_log_likelihood, gamma_negative_log_likelihood, &
    grad_gamma_negative_log_likelihood
  implicit none
  private

  public :: fit_fixed_model, glmnet_exp_fixed, glm_gamma_net_fixed

contains

  subroutine fit_fixed_model(a, b, model, shape, lambda, options, result, start)
    real(dp), intent(in) :: a(:, :), b(:)
    integer, intent(in) :: model
    real(dp), intent(in) :: shape, lambda
    type(enet_options), intent(in) :: options
    type(fit_result), intent(out) :: result
    real(dp), intent(in), optional :: start(:)

    real(dp), allocatable :: x(:), x_old(:), y(:), z(:), gradient(:), objective_work(:)
    real(dp) :: smooth_y, smooth_z, rhs, step, delta(size(a, 2))
    real(dp) :: objective_value, previous_objective, momentum, momentum_next
    integer :: p, iteration, bt, used
    logical :: accepted, penalize_first

    call clear_fit_result(result)
    p = size(a, 2)

    if (size(a, 1) /= size(b) .or. p < 1 .or. size(b) < 1) then
      call fail(result, rpe_invalid_input, 'incompatible or empty design and response')
      return
    end if
    if (any(b <= 0.0_dp) .or. .not. all_finite(b)) then
      call fail(result, rpe_invalid_input, 'responses must be finite and strictly positive')
      return
    end if
    if (options%alpha < 0.0_dp .or. options%alpha > 1.0_dp .or. lambda < 0.0_dp) then
      call fail(result, rpe_invalid_input, 'alpha must be in [0,1] and lambda nonnegative')
      return
    end if
    if (options%max_iter < 1 .or. options%max_backtrack < 1 .or. options%initial_step <= 0.0_dp &
      .or. options%backtrack <= 0.0_dp .or. options%backtrack >= 1.0_dp) then
      call fail(result, rpe_invalid_input, 'invalid iteration or line-search options')
      return
    end if
    if (model == model_gamma .and. shape <= 0.0_dp) then
      call fail(result, rpe_invalid_input, 'Gamma shape must be positive')
      return
    end if
    if (present(start)) then
      if (size(start) /= p) then
        call fail(result, rpe_invalid_input, 'starting vector has wrong size')
        return
      end if
    end if

    allocate(x(p), x_old(p), y(p), z(p), gradient(p), objective_work(options%max_iter))
    x = 0.0_dp
    if (present(start)) x = start
    x_old = x
    y = x
    step = options%initial_step
    momentum = 1.0_dp
    previous_objective = huge(1.0_dp)
    used = 0
    penalize_first = options%penalize_intercept .or. .not. options%has_intercept

    do iteration = 1, options%max_iter
      call smooth_value_gradient(model, y, a, b, shape, options%normalize_gradient, smooth_y, gradient)
      if (.not. ieee_is_finite(smooth_y) .or. .not. all_finite(gradient)) then
        call fail(result, rpe_no_convergence, 'non-finite smooth objective or gradient')
        result%coefficients = x
        return
      end if

      accepted = .false.
      do bt = 1, options%max_backtrack
        z = prox_en(y - step * gradient, step, lambda, options%alpha, &
          penalize_first, options%source_proximal)
        call smooth_value(model, z, a, b, shape, smooth_z)
        delta = z - y
        rhs = smooth_y + dot_product(gradient, delta) &
          + dot_product(delta, delta) / (2.0_dp * step)
        if (ieee_is_finite(smooth_z) .and. smooth_z <= rhs + 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(rhs))) then
          accepted = .true.
          exit
        end if
        step = step * options%backtrack
        if (step <= tiny(1.0_dp)) exit
      end do

      if (.not. accepted) then
        call fail(result, rpe_line_search_failure, 'proximal line search failed')
        result%coefficients = x
        result%iterations = iteration - 1
        return
      end if

      objective_value = smooth_z + lambda * regularizer_en(z, options%alpha, penalize_first)
      used = used + 1
      objective_work(used) = objective_value

      if (iteration > 1) then
        if (abs(objective_value - previous_objective) <= options%abs_tol &
          + options%rel_tol * max(1.0_dp, abs(previous_objective))) then
          x = z
          result%converged = .true.
          exit
        end if
      end if

      previous_objective = objective_value
      if (options%use_fista) then
        momentum_next = 0.5_dp * (1.0_dp + sqrt(1.0_dp + 4.0_dp * momentum * momentum))
        y = z + ((momentum - 1.0_dp) / momentum_next) * (z - x)
        momentum = momentum_next
      else
        y = z
      end if
      x_old = x
      x = z
      if (vector_norm2(x - x_old) <= options%abs_tol &
        + options%rel_tol * max(1.0_dp, vector_norm2(x_old))) then
        result%converged = .true.
        exit
      end if
    end do

    x = z
    result%coefficients = x
    result%objective = objective_work(1:used)
    result%iterations = used
    result%selected_lambda = lambda
    result%shape = shape
    if (result%converged) then
      result%status = rpe_success
      result%message = 'converged'
    else
      result%status = rpe_no_convergence
      result%message = 'maximum iterations reached'
    end if
  end subroutine fit_fixed_model

  subroutine glmnet_exp_fixed(a, b, lambda, result, options, start)
    real(dp), intent(in) :: a(:, :), b(:), lambda
    type(fit_result), intent(out) :: result
    type(enet_options), intent(in), optional :: options
    real(dp), intent(in), optional :: start(:)
    type(enet_options) :: local_options

    local_options = enet_options()
    if (present(options)) local_options = options
    call fit_fixed_model(a, b, model_exponential, 1.0_dp, lambda, local_options, result, start)
  end subroutine glmnet_exp_fixed

  subroutine glm_gamma_net_fixed(a, b, shape, lambda, result, options, start)
    real(dp), intent(in) :: a(:, :), b(:), shape, lambda
    type(fit_result), intent(out) :: result
    type(enet_options), intent(in), optional :: options
    real(dp), intent(in), optional :: start(:)
    type(enet_options) :: local_options

    local_options = enet_options()
    if (present(options)) local_options = options
    call fit_fixed_model(a, b, model_gamma, shape, lambda, local_options, result, start)
  end subroutine glm_gamma_net_fixed

  subroutine smooth_value_gradient(model, x, a, b, shape, normalize, value, gradient)
    integer, intent(in) :: model
    real(dp), intent(in) :: x(:), a(:, :), b(:), shape
    logical, intent(in) :: normalize
    real(dp), intent(out) :: value, gradient(:)

    select case (model)
    case (model_exponential)
      value = exp_negative_log_likelihood(x, a, b)
      gradient = grad_exp_negative_log_likelihood(x, a, b, normalize)
    case (model_gamma)
      value = gamma_negative_log_likelihood(x, a, b, shape)
      gradient = grad_gamma_negative_log_likelihood(x, a, b, shape, normalize)
    case default
      value = huge(1.0_dp)
      gradient = 0.0_dp
    end select
  end subroutine smooth_value_gradient

  subroutine smooth_value(model, x, a, b, shape, value)
    integer, intent(in) :: model
    real(dp), intent(in) :: x(:), a(:, :), b(:), shape
    real(dp), intent(out) :: value

    select case (model)
    case (model_exponential)
      value = exp_negative_log_likelihood(x, a, b)
    case (model_gamma)
      value = gamma_negative_log_likelihood(x, a, b, shape)
    case default
      value = huge(1.0_dp)
    end select
  end subroutine smooth_value

  subroutine fail(result, status, message)
    type(fit_result), intent(inout) :: result
    integer, intent(in) :: status
    character(len=*), intent(in) :: message

    result%status = status
    result%message = message
    result%converged = .false.
  end subroutine fail

end module rpeglmen_solver
