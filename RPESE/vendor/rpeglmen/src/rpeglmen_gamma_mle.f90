! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the RPEGLMEN computational code.

module rpeglmen_gamma_mle
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rpeglmen_kinds, only : dp
  use rpeglmen_types, only : enet_options, fit_result, clear_fit_result, &
    rpe_success, rpe_invalid_input, rpe_no_convergence, rpe_numerical_failure
  use rpeglmen_math, only : vector_norm2, all_finite, approximately_constant_one
  use rpeglmen_likelihood, only : gamma_joint_objective_gradient
  implicit none
  private

  public :: fit_glm_gamma_mle

contains

  subroutine fit_glm_gamma_mle(a, b, result, options, start_coefficients, start_shape)
    real(dp), intent(in) :: a(:, :), b(:)
    type(fit_result), intent(out) :: result
    type(enet_options), intent(in), optional :: options
    real(dp), intent(in), optional :: start_coefficients(:)
    real(dp), intent(in), optional :: start_shape

    type(enet_options) :: local_options
    real(dp), allocatable :: theta(:), theta_new(:), gradient(:), gradient_new(:)
    real(dp), allocatable :: hessian_inverse(:, :), identity(:, :), direction(:)
    real(dp), allocatable :: s(:), y(:), left(:, :), right(:, :), objective_work(:)
    real(dp) :: value, value_new, step, directional, ys, rho
    real(dp) :: mean_b, variance_b, shape0
    integer :: p, q, iteration, ls, i, used
    logical :: accepted

    call clear_fit_result(result)
    local_options = enet_options()
    if (present(options)) local_options = options

    p = size(a, 2)
    q = p + 1
    if (size(a, 1) /= size(b) .or. size(b) < 2 .or. p < 1) then
      call fail(result, rpe_invalid_input, 'incompatible or empty design and response')
      return
    end if
    if (any(b <= 0.0_dp) .or. .not. all_finite(b)) then
      call fail(result, rpe_invalid_input, 'Gamma responses must be finite and strictly positive')
      return
    end if
    if (local_options%max_iter < 1 .or. local_options%max_backtrack < 1 &
      .or. local_options%backtrack <= 0.0_dp .or. local_options%backtrack >= 1.0_dp) then
      call fail(result, rpe_invalid_input, 'invalid iteration or line-search options')
      return
    end if
    if (present(start_coefficients)) then
      if (size(start_coefficients) /= p) then
        call fail(result, rpe_invalid_input, 'starting coefficient vector has wrong size')
        return
      end if
    end if
    if (present(start_shape)) then
      if (start_shape <= 0.0_dp) then
        call fail(result, rpe_invalid_input, 'starting shape must be positive')
        return
      end if
    end if

    allocate(theta(q), theta_new(q), gradient(q), gradient_new(q), direction(q), s(q), y(q))
    allocate(hessian_inverse(q, q), identity(q, q), left(q, q), right(q, q))
    allocate(objective_work(local_options%max_iter))

    theta = 0.0_dp
    if (present(start_coefficients)) then
      theta(1:p) = start_coefficients
    else if (local_options%has_intercept .and. approximately_constant_one(a(:, 1))) then
      theta(1) = log(sum(b) / real(size(b), dp))
    end if

    mean_b = sum(b) / real(size(b), dp)
    variance_b = sum((b - mean_b)**2) / real(size(b) - 1, dp)
    shape0 = max(0.1_dp, min(1.0e4_dp, mean_b * mean_b / max(variance_b, tiny(1.0_dp))))
    if (present(start_shape)) shape0 = start_shape
    theta(q) = log(shape0)

    identity = 0.0_dp
    do i = 1, q
      identity(i, i) = 1.0_dp
    end do
    hessian_inverse = identity

    call gamma_joint_objective_gradient(theta, a, b, value, gradient)
    if (.not. ieee_is_finite(value) .or. .not. all_finite(gradient)) then
      call fail(result, rpe_numerical_failure, 'non-finite initial Gamma likelihood')
      return
    end if

    used = 0
    do iteration = 1, local_options%max_iter
      if (vector_norm2(gradient) <= local_options%abs_tol &
        + local_options%rel_tol * max(1.0_dp, vector_norm2(theta))) then
        result%converged = .true.
        exit
      end if

      direction = -matmul(hessian_inverse, gradient)
      directional = dot_product(gradient, direction)
      if (directional >= -sqrt(epsilon(1.0_dp)) * vector_norm2(gradient) * max(1.0_dp, vector_norm2(direction))) then
        direction = -gradient
        directional = -dot_product(gradient, gradient)
        hessian_inverse = identity
      end if

      step = 1.0_dp
      accepted = .false.
      do ls = 1, local_options%max_backtrack
        theta_new = theta + step * direction
        theta_new(q) = max(-12.0_dp, min(20.0_dp, theta_new(q)))
        call gamma_joint_objective_gradient(theta_new, a, b, value_new, gradient_new)
        if (ieee_is_finite(value_new) .and. all_finite(gradient_new)) then
          if (value_new <= value + 1.0e-4_dp * step * directional) then
            accepted = .true.
            exit
          end if
        end if
        step = step * local_options%backtrack
      end do

      if (.not. accepted) then
        call fail(result, rpe_no_convergence, 'Gamma MLE line search failed')
        exit
      end if

      s = theta_new - theta
      y = gradient_new - gradient
      ys = dot_product(y, s)
      if (ys > sqrt(epsilon(1.0_dp)) * vector_norm2(y) * max(vector_norm2(s), tiny(1.0_dp))) then
        rho = 1.0_dp / ys
        left = identity - rho * outer_product(s, y)
        right = identity - rho * outer_product(y, s)
        hessian_inverse = matmul(matmul(left, hessian_inverse), right) + rho * outer_product(s, s)
      else
        hessian_inverse = identity
      end if

      theta = theta_new
      gradient = gradient_new
      value = value_new
      used = used + 1
      objective_work(used) = value

      if (vector_norm2(s) <= local_options%abs_tol &
        + local_options%rel_tol * max(1.0_dp, vector_norm2(theta))) then
        result%converged = .true.
        exit
      end if
    end do

    result%coefficients = theta(1:p)
    result%shape = exp(theta(q))
    if (used > 0) result%objective = objective_work(1:used)
    result%iterations = used
    if (result%converged) then
      result%status = rpe_success
      result%message = 'converged'
    else if (result%status == rpe_success) then
      result%status = rpe_no_convergence
      result%message = 'maximum iterations reached'
    end if
  end subroutine fit_glm_gamma_mle

  pure function outer_product(x, y) result(matrix)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: matrix(size(x), size(y))
    integer :: j

    do j = 1, size(y)
      matrix(:, j) = x * y(j)
    end do
  end function outer_product

  subroutine fail(result, status, message)
    type(fit_result), intent(inout) :: result
    integer, intent(in) :: status
    character(len=*), intent(in) :: message

    result%status = status
    result%message = message
    result%converged = .false.
  end subroutine fail

end module rpeglmen_gamma_mle
