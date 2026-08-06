! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the RPEGLMEN computational code.

module rpeglmen_cv
  use, intrinsic :: iso_fortran_env, only : int64
  use rpeglmen_kinds, only : dp
  use rpeglmen_types, only : enet_options, fit_result, path_result, clear_fit_result, &
    rpe_success, rpe_invalid_input, rpe_no_convergence, model_exponential, model_gamma
  use rpeglmen_math, only : all_finite
  use rpeglmen_likelihood, only : exp_negative_log_likelihood, &
    grad_exp_negative_log_likelihood, gamma_negative_log_likelihood, &
    grad_gamma_negative_log_likelihood, predict_mean
  use rpeglmen_solver, only : fit_fixed_model
  use rpeglmen_gamma_mle, only : fit_glm_gamma_mle
  implicit none
  private

  public :: compute_lambda_max, generate_lambda_grid
  public :: fit_regularization_path, glmnet_exp, fit_glm_gamma_net

contains

  real(dp) function compute_lambda_max(a, b, model, shape, options) result(lambda_max)
    real(dp), intent(in) :: a(:, :), b(:), shape
    integer, intent(in) :: model
    type(enet_options), intent(in) :: options
    real(dp) :: coefficients(size(a, 2)), gradient(size(a, 2)), divisor
    integer :: first_penalized

    coefficients = 0.0_dp
    select case (model)
    case (model_exponential)
      gradient = grad_exp_negative_log_likelihood(coefficients, a, b)
    case (model_gamma)
      gradient = grad_gamma_negative_log_likelihood(coefficients, a, b, shape)
    case default
      lambda_max = 0.0_dp
      return
    end select

    first_penalized = 1
    if (options%has_intercept .and. .not. options%penalize_intercept .and. size(gradient) > 0) first_penalized = 2
    if (first_penalized > size(gradient)) then
      lambda_max = 1.0_dp
      return
    end if

    divisor = max(options%alpha, 1.0e-6_dp)
    lambda_max = maxval(abs(gradient(first_penalized:))) / divisor
    lambda_max = max(lambda_max, 100.0_dp * epsilon(1.0_dp))
  end function compute_lambda_max

  function generate_lambda_grid(a, b, model, shape, options) result(grid)
    real(dp), intent(in) :: a(:, :), b(:), shape
    integer, intent(in) :: model
    type(enet_options), intent(in) :: options
    real(dp), allocatable :: grid(:)
    real(dp) :: lambda_max, lambda_min, fraction
    integer :: i, nlambda

    nlambda = max(1, options%num_lambda)
    allocate(grid(nlambda))
    lambda_max = compute_lambda_max(a, b, model, shape, options)
    if (options%min_lambda_absolute > 0.0_dp) then
      lambda_min = min(lambda_max, options%min_lambda_absolute)
    else
      lambda_min = lambda_max * max(tiny(1.0_dp), min(1.0_dp, options%min_lambda_ratio))
    end if

    if (nlambda == 1 .or. lambda_max <= lambda_min) then
      grid = lambda_max
      return
    end if

    do i = 1, nlambda
      fraction = real(i - 1, dp) / real(nlambda - 1, dp)
      grid(i) = exp(log(lambda_max) + fraction * (log(lambda_min) - log(lambda_max)))
    end do
  end function generate_lambda_grid

  subroutine fit_regularization_path(a, b, model, shape, path, options, lambda_grid)
    real(dp), intent(in) :: a(:, :), b(:), shape
    integer, intent(in) :: model
    type(path_result), intent(out) :: path
    type(enet_options), intent(in), optional :: options
    real(dp), intent(in), optional :: lambda_grid(:)

    type(enet_options) :: local_options
    type(fit_result) :: fit
    real(dp), allocatable :: grid(:), start(:)
    integer :: i, p

    local_options = enet_options()
    if (present(options)) local_options = options
    p = size(a, 2)

    if (present(lambda_grid)) then
      if (size(lambda_grid) < 1 .or. any(lambda_grid < 0.0_dp)) then
        path%status = rpe_invalid_input
        path%message = 'lambda grid must be nonempty and nonnegative'
        return
      end if
      allocate(grid(size(lambda_grid)))
      grid = lambda_grid
    else
      allocate(grid(max(1, local_options%num_lambda)))
      grid = generate_lambda_grid(a, b, model, shape, local_options)
    end if

    allocate(path%coefficients(p, size(grid)), path%iterations(size(grid)), path%converged(size(grid)))
    path%lambda_grid = grid
    allocate(start(p))
    start = 0.0_dp

    do i = 1, size(grid)
      call fit_fixed_model(a, b, model, shape, grid(i), local_options, fit, start)
      if (.not. allocated(fit%coefficients)) then
        path%status = fit%status
        path%message = fit%message
        return
      end if
      path%coefficients(:, i) = fit%coefficients
      path%iterations(i) = fit%iterations
      path%converged(i) = fit%converged
      start = fit%coefficients
    end do

    path%status = rpe_success
    path%message = 'completed'
  end subroutine fit_regularization_path

  subroutine glmnet_exp(a, b, result, options)
    real(dp), intent(in) :: a(:, :), b(:)
    type(fit_result), intent(out) :: result
    type(enet_options), intent(in), optional :: options
    type(enet_options) :: local_options

    local_options = enet_options()
    if (present(options)) local_options = options
    call cross_validate(a, b, model_exponential, 1.0_dp, result, local_options)
  end subroutine glmnet_exp

  subroutine fit_glm_gamma_net(a, b, result, options)
    real(dp), intent(in) :: a(:, :), b(:)
    type(fit_result), intent(out) :: result
    type(enet_options), intent(in), optional :: options
    type(enet_options) :: local_options
    type(fit_result) :: mle

    local_options = enet_options()
    local_options%penalize_intercept = .true.
    if (present(options)) local_options = options

    call fit_glm_gamma_mle(a, b, mle, local_options)
    if (.not. allocated(mle%coefficients)) then
      result = mle
      return
    end if
    call cross_validate(a, b, model_gamma, mle%shape, result, local_options)
    result%shape = mle%shape
  end subroutine fit_glm_gamma_net

  subroutine cross_validate(a, b, model, shape, result, options)
    real(dp), intent(in) :: a(:, :), b(:), shape
    integer, intent(in) :: model
    type(fit_result), intent(out) :: result
    type(enet_options), intent(in) :: options

    real(dp), allocatable :: grid(:), sum_error(:), sumsq_error(:), errors(:)
    real(dp), allocatable :: a_train(:, :), a_test(:, :), b_train(:), b_test(:)
    integer, allocatable :: permutation(:), fold_id(:)
    logical, allocatable :: train_mask(:), test_mask(:)
    type(fit_result) :: fit
    integer :: n, p, folds, repeats, repeat_index, fold, i, count, best
    real(dp) :: error

    call clear_fit_result(result)
    n = size(b)
    p = size(a, 2)
    if (size(a, 1) /= n .or. n < 2 .or. p < 1 .or. any(b <= 0.0_dp) .or. .not. all_finite(b)) then
      result%status = rpe_invalid_input
      result%message = 'invalid design or positive response vector'
      return
    end if

    folds = min(max(2, options%k_fold), n)
    repeats = max(1, options%k_fold_iter)
    grid = generate_lambda_grid(a, b, model, shape, options)
    allocate(sum_error(size(grid)), sumsq_error(size(grid)), errors(size(grid)))
    allocate(permutation(n), fold_id(n), train_mask(n), test_mask(n))
    sum_error = 0.0_dp
    sumsq_error = 0.0_dp
    count = 0

    do repeat_index = 1, repeats
      call make_permutation(n, options%seed + 104729 * (repeat_index - 1), permutation)
      do i = 1, n
        fold_id(permutation(i)) = 1 + modulo(i - 1, folds)
      end do

      do fold = 1, folds
        test_mask = fold_id == fold
        train_mask = .not. test_mask
        call subset_rows(a, b, train_mask, a_train, b_train)
        call subset_rows(a, b, test_mask, a_test, b_test)

        do i = 1, size(grid)
          call fit_fixed_model(a_train, b_train, model, shape, grid(i), options, fit)
          if (.not. allocated(fit%coefficients)) then
            error = huge(1.0_dp) / 1000.0_dp
          else
            error = validation_error(a_test, b_test, fit%coefficients, model, shape, options%cv_metric)
          end if
          errors(i) = error
        end do
        sum_error = sum_error + errors
        sumsq_error = sumsq_error + errors * errors
        count = count + 1
      end do
    end do

    result%lambda_grid = grid
    result%cv_mean = sum_error / real(count, dp)
    allocate(result%cv_sd(size(grid)))
    if (count > 1) then
      result%cv_sd = sqrt(max(0.0_dp, (sumsq_error - real(count, dp) * result%cv_mean**2) / real(count - 1, dp)))
    else
      result%cv_sd = 0.0_dp
    end if

    best = minloc(result%cv_mean, dim=1)
    result%selected_lambda = grid(best)
    call fit_fixed_model(a, b, model, shape, result%selected_lambda, options, fit)
    if (allocated(fit%coefficients)) result%coefficients = fit%coefficients
    if (allocated(fit%objective)) result%objective = fit%objective
    result%iterations = fit%iterations
    result%converged = fit%converged
    result%status = fit%status
    result%message = fit%message
    result%shape = shape
  end subroutine cross_validate

  real(dp) function validation_error(a, b, coefficients, model, shape, metric) result(error)
    real(dp), intent(in) :: a(:, :), b(:), coefficients(:), shape
    integer, intent(in) :: model
    character(len=*), intent(in) :: metric
    real(dp) :: predicted(size(b))
    character(len=:), allocatable :: metric_lower

    metric_lower = lower_case(trim(metric))
    select case (metric_lower)
    case ('rmse', 'source', 'prediction')
      predicted = predict_mean(a, coefficients)
      error = sqrt(sum((predicted - b)**2) / real(size(b), dp))
    case default
      if (model == model_exponential) then
        error = exp_negative_log_likelihood(coefficients, a, b) / real(size(b), dp)
      else
        error = gamma_negative_log_likelihood(coefficients, a, b, shape)
      end if
    end select
  end function validation_error

  subroutine subset_rows(a, b, mask, a_sub, b_sub)
    real(dp), intent(in) :: a(:, :), b(:)
    logical, intent(in) :: mask(:)
    real(dp), allocatable, intent(out) :: a_sub(:, :), b_sub(:)
    integer :: i, j, row_count

    row_count = count(mask)
    allocate(a_sub(row_count, size(a, 2)), b_sub(row_count))
    j = 0
    do i = 1, size(mask)
      if (mask(i)) then
        j = j + 1
        a_sub(j, :) = a(i, :)
        b_sub(j) = b(i)
      end if
    end do
  end subroutine subset_rows

  subroutine make_permutation(n, seed, permutation)
    integer, intent(in) :: n, seed
    integer, intent(out) :: permutation(n)
    integer(int64) :: state
    integer :: i, j, temporary

    do i = 1, n
      permutation(i) = i
    end do
    state = int(max(1, abs(seed)), int64)
    do i = n, 2, -1
      state = modulo(16807_int64 * state, 2147483647_int64)
      j = 1 + int(modulo(state, int(i, int64)))
      temporary = permutation(i)
      permutation(i) = permutation(j)
      permutation(j) = temporary
    end do
  end subroutine make_permutation

  pure function lower_case(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: i, code

    lower = text
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
    end do
  end function lower_case

end module rpeglmen_cv
