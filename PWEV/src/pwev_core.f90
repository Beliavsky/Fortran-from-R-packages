! SPDX-License-Identifier: GPL-3.0-only
module pwev_core
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use pwev_kinds, only : dp
  use pwev_status
  use pwev_types
  use pwev_pso, only : pso_ensemble_weights
  use pwev_metrics, only : pwev_accuracy_table
  use pwev_models, only : fit_pwev_base_models
  implicit none
  private
  public :: pwev_fit, pwev_fit_from_components
contains

  subroutine pwev_fit(data, split_ratio, result, status, control)
    real(dp), intent(in) :: data(:)
    real(dp), intent(in) :: split_ratio
    type(pwev_result), intent(out) :: result
    integer, intent(out), optional :: status
    type(pwev_control), intent(in), optional :: control
    type(pwev_control) :: local_control
    real(dp), allocatable :: train_models(:, :), test_models(:, :)
    integer :: n, n_train, n_test
    integer :: base_status(PWEV_N_BASE_MODELS)

    local_control = pwev_control()
    if (present(control)) local_control = control
    n = size(data)
    if (n < 8 .or. split_ratio <= 0.0_dp .or. split_ratio >= 1.0_dp .or. &
        any(.not. ieee_is_finite(data))) then
      call set_failure(result, PWEV_INVALID_INPUT, 'data must be finite with at least 8 observations and 0 < split_ratio < 1')
      if (present(status)) status = result%status
      return
    end if
    n_train = int(real(n, dp) * split_ratio)
    n_train = max(4, min(n - 1, n_train))
    n_test = n - n_train
    allocate(train_models(n_train, PWEV_N_BASE_MODELS), test_models(n_test, PWEV_N_BASE_MODELS))
    call fit_pwev_base_models(data(1:n_train), data(n_train + 1:n), local_control, &
      train_models, test_models, base_status)
    call pwev_fit_from_components(data(1:n_train), data(n_train + 1:n), train_models, test_models, &
      result, status, local_control, model_status=base_status)
  end subroutine pwev_fit

  subroutine pwev_fit_from_components(train_actual, test_actual, train_models, test_models, result, &
      status, control, model_status)
    real(dp), intent(in) :: train_actual(:), test_actual(:)
    real(dp), intent(in) :: train_models(:, :), test_models(:, :)
    type(pwev_result), intent(out) :: result
    integer, intent(out), optional :: status
    type(pwev_control), intent(in), optional :: control
    integer, intent(in), optional :: model_status(PWEV_N_BASE_MODELS)
    type(pwev_control) :: local_control
    type(pwev_pso_result) :: pso
    real(dp), allocatable :: train_all(:, :), test_all(:, :)
    logical :: base_failed

    local_control = pwev_control()
    if (present(control)) local_control = control
    if (size(train_actual) <= 0 .or. size(test_actual) <= 0 .or. &
        size(train_models, 1) /= size(train_actual) .or. size(test_models, 1) /= size(test_actual) .or. &
        size(train_models, 2) /= PWEV_N_BASE_MODELS .or. &
        size(test_models, 2) /= PWEV_N_BASE_MODELS .or. &
        any(.not. ieee_is_finite(train_actual)) .or. any(.not. ieee_is_finite(test_actual)) .or. &
        any(.not. ieee_is_finite(train_models)) .or. any(.not. ieee_is_finite(test_models))) then
      call set_failure(result, PWEV_INVALID_INPUT, 'component matrices have invalid dimensions or values')
      if (present(status)) status = result%status
      return
    end if

    result%model_status = PWEV_SUCCESS
    if (present(model_status)) result%model_status = model_status
    base_failed = any(result%model_status /= PWEV_SUCCESS)
    if (base_failed .and. local_control%fail_on_base_model_error) then
      call set_failure(result, PWEV_MODEL_FAILURE, 'a base model failed and fail_on_base_model_error is enabled')
      if (present(status)) status = result%status
      return
    end if

    call pso_ensemble_weights(train_actual, train_models, local_control, pso)
    if (pso%status /= PWEV_SUCCESS) then
      call set_failure(result, PWEV_OPTIMIZER_FAILURE, 'particle-swarm ensemble optimization failed')
      if (present(status)) status = result%status
      return
    end if

    allocate(train_all(size(train_actual), PWEV_N_MODELS), test_all(size(test_actual), PWEV_N_MODELS))
    train_all(:, 1:PWEV_N_BASE_MODELS) = train_models
    test_all(:, 1:PWEV_N_BASE_MODELS) = test_models
    train_all(:, PWEV_N_MODELS) = matmul(train_models, pso%weights)
    test_all(:, PWEV_N_MODELS) = matmul(test_models, pso%weights)

    allocate(result%train_fitted(size(train_actual), PWEV_N_MODELS + 1))
    allocate(result%test_pred(size(test_actual), PWEV_N_MODELS + 1))
    result%train_fitted(:, 1) = train_actual
    result%train_fitted(:, 2:) = train_all
    result%test_pred(:, 1) = test_actual
    result%test_pred(:, 2:) = test_all
    result%weights = pso%weights
    result%train_size = size(train_actual)
    result%test_size = size(test_actual)
    result%pso_iterations = pso%iterations
    result%pso_objective = pso%objective
    call pwev_accuracy_table(train_actual, train_all, test_actual, test_all, result%accuracy, &
      local_control%round_accuracy)
    if (base_failed) then
      result%status = PWEV_MODEL_FAILURE
      result%message = 'completed with one or more base-model fallbacks'
    else
      result%status = PWEV_SUCCESS
      result%message = 'success'
    end if
    if (present(status)) status = result%status
  end subroutine pwev_fit_from_components

  subroutine set_failure(result, code, message)
    type(pwev_result), intent(out) :: result
    integer, intent(in) :: code
    character(len=*), intent(in) :: message
    result%status = code
    result%message = message
    allocate(result%train_fitted(0, 0), result%test_pred(0, 0), result%accuracy(0, 0), result%weights(0))
  end subroutine set_failure

end module pwev_core
