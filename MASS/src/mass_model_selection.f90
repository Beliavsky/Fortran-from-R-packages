! SPDX-License-Identifier: GPL-3.0-only
module mass_model_selection
  use rrcov_kinds, only : dp
  use mass_types, only : regression_result, model_selection_result, &
    mass_success, mass_invalid_argument, mass_dimension_error
  use mass_regression, only : linear_model_fit
  implicit none
  private
  public :: step_aic_linear, addterm_linear, dropterm_linear
contains

  subroutine step_aic_linear(x, y, result, start, direction, penalty, &
      max_steps, intercept)
    real(dp), intent(in) :: x(:, :), y(:)
    type(model_selection_result), intent(out) :: result
    logical, intent(in), optional :: start(:), intercept
    character(len=*), intent(in), optional :: direction
    real(dp), intent(in), optional :: penalty
    integer, intent(in), optional :: max_steps
    logical, allocatable :: selected(:), candidate(:), best_selected(:)
    real(dp), allocatable :: path(:)
    integer, allocatable :: changed(:)
    real(dp) :: current_aic, best_aic, candidate_aic, kpen
    integer :: p, step, limit, j, best_change, status
    logical :: use_intercept, allow_add, allow_drop
    character(len=12) :: mode

    p = size(x, 2)
    if (size(x, 1) /= size(y) .or. size(y) < 2 .or. p < 1) then
      result%status = mass_dimension_error
      return
    end if
    allocate(selected(p), candidate(p), best_selected(p))
    selected = .false.
    if (present(start)) then
      if (size(start) /= p) then
        result%status = mass_invalid_argument
        return
      end if
      selected = start
    end if
    use_intercept = .true.
    if (present(intercept)) use_intercept = intercept
    mode = "both"
    if (present(direction)) mode = adjustl(direction)
    allow_add = trim(mode) == "both" .or. trim(mode) == "forward"
    allow_drop = trim(mode) == "both" .or. trim(mode) == "backward"
    if (.not. allow_add .and. .not. allow_drop) then
      result%status = mass_invalid_argument
      return
    end if
    kpen = 2.0_dp
    if (present(penalty)) kpen = penalty
    limit = 1000
    if (present(max_steps)) limit = max_steps
    allocate(path(limit + 1), changed(limit))
    call subset_aic(x, y, selected, use_intercept, kpen, current_aic, status)
    if (status /= mass_success) then
      result%status = status
      return
    end if
    path(1) = current_aic
    changed = 0
    do step = 1, limit
      best_aic = current_aic
      best_change = 0
      best_selected = selected
      do j = 1, p
        if ((.not. selected(j) .and. allow_add) .or. &
            (selected(j) .and. allow_drop)) then
          candidate = selected
          candidate(j) = .not. candidate(j)
          call subset_aic(x, y, candidate, use_intercept, kpen, &
            candidate_aic, status)
          if (status == mass_success .and. candidate_aic < best_aic - 1.0e-7_dp) then
            best_aic = candidate_aic
            best_change = merge(j, -j, candidate(j))
            best_selected = candidate
          end if
        end if
      end do
      if (best_change == 0) exit
      selected = best_selected
      current_aic = best_aic
      path(step + 1) = current_aic
      changed(step) = best_change
    end do
    result%selected = selected
    result%steps = step - 1
    result%aic_path = path(1:result%steps + 1)
    result%changed_column = changed(1:result%steps)
    call subset_coefficients(x, y, selected, use_intercept, &
      result%coefficients, status)
    result%status = status
  end subroutine step_aic_linear

  subroutine addterm_linear(x, y, selected, aic, penalty, intercept)
    real(dp), intent(in) :: x(:, :), y(:)
    logical, intent(in) :: selected(:)
    real(dp), allocatable, intent(out) :: aic(:)
    real(dp), intent(in), optional :: penalty
    logical, intent(in), optional :: intercept
    logical, allocatable :: candidate(:)
    real(dp) :: kpen
    integer :: j, status
    logical :: use_intercept
    allocate(aic(size(selected)), candidate(size(selected)))
    aic = huge(1.0_dp)
    kpen = 2.0_dp
    if (present(penalty)) kpen = penalty
    use_intercept = .true.
    if (present(intercept)) use_intercept = intercept
    do j = 1, size(selected)
      if (.not. selected(j)) then
        candidate = selected
        candidate(j) = .true.
        call subset_aic(x, y, candidate, use_intercept, kpen, aic(j), status)
      end if
    end do
  end subroutine addterm_linear

  subroutine dropterm_linear(x, y, selected, aic, penalty, intercept)
    real(dp), intent(in) :: x(:, :), y(:)
    logical, intent(in) :: selected(:)
    real(dp), allocatable, intent(out) :: aic(:)
    real(dp), intent(in), optional :: penalty
    logical, intent(in), optional :: intercept
    logical, allocatable :: candidate(:)
    real(dp) :: kpen
    integer :: j, status
    logical :: use_intercept
    allocate(aic(size(selected)), candidate(size(selected)))
    aic = huge(1.0_dp)
    kpen = 2.0_dp
    if (present(penalty)) kpen = penalty
    use_intercept = .true.
    if (present(intercept)) use_intercept = intercept
    do j = 1, size(selected)
      if (selected(j)) then
        candidate = selected
        candidate(j) = .false.
        call subset_aic(x, y, candidate, use_intercept, kpen, aic(j), status)
      end if
    end do
  end subroutine dropterm_linear

  subroutine subset_aic(x, y, selected, intercept, penalty, aic, status)
    real(dp), intent(in) :: x(:, :), y(:), penalty
    logical, intent(in) :: selected(:), intercept
    real(dp), intent(out) :: aic
    integer, intent(out) :: status
    type(regression_result) :: fit
    real(dp), allocatable :: design(:, :)
    integer :: rank
    call make_design(x, selected, intercept, design)
    if (size(design, 2) == 0) then
      aic = real(size(y), dp) * log(sum(y**2) / real(size(y), dp))
      status = mass_success
      return
    end if
    call linear_model_fit(design, y, fit)
    status = fit%status
    if (status /= mass_success) then
      aic = huge(1.0_dp)
      return
    end if
    rank = fit%rank
    aic = real(size(y), dp) * log(max(sum(fit%residuals**2) / &
      real(size(y), dp), tiny(1.0_dp))) + penalty * real(rank, dp)
  end subroutine subset_aic

  subroutine subset_coefficients(x, y, selected, intercept, coefficients, status)
    real(dp), intent(in) :: x(:, :), y(:)
    logical, intent(in) :: selected(:), intercept
    real(dp), allocatable, intent(out) :: coefficients(:)
    integer, intent(out) :: status
    type(regression_result) :: fit
    real(dp), allocatable :: design(:, :)
    call make_design(x, selected, intercept, design)
    if (size(design, 2) == 0) then
      allocate(coefficients(0))
      status = mass_success
      return
    end if
    call linear_model_fit(design, y, fit)
    coefficients = fit%coefficients
    status = fit%status
  end subroutine subset_coefficients

  subroutine make_design(x, selected, intercept, design)
    real(dp), intent(in) :: x(:, :)
    logical, intent(in) :: selected(:), intercept
    real(dp), allocatable, intent(out) :: design(:, :)
    integer :: n_columns, j, column
    n_columns = count(selected) + merge(1, 0, intercept)
    allocate(design(size(x, 1), n_columns))
    column = 0
    if (intercept) then
      column = 1
      design(:, column) = 1.0_dp
    end if
    do j = 1, size(selected)
      if (selected(j)) then
        column = column + 1
        design(:, column) = x(:, j)
      end if
    end do
  end subroutine make_design

end module mass_model_selection
