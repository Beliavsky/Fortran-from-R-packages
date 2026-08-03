! SPDX-License-Identifier: GPL-3.0-only
module rumidas_weights
  use rumidas_kinds, only: dp
  use rumidas_status
  use rumidas_types, only: RUMIDAS_BETA_LAG, RUMIDAS_ALMON_LAG
  implicit none
  private
  public :: beta_weights, exponential_almon_weights, rumidas_lag_weights
  public :: beta_function, exp_almon, mv_into_mat
  public :: midas_weighted_component, lag_matrix_from_period_index

contains

  function beta_weights(k, w1, w2, status) result(weights)
    integer, intent(in) :: k
    real(dp), intent(in) :: w1, w2
    integer, intent(out), optional :: status
    real(dp), allocatable :: weights(:)
    real(dp), allocatable :: logw(:)
    real(dp) :: x, maximum, denominator
    integer :: j, stat

    stat = RUMIDAS_SUCCESS
    allocate(weights(max(k, 0)))
    if (k <= 0 .or. w1 <= 0.0_dp .or. w2 <= 0.0_dp) then
      if (k > 0) weights = 0.0_dp
      stat = RUMIDAS_INVALID_PARAMETER
      if (present(status)) status = stat
      return
    end if

    allocate(logw(k))
    do j = 1, k
      x = real(j, dp) / real(k, dp)
      if (x >= 1.0_dp) then
        if (w2 > 1.0_dp) then
          logw(j) = -huge(1.0_dp)
        else if (abs(w2 - 1.0_dp) <= 10.0_dp * epsilon(1.0_dp)) then
          logw(j) = (w1 - 1.0_dp) * log(x)
        else
          logw(j) = huge(1.0_dp)
        end if
      else
        logw(j) = (w1 - 1.0_dp) * log(max(x, tiny(1.0_dp))) + &
          (w2 - 1.0_dp) * log(max(1.0_dp - x, tiny(1.0_dp)))
      end if
    end do
    maximum = maxval(logw, mask=logw < huge(1.0_dp) / 2.0_dp)
    weights = 0.0_dp
    do j = 1, k
      if (logw(j) < huge(1.0_dp) / 2.0_dp .and. logw(j) > -huge(1.0_dp) / 2.0_dp) &
        weights(j) = exp(logw(j) - maximum)
    end do
    denominator = sum(weights)
    if (denominator <= tiny(1.0_dp)) then
      weights = 0.0_dp
      stat = RUMIDAS_NUMERICAL_ERROR
    else
      weights = weights / denominator
    end if
    if (present(status)) status = stat
  end function beta_weights

  function exponential_almon_weights(k, w1, w2, status) result(weights)
    integer, intent(in) :: k
    real(dp), intent(in) :: w1, w2
    integer, intent(out), optional :: status
    real(dp), allocatable :: weights(:)
    real(dp), allocatable :: exponent(:)
    real(dp) :: maximum, denominator, x
    integer :: j, stat

    stat = RUMIDAS_SUCCESS
    allocate(weights(max(k, 0)))
    if (k <= 0) then
      stat = RUMIDAS_INVALID_PARAMETER
      if (present(status)) status = stat
      return
    end if
    allocate(exponent(k))
    do j = 1, k
      x = real(j, dp)
      exponent(j) = w1 * x + w2 * x * x
    end do
    maximum = maxval(exponent)
    weights = exp(exponent - maximum)
    denominator = sum(weights)
    if (denominator <= tiny(1.0_dp)) then
      weights = 0.0_dp
      stat = RUMIDAS_NUMERICAL_ERROR
    else
      weights = weights / denominator
    end if
    if (present(status)) status = stat
  end function exponential_almon_weights

  function rumidas_lag_weights(k, w2, lag_function, status) result(weights)
    integer, intent(in) :: k, lag_function
    real(dp), intent(in) :: w2
    integer, intent(out), optional :: status
    real(dp), allocatable :: weights(:)
    real(dp) :: raw(k + 1)
    real(dp) :: w1
    integer :: j, stat

    stat = RUMIDAS_SUCCESS
    if (k <= 0) then
      allocate(weights(0))
      stat = RUMIDAS_INVALID_PARAMETER
      if (present(status)) status = stat
      return
    end if
    w1 = merge(1.0_dp, 0.0_dp, lag_function == RUMIDAS_BETA_LAG)
    select case (lag_function)
    case (RUMIDAS_BETA_LAG)
      raw = beta_weights(k + 1, w1, w2, stat)
    case (RUMIDAS_ALMON_LAG)
      raw = exponential_almon_weights(k + 1, w1, w2, stat)
    case default
      allocate(weights(k + 1))
      weights = 0.0_dp
      stat = RUMIDAS_INVALID_INPUT
      if (present(status)) status = stat
      return
    end select
    allocate(weights(k + 1))
    weights(k + 1) = 0.0_dp
    do j = 1, k
      weights(j) = raw(k + 1 - j)
    end do
    if (present(status)) status = stat
  end function rumidas_lag_weights

  subroutine midas_weighted_component(mv_matrix, k, w2, lag_function, component, status, split_sign)
    real(dp), intent(in) :: mv_matrix(:, :)
    integer, intent(in) :: k, lag_function
    real(dp), intent(in) :: w2
    real(dp), intent(out) :: component(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: split_sign
    real(dp), allocatable :: weights(:)
    real(dp) :: value
    integer :: t, j, sign_mode

    status = RUMIDAS_SUCCESS
    sign_mode = 0
    if (present(split_sign)) sign_mode = split_sign
    if (size(mv_matrix, 1) /= k + 1 .or. size(component) /= size(mv_matrix, 2)) then
      component = 0.0_dp
      status = RUMIDAS_DIMENSION_ERROR
      return
    end if
    weights = rumidas_lag_weights(k, w2, lag_function, status)
    if (status /= RUMIDAS_SUCCESS) then
      component = 0.0_dp
      return
    end if
    do t = 1, size(mv_matrix, 2)
      value = 0.0_dp
      do j = 1, k + 1
        select case (sign_mode)
        case (1)
          value = value + weights(j) * max(mv_matrix(j, t), 0.0_dp)
        case (-1)
          value = value + weights(j) * min(mv_matrix(j, t), 0.0_dp)
        case default
          value = value + weights(j) * mv_matrix(j, t)
        end select
      end do
      component(t) = value
    end do
  end subroutine midas_weighted_component

  subroutine lag_matrix_from_period_index(mv, current_period_index, k, matrix, status)
    real(dp), intent(in) :: mv(:)
    integer, intent(in) :: current_period_index(:)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: matrix(:, :)
    integer, intent(out) :: status
    integer :: t, first_index, current

    status = RUMIDAS_SUCCESS
    allocate(matrix(k + 1, size(current_period_index)))
    matrix = 0.0_dp
    if (k < 0) then
      status = RUMIDAS_INVALID_PARAMETER
      return
    end if
    do t = 1, size(current_period_index)
      current = current_period_index(t)
      first_index = current - k
      if (first_index < 1 .or. current > size(mv)) then
        status = RUMIDAS_DIMENSION_ERROR
        matrix = 0.0_dp
        return
      end if
      matrix(:, t) = mv(first_index:current)
    end do
  end subroutine lag_matrix_from_period_index

  function beta_function(k_values, k, w1, w2, status) result(values)
    integer, intent(in) :: k_values(:), k
    real(dp), intent(in) :: w1, w2
    integer, intent(out), optional :: status
    real(dp), allocatable :: values(:), all_weights(:)
    integer :: i, stat
    all_weights = beta_weights(k, w1, w2, stat)
    allocate(values(size(k_values)))
    values = 0.0_dp
    if (stat == RUMIDAS_SUCCESS) then
      do i = 1, size(k_values)
        if (k_values(i) < 1 .or. k_values(i) > k) then
          stat = RUMIDAS_INVALID_INPUT
          exit
        end if
        values(i) = all_weights(k_values(i))
      end do
    end if
    if (present(status)) status = stat
  end function beta_function

  function exp_almon(k_values, k, w1, w2, status) result(values)
    integer, intent(in) :: k_values(:), k
    real(dp), intent(in) :: w1, w2
    integer, intent(out), optional :: status
    real(dp), allocatable :: values(:), all_weights(:)
    integer :: i, stat
    all_weights = exponential_almon_weights(k, w1, w2, stat)
    allocate(values(size(k_values)))
    values = 0.0_dp
    if (stat == RUMIDAS_SUCCESS) then
      do i = 1, size(k_values)
        if (k_values(i) < 1 .or. k_values(i) > k) then
          stat = RUMIDAS_INVALID_INPUT
          exit
        end if
        values(i) = all_weights(k_values(i))
      end do
    end if
    if (present(status)) status = stat
  end function exp_almon

  subroutine mv_into_mat(mv, current_period_index, k, matrix, status)
    real(dp), intent(in) :: mv(:)
    integer, intent(in) :: current_period_index(:), k
    real(dp), allocatable, intent(out) :: matrix(:, :)
    integer, intent(out) :: status
    call lag_matrix_from_period_index(mv, current_period_index, k, matrix, status)
  end subroutine mv_into_mat

end module rumidas_weights
