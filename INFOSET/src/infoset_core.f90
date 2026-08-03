! SPDX-License-Identifier: GPL-2.0-or-later
module infoset_core
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use infoset_kinds, only : dp
  use infoset_status
  use infoset_types
  use infoset_stats, only : sort_real, median_real, quantile_real, left_histogram_risk
  use infoset_mixture, only : tail_mixture
  implicit none
  private
  public :: g_ret, create_overlapping_windows, infoset_estimate, lr_cp

  interface g_ret
    module procedure gross_returns_vector
    module procedure gross_returns_matrix
  end interface g_ret
contains
  function gross_returns_vector(prices) result(returns)
    real(dp), intent(in) :: prices(:)
    real(dp), allocatable :: returns(:)
    integer :: n
    n = size(prices)
    if (n < 2 .or. any(prices <= 0.0_dp) .or. .not. all(ieee_is_finite(prices))) then
      allocate(returns(0))
      return
    end if
    allocate(returns(n - 1))
    returns = prices(2:n) / prices(1:n - 1)
    call sort_real(returns)
  end function gross_returns_vector

  function gross_returns_matrix(prices) result(returns)
    real(dp), intent(in) :: prices(:,:)
    real(dp), allocatable :: returns(:,:)
    real(dp), allocatable :: one_asset(:)
    integer :: j, n, p
    n = size(prices, 1)
    p = size(prices, 2)
    if (n < 2 .or. p < 1 .or. any(prices <= 0.0_dp) &
        .or. .not. all(ieee_is_finite(prices))) then
      allocate(returns(0, 0))
      return
    end if
    allocate(returns(n - 1, p))
    do j = 1, p
      one_asset = gross_returns_vector(prices(:, j))
      returns(:, j) = one_asset
    end do
  end function gross_returns_matrix

  function create_overlapping_windows(data, window_size, overlap) result(windows)
    real(dp), intent(in) :: data(:,:)
    integer, intent(in) :: window_size, overlap
    type(window_collection) :: windows
    integer :: number, t, first
    windows%window_size = window_size
    windows%step = overlap
    if (window_size <= 0 .or. overlap <= 0 .or. size(data, 1) < window_size &
        .or. size(data, 2) < 1 .or. .not. all(ieee_is_finite(data))) then
      allocate(windows%values(0, 0, 0))
      windows%status = infoset_invalid_argument
      return
    end if
    number = (size(data, 1) - window_size) / overlap + 1
    allocate(windows%values(window_size, size(data, 2), number))
    do t = 1, number
      first = 1 + (t - 1) * overlap
      windows%values(:, :, t) = data(first:first + window_size - 1, :)
    end do
    windows%status = infoset_success
  end function create_overlapping_windows

  subroutine infoset_estimate(y, result, control)
    real(dp), intent(in) :: y(:)
    type(information_set_result), intent(out) :: result
    type(mixture_control), intent(in), optional :: control
    type(tail_mixture_result) :: fit
    type(mixture_control) :: ctl
    real(dp) :: shift, sample_median
    integer :: iteration
    result = information_set_result()
    ctl = mixture_control()
    if (present(control)) ctl = control
    if (size(y) < 6 .or. any(y <= 0.0_dp) .or. .not. all(ieee_is_finite(y))) then
      result%status = infoset_invalid_argument
      return
    end if
    shift = 0.0_dp
    sample_median = median_real(y)
    do iteration = 1, 2
      call tail_mixture(y, shift, iteration, fit, ctl)
      if (fit%status /= infoset_success .and. fit%status /= infoset_not_converged) then
        if (result%n_change_points == 0) result%status = fit%status
        exit
      end if
      if (fit%flag /= 0 .or. fit%change_point <= shift) exit
      if (fit%change_point >= sample_median) then
        if (result%n_change_points == 0) result%status = infoset_no_split
        exit
      end if
      result%n_change_points = result%n_change_points + 1
      result%change_points(result%n_change_points) = fit%change_point
      result%prior_probability(result%n_change_points) = fit%left_probability
      result%first_type_error(result%n_change_points) = fit%first_type_error
      result%second_type_error(result%n_change_points) = fit%second_type_error
      result%left_mean(result%n_change_points) = fit%left_mean
      result%left_sd(result%n_change_points) = fit%left_sd
      if (fit%status == infoset_not_converged) result%status = infoset_not_converged
      shift = fit%change_point
    end do
    if (result%n_change_points > 0 .and. result%status == infoset_success) then
      result%status = infoset_success
    else if (result%n_change_points == 0 .and. result%status == infoset_success) then
      result%status = infoset_no_split
    end if
  end subroutine infoset_estimate

  subroutine lr_cp(data, window_size, overlap, result, control)
    real(dp), intent(in) :: data(:,:)
    integer, intent(in) :: window_size, overlap
    type(left_risk_result), intent(out) :: result
    type(mixture_control), intent(in), optional :: control
    type(window_collection) :: windows
    type(information_set_result) :: fit
    type(mixture_control) :: ctl
    real(dp), allocatable :: gross(:), log_returns(:)
    real(dp) :: risk, change_point
    integer :: p, number, i, t, status
    logical :: used_fallback

    result = left_risk_result()
    ctl = mixture_control()
    if (present(control)) ctl = control
    windows = create_overlapping_windows(data, window_size, overlap)
    if (windows%status /= infoset_success) then
      allocate(result%values(0, 0), result%first_change_point(0))
      result%status = windows%status
      return
    end if
    p = size(data, 2)
    number = size(windows%values, 3)
    allocate(result%values(p, number), result%first_change_point(p))
    result%values = 0.0_dp
    result%first_change_point = 0.0_dp
    result%status = infoset_success
    used_fallback = .false.
    do i = 1, p
      gross = gross_returns_vector(data(:, i))
      call infoset_estimate(gross, fit, ctl)
      if (fit%n_change_points >= 1) then
        change_point = fit%change_points(1)
      else
        change_point = quantile_real(gross, 0.10_dp)
        used_fallback = .true.
      end if
      result%first_change_point(i) = change_point
      do t = 1, number
        allocate(log_returns(window_size - 1))
        log_returns = log(windows%values(2:window_size, i, t) &
          / windows%values(1:window_size - 1, i, t))
        call left_histogram_risk(log_returns, log(change_point), risk, status)
        if (status /= infoset_success) then
          risk = -sum(log_returns) / real(size(log_returns), dp)
          used_fallback = .true.
        end if
        result%values(i, t) = risk
        deallocate(log_returns)
      end do
    end do
    if (used_fallback) result%status = infoset_no_split
  end subroutine lr_cp
end module infoset_core
