! SPDX-License-Identifier: MIT
! Based on bidask 2.1.5, Copyright (c) 2024 Emanuele Guidotti.
module bidask_windows
  use bidask_kinds, only: dp
  use bidask_types, only: ohlc_data, spread_result, spread_series_result
  use bidask_statistics, only: nan_dp
  use bidask_estimators, only: edge_estimate, estimate_method
  implicit none
  private

  public :: edge_rolling, edge_expanding
  public :: spread, spread_expanding, spread_endpoints

  interface edge_rolling
    module procedure edge_rolling_fixed
    module procedure edge_rolling_adaptive
  end interface edge_rolling

  interface spread
    module procedure spread_full
    module procedure spread_fixed
    module procedure spread_adaptive
  end interface spread

contains

  function edge_rolling_fixed(open, high, low, close, width, signed, na_rm) result(values)
    real(dp), intent(in) :: open(:), high(:), low(:), close(:)
    integer, intent(in) :: width
    logical, intent(in), optional :: signed, na_rm
    real(dp), allocatable :: values(:)
    logical :: keep_sign, remove_na
    integer :: n, i, first

    keep_sign = .false.
    remove_na = .false.
    if (present(signed)) keep_sign = signed
    if (present(na_rm)) remove_na = na_rm
    n = size(open)
    allocate(values(n), source=nan_dp())
    if (size(high) /= n .or. size(low) /= n .or. size(close) /= n) return
    if (width < 3 .or. width > n) return
    do i = width, n
      first = i - width + 1
      values(i) = edge_estimate(open(first:i), high(first:i), low(first:i), &
        close(first:i), keep_sign, remove_na)
    end do
  end function edge_rolling_fixed

  function edge_rolling_adaptive(open, high, low, close, width, signed, na_rm) result(values)
    real(dp), intent(in) :: open(:), high(:), low(:), close(:)
    integer, intent(in) :: width(:)
    logical, intent(in), optional :: signed, na_rm
    real(dp), allocatable :: values(:)
    logical :: keep_sign, remove_na
    integer :: n, i, first, last

    keep_sign = .false.
    remove_na = .false.
    if (present(signed)) keep_sign = signed
    if (present(na_rm)) remove_na = na_rm
    n = size(open)
    allocate(values(n), source=nan_dp())
    if (size(high) /= n .or. size(low) /= n .or. size(close) /= n) return
    if (size(width) == n) then
      do i = 1, n
        if (width(i) >= 3 .and. width(i) <= i) then
          first = i - width(i) + 1
          values(i) = edge_estimate(open(first:i), high(first:i), low(first:i), &
            close(first:i), keep_sign, remove_na)
        end if
      end do
    else if (size(width) >= 2) then
      do i = 2, size(width)
        first = max(1, width(i - 1))
        last = width(i)
        if (last >= first .and. last <= n .and. last - first + 1 >= 3) then
          values(last) = edge_estimate(open(first:last), high(first:last), &
            low(first:last), close(first:last), keep_sign, remove_na)
        end if
      end do
    end if
  end function edge_rolling_adaptive

  function edge_expanding(open, high, low, close, signed, na_rm) result(values)
    real(dp), intent(in) :: open(:), high(:), low(:), close(:)
    logical, intent(in), optional :: signed, na_rm
    real(dp), allocatable :: values(:)
    logical :: keep_sign, remove_na
    integer :: n, i

    keep_sign = .false.
    remove_na = .true.
    if (present(signed)) keep_sign = signed
    if (present(na_rm)) remove_na = na_rm
    n = size(open)
    allocate(values(n), source=nan_dp())
    if (size(high) /= n .or. size(low) /= n .or. size(close) /= n) return
    do i = 3, n
      values(i) = edge_estimate(open(1:i), high(1:i), low(1:i), close(1:i), &
        keep_sign, remove_na)
    end do
  end function edge_expanding

  function spread_full(data, methods, signed, na_rm) result(result)
    type(ohlc_data), intent(in) :: data
    character(len=*), intent(in) :: methods(:)
    logical, intent(in), optional :: signed, na_rm
    type(spread_result) :: result
    logical :: keep_sign, remove_na
    integer :: j

    keep_sign = .false.
    remove_na = .false.
    if (present(signed)) keep_sign = signed
    if (present(na_rm)) remove_na = na_rm
    allocate(result%method(size(methods)), result%value(size(methods)))
    result%method = methods
    result%value = nan_dp()
    if (.not. data%valid_dimensions()) then
      result%ok = .false.
      result%message = 'OHLC arrays must be allocated and have equal lengths'
      return
    end if
    do j = 1, size(methods)
      result%value(j) = estimate_method(data%open, data%high, data%low, data%close, &
        trim(methods(j)), keep_sign, remove_na)
    end do
  end function spread_full

  function spread_fixed(data, width, methods, signed, na_rm) result(result)
    type(ohlc_data), intent(in) :: data
    integer, intent(in) :: width
    character(len=*), intent(in) :: methods(:)
    logical, intent(in), optional :: signed, na_rm
    type(spread_series_result) :: result
    logical :: keep_sign, remove_na
    integer :: n, i, j, first

    keep_sign = .false.
    remove_na = .false.
    if (present(signed)) keep_sign = signed
    if (present(na_rm)) remove_na = na_rm
    n = data%size()
    allocate(result%method(size(methods)), result%value(n, size(methods)))
    result%method = methods
    result%value = nan_dp()
    if (.not. data%valid_dimensions()) then
      result%ok = .false.
      result%message = 'OHLC arrays must be allocated and have equal lengths'
      return
    end if
    if (width < 3 .or. width > n) return
    do i = width, n
      first = i - width + 1
      do j = 1, size(methods)
        result%value(i, j) = estimate_method(data%open(first:i), data%high(first:i), &
          data%low(first:i), data%close(first:i), trim(methods(j)), keep_sign, remove_na)
      end do
    end do
  end function spread_fixed

  function spread_adaptive(data, width, methods, signed, na_rm) result(result)
    type(ohlc_data), intent(in) :: data
    integer, intent(in) :: width(:)
    character(len=*), intent(in) :: methods(:)
    logical, intent(in), optional :: signed, na_rm
    type(spread_series_result) :: result
    logical :: keep_sign, remove_na
    integer :: n, i, j, first, last

    keep_sign = .false.
    remove_na = .false.
    if (present(signed)) keep_sign = signed
    if (present(na_rm)) remove_na = na_rm
    n = data%size()
    allocate(result%method(size(methods)), result%value(n, size(methods)))
    result%method = methods
    result%value = nan_dp()
    if (.not. data%valid_dimensions()) then
      result%ok = .false.
      result%message = 'OHLC arrays must be allocated and have equal lengths'
      return
    end if
    if (size(width) == n) then
      do i = 1, n
        if (width(i) < 3 .or. width(i) > i) cycle
        first = i - width(i) + 1
        do j = 1, size(methods)
          result%value(i, j) = estimate_method(data%open(first:i), data%high(first:i), &
            data%low(first:i), data%close(first:i), trim(methods(j)), keep_sign, remove_na)
        end do
      end do
    else if (size(width) >= 2) then
      do i = 2, size(width)
        first = max(1, width(i - 1))
        last = width(i)
        if (last < first .or. last > n .or. last - first + 1 < 3) cycle
        do j = 1, size(methods)
          result%value(last, j) = estimate_method(data%open(first:last), data%high(first:last), &
            data%low(first:last), data%close(first:last), trim(methods(j)), keep_sign, remove_na)
        end do
      end do
    end if
  end function spread_adaptive

  function spread_expanding(data, methods, signed, na_rm) result(result)
    type(ohlc_data), intent(in) :: data
    character(len=*), intent(in) :: methods(:)
    logical, intent(in), optional :: signed, na_rm
    type(spread_series_result) :: result
    integer, allocatable :: widths(:)
    logical :: keep_sign, remove_na
    integer :: i, n

    keep_sign = .false.
    remove_na = .true.
    if (present(signed)) keep_sign = signed
    if (present(na_rm)) remove_na = na_rm
    n = data%size()
    allocate(widths(n))
    widths = [(i, i = 1, n)]
    result = spread_adaptive(data, widths, methods, keep_sign, remove_na)
  end function spread_expanding

  function spread_endpoints(data, endpoints, methods, signed, na_rm) result(result)
    type(ohlc_data), intent(in) :: data
    integer, intent(in) :: endpoints(:)
    character(len=*), intent(in) :: methods(:)
    logical, intent(in), optional :: signed, na_rm
    type(spread_series_result) :: result
    logical :: keep_sign, remove_na

    keep_sign = .false.
    remove_na = .false.
    if (present(signed)) keep_sign = signed
    if (present(na_rm)) remove_na = na_rm
    result = spread_adaptive(data, endpoints, methods, keep_sign, remove_na)
  end function spread_endpoints

end module bidask_windows
