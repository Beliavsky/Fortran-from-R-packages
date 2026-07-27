! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2010-2023 Yang Lu and David Kane
! Copyright (C) 2026 Modern Fortran translation contributors
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License version 2 only.
module pa_exposure
  use pa_kinds, only: dp
  use pa_types, only: exposure_result, exposure_multi_result
  use pa_utils, only: unique_int, quintile_groups
  implicit none
  private
  public :: categorical_exposure, continuous_exposure
  public :: categorical_exposure_multi, continuous_exposure_multi

contains

  subroutine categorical_exposure(group, portfolio_weight, benchmark_weight, result, status)
    integer, intent(in) :: group(:)
    real(dp), intent(in) :: portfolio_weight(:), benchmark_weight(:)
    type(exposure_result), intent(out) :: result
    integer, intent(out) :: status
    integer :: j

    status = 0
    if (size(portfolio_weight) /= size(group) .or. size(benchmark_weight) /= size(group)) then
      status = 1
      allocate(result%group(0), result%portfolio(0), result%benchmark(0), result%difference(0))
      return
    end if
    call unique_int(group, result%group)
    allocate(result%portfolio(size(result%group)), result%benchmark(size(result%group)), &
             result%difference(size(result%group)))
    do j = 1, size(result%group)
      result%portfolio(j) = sum(portfolio_weight, mask=group == result%group(j))
      result%benchmark(j) = sum(benchmark_weight, mask=group == result%group(j))
    end do
    result%difference = result%portfolio - result%benchmark
  end subroutine categorical_exposure

  subroutine continuous_exposure(value, portfolio_weight, benchmark_weight, result, status)
    real(dp), intent(in) :: value(:), portfolio_weight(:), benchmark_weight(:)
    type(exposure_result), intent(out) :: result
    integer, intent(out) :: status
    integer, allocatable :: group(:)
    integer :: j

    status = 0
    if (size(portfolio_weight) /= size(value) .or. size(benchmark_weight) /= size(value)) then
      status = 1
      allocate(result%group(0), result%portfolio(0), result%benchmark(0), result%difference(0))
      return
    end if
    call quintile_groups(value, group)
    allocate(result%group(5), result%portfolio(5), result%benchmark(5), result%difference(5))
    result%group = [(j, j=1,5)]
    do j = 1, 5
      result%portfolio(j) = sum(portfolio_weight, mask=group == j)
      result%benchmark(j) = sum(benchmark_weight, mask=group == j)
    end do
    result%difference = result%portfolio - result%benchmark
  end subroutine continuous_exposure

  subroutine categorical_exposure_multi(period, group, portfolio_weight, benchmark_weight, result, status)
    integer, intent(in) :: period(:), group(:)
    real(dp), intent(in) :: portfolio_weight(:), benchmark_weight(:)
    type(exposure_multi_result), intent(out) :: result
    integer, intent(out) :: status
    integer :: i, j

    status = 0
    if (size(group) /= size(period) .or. size(portfolio_weight) /= size(period) .or. &
        size(benchmark_weight) /= size(period)) then
      status = 1
      allocate(result%group(0), result%period(0), result%portfolio(0,0), &
               result%benchmark(0,0), result%difference(0,0))
      return
    end if
    call unique_int(period, result%period)
    call unique_int(group, result%group)
    allocate(result%portfolio(size(result%group), size(result%period)), &
             result%benchmark(size(result%group), size(result%period)), &
             result%difference(size(result%group), size(result%period)))
    do i = 1, size(result%period)
      do j = 1, size(result%group)
        result%portfolio(j, i) = sum(portfolio_weight, mask=(period == result%period(i) .and. group == result%group(j)))
        result%benchmark(j, i) = sum(benchmark_weight, mask=(period == result%period(i) .and. group == result%group(j)))
      end do
    end do
    result%difference = result%portfolio - result%benchmark
  end subroutine categorical_exposure_multi

  subroutine continuous_exposure_multi(period, value, portfolio_weight, benchmark_weight, result, status)
    integer, intent(in) :: period(:)
    real(dp), intent(in) :: value(:), portfolio_weight(:), benchmark_weight(:)
    type(exposure_multi_result), intent(out) :: result
    integer, intent(out) :: status
    integer, allocatable :: idx(:), group(:)
    real(dp), allocatable :: subvalue(:), subp(:), subb(:)
    integer :: i, j, k, nsub

    status = 0
    if (size(value) /= size(period) .or. size(portfolio_weight) /= size(period) .or. &
        size(benchmark_weight) /= size(period)) then
      status = 1
      allocate(result%group(0), result%period(0), result%portfolio(0,0), &
               result%benchmark(0,0), result%difference(0,0))
      return
    end if
    call unique_int(period, result%period)
    allocate(result%group(5), result%portfolio(5, size(result%period)), &
             result%benchmark(5, size(result%period)), result%difference(5, size(result%period)))
    result%group = [(j, j=1,5)]
    result%portfolio = 0.0_dp
    result%benchmark = 0.0_dp

    do i = 1, size(result%period)
      nsub = count(period == result%period(i))
      allocate(idx(nsub), subvalue(nsub), subp(nsub), subb(nsub))
      k = 0
      do j = 1, size(period)
        if (period(j) == result%period(i)) then
          k = k + 1
          idx(k) = j
          subvalue(k) = value(j)
          subp(k) = portfolio_weight(j)
          subb(k) = benchmark_weight(j)
        end if
      end do
      call quintile_groups(subvalue, group)
      do j = 1, 5
        result%portfolio(j, i) = sum(subp, mask=group == j)
        result%benchmark(j, i) = sum(subb, mask=group == j)
      end do
      deallocate(idx, subvalue, subp, subb, group)
    end do
    result%difference = result%portfolio - result%benchmark
  end subroutine continuous_exposure_multi

end module pa_exposure
