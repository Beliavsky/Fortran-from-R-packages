! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2010-2023 Yang Lu and David Kane
! Copyright (C) 2026 Modern Fortran translation contributors
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License version 2 only.
module pa_regression
  use pa_kinds, only: dp
  use pa_types, only: regression_period_result, regression_multi_result, regression_summary
  use pa_linalg, only: least_squares
  use pa_utils, only: unique_int, geometric_product
  implicit none
  private
  public :: fit_regression_period, fit_regression_multi
  public :: summarize_regression_period, summarize_regression_multi

contains

  subroutine fit_regression_period(asset_return, design, benchmark_weight, portfolio_weight, result)
    real(dp), intent(in) :: asset_return(:), design(:, :), benchmark_weight(:), portfolio_weight(:)
    type(regression_period_result), intent(out) :: result
    real(dp), allocatable :: active_weight(:)

    result%status = 0
    if (size(design,1) /= size(asset_return) .or. size(benchmark_weight) /= size(asset_return) .or. &
        size(portfolio_weight) /= size(asset_return)) then
      result%status = 1
      allocate(result%coefficients(0), result%active_exposure(0), result%contribution(0))
      return
    end if

    call least_squares(design, asset_return, result%coefficients, result%rank, result%status)
    if (result%status /= 0) then
      allocate(result%active_exposure(size(design,2)), result%contribution(size(design,2)))
      result%active_exposure = 0.0_dp
      result%contribution = 0.0_dp
      return
    end if

    allocate(active_weight(size(asset_return)), result%active_exposure(size(design,2)), &
             result%contribution(size(design,2)))
    active_weight = portfolio_weight - benchmark_weight
    result%active_exposure = matmul(transpose(design), active_weight)
    result%contribution = result%active_exposure * result%coefficients
    result%benchmark_return = dot_product(benchmark_weight, asset_return)
    result%portfolio_return = dot_product(portfolio_weight, asset_return)
    result%active_return = result%portfolio_return - result%benchmark_return
  end subroutine fit_regression_period

  subroutine fit_regression_multi(period, asset_return, design, benchmark_weight, portfolio_weight, result)
    integer, intent(in) :: period(:)
    real(dp), intent(in) :: asset_return(:), design(:, :), benchmark_weight(:), portfolio_weight(:)
    type(regression_multi_result), intent(out) :: result
    type(regression_period_result) :: one = regression_period_result()
    integer, allocatable :: idx(:)
    real(dp), allocatable :: subreturn(:), subdesign(:, :), subb(:), subp(:)
    integer :: i, j, k, nsub

    result%status = 0
    if (size(asset_return) /= size(period) .or. size(design,1) /= size(period) .or. &
        size(benchmark_weight) /= size(period) .or. size(portfolio_weight) /= size(period)) then
      result%status = 1
      allocate(result%period(0), result%coefficients(0,0), result%active_exposure(0,0), &
               result%contribution(0,0), result%portfolio_return(0), result%benchmark_return(0), &
               result%active_return(0), result%rank(0))
      return
    end if

    call unique_int(period, result%period)
    allocate(result%coefficients(size(design,2), size(result%period)), &
             result%active_exposure(size(design,2), size(result%period)), &
             result%contribution(size(design,2), size(result%period)), &
             result%portfolio_return(size(result%period)), result%benchmark_return(size(result%period)), &
             result%active_return(size(result%period)), result%rank(size(result%period)))

    do i = 1, size(result%period)
      nsub = count(period == result%period(i))
      allocate(idx(nsub), subreturn(nsub), subdesign(nsub,size(design,2)), subb(nsub), subp(nsub))
      k = 0
      do j = 1, size(period)
        if (period(j) == result%period(i)) then
          k = k + 1
          idx(k) = j
          subreturn(k) = asset_return(j)
          subdesign(k,:) = design(j,:)
          subb(k) = benchmark_weight(j)
          subp(k) = portfolio_weight(j)
        end if
      end do
      call fit_regression_period(subreturn, subdesign, subb, subp, one)
      if (one%status /= 0) result%status = one%status
      result%coefficients(:,i) = one%coefficients
      result%active_exposure(:,i) = one%active_exposure
      result%contribution(:,i) = one%contribution
      result%portfolio_return(i) = one%portfolio_return
      result%benchmark_return(i) = one%benchmark_return
      result%active_return(i) = one%active_return
      result%rank(i) = one%rank
      deallocate(idx, subreturn, subdesign, subb, subp)
    end do
  end subroutine fit_regression_multi

  subroutine summarize_regression_period(model, group_start, group_end, values, status)
    type(regression_period_result), intent(in) :: model
    integer, intent(in) :: group_start(:), group_end(:)
    real(dp), allocatable, intent(out) :: values(:)
    integer, intent(out) :: status
    integer :: i, nvar

    status = 0
    nvar = size(group_start)
    allocate(values(nvar + 4))
    values = 0.0_dp
    if (size(group_end) /= nvar) then
      status = 1
      return
    end if
    do i = 1, nvar
      if (group_start(i) < 1 .or. group_end(i) > size(model%contribution) .or. group_start(i) > group_end(i)) then
        status = 2
        return
      end if
      values(i) = sum(model%contribution(group_start(i):group_end(i)))
    end do
    values(nvar + 1) = model%active_return - sum(model%contribution)
    values(nvar + 2) = model%portfolio_return
    values(nvar + 3) = model%benchmark_return
    values(nvar + 4) = model%active_return
  end subroutine summarize_regression_period

  subroutine summarize_regression_multi(model, group_start, group_end, method, summary)
    type(regression_multi_result), intent(in) :: model
    integer, intent(in) :: group_start(:), group_end(:)
    character(len=*), intent(in) :: method
    type(regression_summary), intent(out) :: summary
    real(dp), allocatable :: full_raw(:, :), one(:)
    real(dp) :: gross_port, gross_bench, active_total, ascale, cvalue, denom
    integer :: i, nvar, nrow, status_one

    summary%status = 0
    nvar = size(group_start)
    nrow = nvar + 4
    allocate(full_raw(nrow, size(model%period)))
    do i = 1, size(model%period)
      call summarize_regression_period(regression_period_result( &
        coefficients=model%coefficients(:,i), active_exposure=model%active_exposure(:,i), &
        contribution=model%contribution(:,i), portfolio_return=model%portfolio_return(i), &
        benchmark_return=model%benchmark_return(i), active_return=model%active_return(i), &
        rank=model%rank(i), status=0), group_start, group_end, one, status_one)
      if (status_one /= 0) summary%status = status_one
      full_raw(:,i) = one
      deallocate(one)
    end do

    select case (trim(method))
    case ('arithmetic')
      allocate(summary%raw(nrow,size(model%period)), summary%aggregate(nrow), summary%linking_coefficient(0))
      summary%raw = full_raw
      summary%aggregate = sum(full_raw, dim=2)

    case ('geometric')
      allocate(summary%raw(nrow,size(model%period)), summary%aggregate(nrow), summary%linking_coefficient(0))
      summary%raw = full_raw
      do i = 1, nrow
        summary%aggregate(i) = geometric_product(full_raw(i,:)) - 1.0_dp
      end do
      summary%aggregate(nvar + 4) = summary%aggregate(nvar + 2) - summary%aggregate(nvar + 3)

    case ('linking')
      allocate(summary%raw(nvar+2,size(model%period)), summary%aggregate(nvar+2), &
               summary%linking_coefficient(size(model%period)))
      gross_port = geometric_product(model%portfolio_return)
      gross_bench = geometric_product(model%benchmark_return)
      active_total = gross_port - gross_bench
      denom = gross_port**(1.0_dp/real(size(model%period),dp)) - &
              gross_bench**(1.0_dp/real(size(model%period),dp))
      if (abs(denom) <= sqrt(epsilon(1.0_dp))) then
        ascale = 1.0_dp
        summary%status = 2
      else
        ascale = active_total / real(size(model%period),dp) / denom
      end if
      denom = sum(model%active_return**2)
      if (abs(denom) <= sqrt(epsilon(1.0_dp))) then
        cvalue = 0.0_dp
        if (abs(active_total - ascale*sum(model%active_return)) > 1.0e-10_dp) summary%status = 3
      else
        cvalue = (active_total - ascale*sum(model%active_return)) / denom
      end if
      summary%linking_coefficient = ascale + cvalue*model%active_return
      summary%raw(1:nvar+1,:) = full_raw(1:nvar+1,:) * &
        spread(summary%linking_coefficient, dim=1, ncopies=nvar+1)
      summary%raw(nvar+2,:) = full_raw(nvar+4,:) * summary%linking_coefficient
      summary%aggregate = sum(summary%raw, dim=2)

    case default
      allocate(summary%raw(0,0), summary%aggregate(0), summary%linking_coefficient(0))
      summary%status = 1
    end select
  end subroutine summarize_regression_multi

end module pa_regression
