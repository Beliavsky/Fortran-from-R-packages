! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2010-2023 Yang Lu and David Kane
! Copyright (C) 2026 Modern Fortran translation contributors
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License version 2 only.
module pa_brinson
  use pa_kinds, only: dp
  use pa_types, only: brinson_period_result, brinson_multi_result, attribution_summary
  use pa_utils, only: unique_int, geometric_product
  implicit none
  private
  public :: fit_brinson_period, fit_brinson_multi, summarize_brinson_multi

contains

  subroutine fit_brinson_period(category, benchmark_weight, portfolio_weight, asset_return, result)
    integer, intent(in) :: category(:)
    real(dp), intent(in) :: benchmark_weight(:), portfolio_weight(:), asset_return(:)
    type(brinson_period_result), intent(out) :: result
    real(dp), allocatable :: weighted_port(:), weighted_bench(:)
    integer :: j

    result%status = 0
    if (size(benchmark_weight) /= size(category) .or. size(portfolio_weight) /= size(category) .or. &
        size(asset_return) /= size(category)) then
      result%status = 1
      allocate(result%category(0), result%weight_portfolio(0), result%weight_benchmark(0), &
               result%return_portfolio(0), result%return_benchmark(0), result%category_effect(0,3))
      return
    end if

    call unique_int(category, result%category)
    allocate(result%weight_portfolio(size(result%category)), result%weight_benchmark(size(result%category)), &
             result%return_portfolio(size(result%category)), result%return_benchmark(size(result%category)), &
             result%category_effect(size(result%category), 3), weighted_port(size(result%category)), &
             weighted_bench(size(result%category)))

    do j = 1, size(result%category)
      result%weight_portfolio(j) = sum(portfolio_weight, mask=category == result%category(j))
      result%weight_benchmark(j) = sum(benchmark_weight, mask=category == result%category(j))
      weighted_port(j) = sum(portfolio_weight * asset_return, mask=category == result%category(j))
      weighted_bench(j) = sum(benchmark_weight * asset_return, mask=category == result%category(j))
      if (abs(result%weight_portfolio(j)) > tiny(1.0_dp)) then
        result%return_portfolio(j) = weighted_port(j) / result%weight_portfolio(j)
      else
        result%return_portfolio(j) = 0.0_dp
      end if
      if (abs(result%weight_benchmark(j)) > tiny(1.0_dp)) then
        result%return_benchmark(j) = weighted_bench(j) / result%weight_benchmark(j)
      else
        result%return_benchmark(j) = 0.0_dp
      end if
    end do

    result%q(1) = dot_product(result%return_portfolio, result%weight_portfolio)
    result%q(2) = dot_product(result%return_portfolio, result%weight_benchmark)
    result%q(3) = dot_product(result%return_benchmark, result%weight_portfolio)
    result%q(4) = dot_product(result%return_benchmark, result%weight_benchmark)

    result%category_effect(:, 1) = (result%weight_portfolio - result%weight_benchmark) * result%return_benchmark
    result%category_effect(:, 2) = (result%return_portfolio - result%return_benchmark) * result%weight_benchmark
    result%category_effect(:, 3) = (result%return_portfolio - result%return_benchmark) * &
                                   (result%weight_portfolio - result%weight_benchmark)

    result%aggregate(1) = result%q(3) - result%q(4)
    result%aggregate(2) = result%q(2) - result%q(4)
    result%aggregate(3) = result%q(1) - result%q(2) - result%q(3) + result%q(4)
    result%aggregate(4) = result%q(1) - result%q(4)
  end subroutine fit_brinson_period

  subroutine fit_brinson_multi(period, category, benchmark_weight, portfolio_weight, asset_return, result)
    integer, intent(in) :: period(:), category(:)
    real(dp), intent(in) :: benchmark_weight(:), portfolio_weight(:), asset_return(:)
    type(brinson_multi_result), intent(out) :: result
    type(brinson_period_result) :: one = brinson_period_result()
    integer, allocatable :: idx(:), subcat(:)
    real(dp), allocatable :: subb(:), subp(:), subr(:)
    integer :: i, j, k, g, nsub

    result%status = 0
    if (size(category) /= size(period) .or. size(benchmark_weight) /= size(period) .or. &
        size(portfolio_weight) /= size(period) .or. size(asset_return) /= size(period)) then
      result%status = 1
      allocate(result%period(0), result%category(0), result%weight_portfolio(0,0), &
               result%weight_benchmark(0,0), result%return_portfolio(0,0), &
               result%return_benchmark(0,0), result%category_effect(0,3,0), &
               result%q(4,0), result%raw(4,0))
      return
    end if

    call unique_int(period, result%period)
    call unique_int(category, result%category)
    allocate(result%weight_portfolio(size(result%category), size(result%period)), &
             result%weight_benchmark(size(result%category), size(result%period)), &
             result%return_portfolio(size(result%category), size(result%period)), &
             result%return_benchmark(size(result%category), size(result%period)), &
             result%category_effect(size(result%category), 3, size(result%period)), &
             result%q(4, size(result%period)), result%raw(4, size(result%period)))
    result%weight_portfolio = 0.0_dp
    result%weight_benchmark = 0.0_dp
    result%return_portfolio = 0.0_dp
    result%return_benchmark = 0.0_dp
    result%category_effect = 0.0_dp

    do i = 1, size(result%period)
      nsub = count(period == result%period(i))
      allocate(idx(nsub), subcat(nsub), subb(nsub), subp(nsub), subr(nsub))
      k = 0
      do j = 1, size(period)
        if (period(j) == result%period(i)) then
          k = k + 1
          idx(k) = j
          subcat(k) = category(j)
          subb(k) = benchmark_weight(j)
          subp(k) = portfolio_weight(j)
          subr(k) = asset_return(j)
        end if
      end do
      call fit_brinson_period(subcat, subb, subp, subr, one)
      if (one%status /= 0) result%status = one%status
      result%q(:, i) = one%q
      result%raw(:, i) = one%aggregate
      do j = 1, size(one%category)
        g = findloc(result%category, one%category(j), dim=1)
        result%weight_portfolio(g, i) = one%weight_portfolio(j)
        result%weight_benchmark(g, i) = one%weight_benchmark(j)
        result%return_portfolio(g, i) = one%return_portfolio(j)
        result%return_benchmark(g, i) = one%return_benchmark(j)
        result%category_effect(g, :, i) = one%category_effect(j, :)
      end do
      deallocate(idx, subcat, subb, subp, subr)
    end do
  end subroutine fit_brinson_multi

  subroutine summarize_brinson_multi(model, method, summary)
    type(brinson_multi_result), intent(in) :: model
    character(len=*), intent(in) :: method
    type(attribution_summary), intent(out) :: summary
    real(dp) :: gross_port, gross_bench, active_total, natural_scale, cvalue, denom
    real(dp) :: compounded_q(4)
    integer :: t

    summary%status = 0
    allocate(summary%raw(size(model%raw,1), size(model%raw,2)), summary%aggregate(4))
    summary%raw = model%raw

    select case (trim(method))
    case ('arithmetic')
      summary%aggregate = sum(model%raw, dim=2)
      allocate(summary%linking_coefficient(0))

    case ('geometric')
      do t = 1, 4
        compounded_q(t) = geometric_product(model%q(t, :)) - 1.0_dp
      end do
      summary%aggregate(1) = compounded_q(3) - compounded_q(4)
      summary%aggregate(2) = compounded_q(2) - compounded_q(4)
      summary%aggregate(3) = compounded_q(1) - compounded_q(2) - compounded_q(3) + compounded_q(4)
      summary%aggregate(4) = compounded_q(1) - compounded_q(4)
      allocate(summary%linking_coefficient(0))

    case ('linking')
      allocate(summary%linking_coefficient(size(model%raw, 2)))
      gross_port = geometric_product(model%q(1, :))
      gross_bench = geometric_product(model%q(4, :))
      active_total = gross_port - gross_bench
      t = size(model%raw, 2)
      denom = gross_port**(1.0_dp / real(t, dp)) - gross_bench**(1.0_dp / real(t, dp))
      if (abs(denom) <= sqrt(epsilon(1.0_dp))) then
        natural_scale = 1.0_dp
        summary%status = 2
      else
        natural_scale = active_total / real(t, dp) / denom
      end if
      denom = sum(model%raw(4, :)**2)
      if (abs(denom) <= sqrt(epsilon(1.0_dp))) then
        cvalue = 0.0_dp
        if (abs(active_total - natural_scale * sum(model%raw(4, :))) > 1.0e-10_dp) summary%status = 3
      else
        cvalue = (active_total - natural_scale * sum(model%raw(4, :))) / denom
      end if
      summary%linking_coefficient = natural_scale + cvalue * model%raw(4, :)
      summary%raw = model%raw * spread(summary%linking_coefficient, dim=1, ncopies=4)
      summary%aggregate = sum(summary%raw, dim=2)

    case default
      summary%status = 1
      summary%aggregate = 0.0_dp
      allocate(summary%linking_coefficient(0))
    end select
  end subroutine summarize_brinson_multi

end module pa_brinson
