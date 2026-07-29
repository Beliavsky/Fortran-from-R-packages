! stockAnalyst-fortran
! Copyright (C) 2022 MaheshP Kumar (original R package)
! Copyright (C) 2026 Fortran port contributors
! SPDX-License-Identifier: GPL-3.0-only
!
! This file is part of stockAnalyst-fortran.
! It is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, version 3 of the License.

module stock_analyst
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, ieee_value
  implicit none
  private

  integer, parameter, public :: dp = kind(1.0d0)

  public :: round_to
  public :: share_value_using_ddm_1yr, share_value_using_ddm_n_years
  public :: share_value_ggm_constant_growth, share_value_preferred_stock
  public :: share_value_ggm_negative_growth, computing_g_using_ggm
  public :: justified_leading_pe, justified_trailing_pe, computing_r_with_ggm
  public :: share_val_using_two_stage_ddm, share_val_using_three_stage_ddm
  public :: share_val_using_two_stage_hmodel, share_value_no_current_dividend
  public :: annualized_hpr, firm_value_using_disc_fcff, equity_value_given_debt_mv
  public :: share_value_given_debt_mv, share_value_using_disc_fcfe
  public :: firm_value_constant_g, equity_value_constant_g, share_val_constant_g
  public :: share_val_two_stage, share_val_three_stage
  public :: share_value_ri, share_value_computed_ri, computing_abs_ri, computing_ri
  public :: share_value_roe, single_stage_r, share_value_ri_multi_stage_eps
  public :: share_value_ri_multi_stage_roe, share_value_ri_plus_pvtv
  public :: trailing_pe_basic_eps, trailing_pe_diluted_eps, earning_yield_ep
  public :: leading_pe_next_4qs, leading_fy1_pe, leading_fy2_pe
  public :: predicted_pe_on_csr, forward_peg, predicted_pe_by_fed_model
  public :: implied_pe_by_yardeni_model, share_price_using_past_pe
  public :: pe_for_pass_through_inflation, terminal_value_using_pe
  public :: computing_bv_per_share, computing_pb, computing_ps
  public :: computing_ev_dollar_val, computing_ev_multiple
  public :: computing_sustainable_g, computing_r_with_capm, computing_wacc
  public :: computing_r_with_ffm, computing_r_with_hmodel

  ! Compatibility generics corresponding to the original exported R names.
  public :: sharevalueusingddm1yr, sharevalueusingddmnyrs
  public :: sharevalueggmconstantgrowth, sharevaluepreferredstock
  public :: sharevalueggmnegativegrowth, computinggusingggm
  public :: justifiedleadingpe, justifiedtrailingpe, computingrwithggm
  public :: sharevalusingtwostageddm, sharevalusingthreestageddm
  public :: sharevalusingtwostagehmodel, sharevaluenocurrentdivdend
  public :: annulizedhpr, firmvalueusingdiscfcff, equityvaluegivendebtmv
  public :: sharevaluegivendebtmv, sharevalueusingdiscfcfe
  public :: firmvalueconstantg, equityvalueconstantg, sharevalconstantg
  public :: sharevaltwostage, sharevalthreestg
  public :: sharevalueri, sharevaluecomputedri, computingabsri, computingri
  public :: sharevalueroe, singlestager, sharevaluerimultistageeps
  public :: sharevaluerimultistg, sharevalueripluspvtv
  public :: trailingpebasiceps, trailingpedilutedeps, earningyieldep
  public :: leadingpenext4qs, leadingfy1pe, leadingfy2pe
  public :: predictedpeoncsr, forwardpeg, predictedpebyfedmodel
  public :: impliedpebyyardenimodel, sharepriceusingpastpe
  public :: peforpassthroughinflation, terminalvalueusingpe
  public :: computingbvpershare, computingpb, computingps
  public :: computingevdollarval, computingevmultiple
  public :: computingsustainableg, computingrwithcapm, computingwacc
  public :: computingrwithffm, computingrwithhmodel

  interface sharevalueusingddm1yr
    module procedure share_value_using_ddm_1yr
  end interface
  interface sharevalueusingddmnyrs
    module procedure share_value_using_ddm_n_years
  end interface
  interface sharevalueggmconstantgrowth
    module procedure share_value_ggm_constant_growth
  end interface
  interface sharevaluepreferredstock
    module procedure share_value_preferred_stock
  end interface
  interface sharevalueggmnegativegrowth
    module procedure share_value_ggm_negative_growth
  end interface
  interface computinggusingggm
    module procedure computing_g_using_ggm
  end interface
  interface justifiedleadingpe
    module procedure justified_leading_pe
  end interface
  interface justifiedtrailingpe
    module procedure justified_trailing_pe
  end interface
  interface computingrwithggm
    module procedure computing_r_with_ggm
  end interface
  interface sharevalusingtwostageddm
    module procedure share_val_using_two_stage_ddm
  end interface
  interface sharevalusingthreestageddm
    module procedure share_val_using_three_stage_ddm
  end interface
  interface sharevalusingtwostagehmodel
    module procedure share_val_using_two_stage_hmodel
  end interface
  interface sharevaluenocurrentdivdend
    module procedure share_value_no_current_dividend
  end interface
  interface annulizedhpr
    module procedure annualized_hpr
  end interface
  interface firmvalueusingdiscfcff
    module procedure firm_value_using_disc_fcff
  end interface
  interface equityvaluegivendebtmv
    module procedure equity_value_given_debt_mv
  end interface
  interface sharevaluegivendebtmv
    module procedure share_value_given_debt_mv
  end interface
  interface sharevalueusingdiscfcfe
    module procedure share_value_using_disc_fcfe
  end interface
  interface firmvalueconstantg
    module procedure firm_value_constant_g
  end interface
  interface equityvalueconstantg
    module procedure equity_value_constant_g
  end interface
  interface sharevalconstantg
    module procedure share_val_constant_g
  end interface
  interface sharevaltwostage
    module procedure share_val_two_stage
  end interface
  interface sharevalthreestg
    module procedure share_val_three_stage
  end interface
  interface sharevalueri
    module procedure share_value_ri
  end interface
  interface sharevaluecomputedri
    module procedure share_value_computed_ri
  end interface
  interface computingabsri
    module procedure computing_abs_ri
  end interface
  interface computingri
    module procedure computing_ri
  end interface
  interface sharevalueroe
    module procedure share_value_roe
  end interface
  interface singlestager
    module procedure single_stage_r
  end interface
  interface sharevaluerimultistageeps
    module procedure share_value_ri_multi_stage_eps
  end interface
  interface sharevaluerimultistg
    module procedure share_value_ri_multi_stage_roe
  end interface
  interface sharevalueripluspvtv
    module procedure share_value_ri_plus_pvtv
  end interface
  interface trailingpebasiceps
    module procedure trailing_pe_basic_eps
  end interface
  interface trailingpedilutedeps
    module procedure trailing_pe_diluted_eps
  end interface
  interface earningyieldep
    module procedure earning_yield_ep
  end interface
  interface leadingpenext4qs
    module procedure leading_pe_next_4qs
  end interface
  interface leadingfy1pe
    module procedure leading_fy1_pe
  end interface
  interface leadingfy2pe
    module procedure leading_fy2_pe
  end interface
  interface predictedpeoncsr
    module procedure predicted_pe_on_csr
  end interface
  interface forwardpeg
    module procedure forward_peg
  end interface
  interface predictedpebyfedmodel
    module procedure predicted_pe_by_fed_model
  end interface
  interface impliedpebyyardenimodel
    module procedure implied_pe_by_yardeni_model
  end interface
  interface sharepriceusingpastpe
    module procedure share_price_using_past_pe
  end interface
  interface peforpassthroughinflation
    module procedure pe_for_pass_through_inflation
  end interface
  interface terminalvalueusingpe
    module procedure terminal_value_using_pe
  end interface
  interface computingbvpershare
    module procedure computing_bv_per_share
  end interface
  interface computingpb
    module procedure computing_pb
  end interface
  interface computingps
    module procedure computing_ps
  end interface
  interface computingevdollarval
    module procedure computing_ev_dollar_val
  end interface
  interface computingevmultiple
    module procedure computing_ev_multiple
  end interface
  interface computingsustainableg
    module procedure computing_sustainable_g
  end interface
  interface computingrwithcapm
    module procedure computing_r_with_capm
  end interface
  interface computingwacc
    module procedure computing_wacc
  end interface
  interface computingrwithffm
    module procedure computing_r_with_ffm
  end interface
  interface computingrwithhmodel
    module procedure computing_r_with_hmodel
  end interface

contains

  pure elemental real(dp) function round_to(x, digits) result(value)
    real(dp), intent(in) :: x
    integer, intent(in) :: digits
    real(dp) :: scale, ax, lower, fraction, rounded, tolerance

    if (.not. ieee_is_finite(x)) then
      value = x
      return
    end if

    scale = 10.0_dp**digits
    ax = abs(x*scale)
    lower = floor(ax)
    fraction = ax - lower
    tolerance = 16.0_dp*epsilon(max(1.0_dp, ax))

    if (fraction > 0.5_dp + tolerance) then
      rounded = lower + 1.0_dp
    else if (fraction < 0.5_dp - tolerance) then
      rounded = lower
    else if (int(modulo(lower, 2.0_dp)) == 0) then
      rounded = lower
    else
      rounded = lower + 1.0_dp
    end if

    value = sign(rounded/scale, x)
  end function round_to

  pure real(dp) function nan_value() result(value)
    value = ieee_value(0.0_dp, ieee_quiet_nan)
  end function nan_value

  pure function lowercase(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: i, code

    lower = text
    do i = 1, len(text)
      code = iachar(lower(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lower(i:i) = achar(code + iachar('a') - iachar('A'))
      end if
    end do
  end function lowercase

  pure real(dp) function median_value(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: work(size(x)), key
    integer :: i, j, n

    n = size(x)
    if (n == 0) then
      value = nan_value()
      return
    end if

    work = x
    do i = 2, n
      key = work(i)
      j = i - 1
      do while (j >= 1)
        if (work(j) <= key) exit
        work(j + 1) = work(j)
        j = j - 1
      end do
      work(j + 1) = key
    end do

    if (mod(n, 2) == 1) then
      value = work((n + 1)/2)
    else
      value = 0.5_dp*(work(n/2) + work(n/2 + 1))
    end if
  end function median_value

  pure real(dp) function discounted_sum(cash_flow, times, rate) result(value)
    real(dp), intent(in) :: cash_flow(:), times(:), rate
    if (size(cash_flow) /= size(times)) then
      value = nan_value()
    else
      value = sum(cash_flow/(1.0_dp + rate)**times)
    end if
  end function discounted_sum

  pure real(dp) function share_value_using_ddm_1yr(dividend1yr, exp_share_price_in_1yr, n, r) result(value)
    real(dp), intent(in) :: dividend1yr, exp_share_price_in_1yr, n, r
    value = round_to((dividend1yr + exp_share_price_in_1yr)/(1.0_dp + r)**n, 2)
  end function share_value_using_ddm_1yr

  pure real(dp) function share_value_using_ddm_n_years(dividend, exp_share_price_nyr, times, n, r) result(value)
    real(dp), intent(in) :: dividend(:), exp_share_price_nyr, times(:), n, r
    value = round_to(discounted_sum(dividend, times, r) + exp_share_price_nyr/(1.0_dp + r)**n, 2)
  end function share_value_using_ddm_n_years

  pure real(dp) function share_value_ggm_constant_growth(dividend, r, g, div_n) result(value)
    real(dp), intent(in) :: dividend, r, g
    integer, intent(in) :: div_n
    if (div_n /= 1) then
      value = round_to(dividend*(1.0_dp + g)/(r - g), 2)
    else
      value = round_to(dividend/(r - g), 2)
    end if
  end function share_value_ggm_constant_growth

  pure real(dp) function share_value_preferred_stock(dividend, r) result(value)
    real(dp), intent(in) :: dividend, r
    value = round_to(dividend/r, 2)
  end function share_value_preferred_stock

  pure real(dp) function share_value_ggm_negative_growth(dividend, r, neg_g) result(value)
    real(dp), intent(in) :: dividend, r, neg_g
    if (neg_g <= 0.0_dp) then
      value = round_to(dividend/(r - neg_g), 2)
    else
      value = round_to(dividend/(r + neg_g), 2)
    end if
  end function share_value_ggm_negative_growth

  pure real(dp) function computing_g_using_ggm(div_not, r, share_price) result(value)
    real(dp), intent(in) :: div_not, r, share_price
    value = round_to((r*share_price - div_not)/(div_not + share_price), 4)
  end function computing_g_using_ggm

  pure real(dp) function justified_leading_pe(r_capm, payout_ratio, g) result(value)
    real(dp), intent(in) :: r_capm, payout_ratio, g
    value = round_to(payout_ratio/(r_capm - g), 1)
  end function justified_leading_pe

  pure real(dp) function justified_trailing_pe(r_capm, payout_ratio, g) result(value)
    real(dp), intent(in) :: r_capm, payout_ratio, g
    value = round_to(payout_ratio*(1.0_dp + g)/(r_capm - g), 1)
  end function justified_trailing_pe

  pure real(dp) function computing_r_with_ggm(div_n1, g, sp_not) result(value)
    real(dp), intent(in) :: div_n1, g, sp_not
    value = round_to(div_n1/sp_not + g, 4)
  end function computing_r_with_ggm

  pure real(dp) function share_val_using_two_stage_ddm(div_not, r, n, g_s, g_l) result(value)
    real(dp), intent(in) :: div_not, r, n, g_s, g_l
    value = round_to(div_not*(1.0_dp + g_s)**n*(1.0_dp + g_l)/(r - g_l), 4)
  end function share_val_using_two_stage_ddm

  pure real(dp) function share_val_using_three_stage_ddm(div_not, r, n1, n2, g1, g2, g3) result(value)
    real(dp), intent(in) :: div_not, r, n1, n2, g1, g2, g3
    value = round_to(div_not*(1.0_dp + g1)**n1*(1.0_dp + g2)**n2*(1.0_dp + g3)/(r - g3), 4)
  end function share_val_using_three_stage_ddm

  pure real(dp) function share_val_using_two_stage_hmodel(div_not, r, n, h, g_s, g_l) result(value)
    real(dp), intent(in) :: div_not, r, n, h, g_s, g_l
    if (.not. ieee_is_finite(n)) then
      value = nan_value()
    else
      value = round_to(div_not*(1.0_dp + g_l)/(r - g_l) + div_not*h*(g_s - g_l)/(r - g_l), 2)
    end if
  end function share_val_using_two_stage_hmodel

  pure real(dp) function share_value_no_current_dividend(div_n, t, g, r) result(value)
    real(dp), intent(in) :: div_n, t, g, r
    value = round_to((div_n/(r - g))/(1.0_dp + r)**(t - 1.0_dp), 2)
  end function share_value_no_current_dividend

  pure real(dp) function annualized_hpr(total_per_share_dividend_hp, sp_h, sp_not, n) result(value)
    real(dp), intent(in) :: total_per_share_dividend_hp, sp_h, sp_not, n
    real(dp) :: hpr
    hpr = total_per_share_dividend_hp/sp_not + (sp_h - sp_not)/sp_not
    value = round_to((1.0_dp + hpr)**(1.0_dp/n) - 1.0_dp, 4)
  end function annualized_hpr

  pure real(dp) function firm_value_using_disc_fcff(fcff, times, wacc) result(value)
    real(dp), intent(in) :: fcff(:), times(:), wacc
    value = round_to(discounted_sum(fcff, times, wacc), 6)
  end function firm_value_using_disc_fcff

  pure real(dp) function equity_value_given_debt_mv(fcff, t, wacc, debt_mv) result(value)
    real(dp), intent(in) :: fcff(:), t(:), wacc, debt_mv
    value = round_to(discounted_sum(fcff, t, wacc) - debt_mv, 6)
  end function equity_value_given_debt_mv

  pure real(dp) function share_value_given_debt_mv(fcff, t, wacc, debt_mv, shares) result(value)
    real(dp), intent(in) :: fcff(:), t(:), wacc, debt_mv, shares
    value = round_to((discounted_sum(fcff, t, wacc) - debt_mv)/shares, 2)
  end function share_value_given_debt_mv

  pure real(dp) function share_value_using_disc_fcfe(fcfe, t, r, shares) result(value)
    real(dp), intent(in) :: fcfe(:), t(:), r, shares
    value = round_to(discounted_sum(fcfe, t, r)/shares, 2)
  end function share_value_using_disc_fcfe

  pure real(dp) function firm_value_constant_g(fcff0, g, wacc) result(value)
    real(dp), intent(in) :: fcff0, g, wacc
    value = round_to(fcff0*(1.0_dp + g)/(wacc - g), 2)
  end function firm_value_constant_g

  pure real(dp) function equity_value_constant_g(fcff0, g, wacc, debt_val) result(value)
    real(dp), intent(in) :: fcff0, g, wacc, debt_val
    value = round_to(fcff0*(1.0_dp + g)/(wacc - g) - debt_val, 2)
  end function equity_value_constant_g

  pure real(dp) function share_val_constant_g(fcfe0, g, wacc, shares) result(value)
    real(dp), intent(in) :: fcfe0, g, wacc, shares
    value = round_to((fcfe0*(1.0_dp + g)/(wacc - g))/shares, 2)
  end function share_val_constant_g

  pure real(dp) function share_val_two_stage(fcfe, t, growth, r, shares) result(value)
    real(dp), intent(in) :: fcfe(:), t(:), growth(:), r, shares
    if (size(fcfe) /= size(growth)) then
      value = nan_value()
    else
      value = round_to(discounted_sum(fcfe*growth, t, r)/shares, 2)
    end if
  end function share_val_two_stage

  pure real(dp) function share_val_three_stage(fcfe, t, growth, r, shares) result(value)
    real(dp), intent(in) :: fcfe(:), t(:), growth(:), r, shares
    if (size(fcfe) /= size(growth)) then
      value = nan_value()
    else
      value = round_to(discounted_sum(fcfe*growth, t, r)/shares, 2)
    end if
  end function share_val_three_stage

  pure real(dp) function share_value_ri(bgn_bvps, ri, r, times) result(value)
    real(dp), intent(in) :: bgn_bvps, ri(:), r, times(:)
    value = round_to(bgn_bvps + discounted_sum(ri, times, r), 2)
  end function share_value_ri

  pure real(dp) function share_value_computed_ri(bgn_bvps, eps, r, times) result(value)
    real(dp), intent(in) :: bgn_bvps(:), eps(:), r, times(:)
    if (size(bgn_bvps) == 0 .or. size(bgn_bvps) /= size(eps)) then
      value = nan_value()
    else
      value = round_to(bgn_bvps(1) + discounted_sum(eps - r*bgn_bvps, times, r), 2)
    end if
  end function share_value_computed_ri

  pure function computing_abs_ri(ebit, debt, equity, r, rd, tax_rate) result(ri)
    real(dp), intent(in) :: ebit(:), debt(:), equity(:), r, rd, tax_rate
    real(dp) :: ri(size(ebit)), pre_tax_income(size(ebit)), net_income(size(ebit))
    if (size(debt) /= size(ebit) .or. size(equity) /= size(ebit)) then
      ri = nan_value()
      return
    end if
    pre_tax_income = ebit - rd*debt
    net_income = pre_tax_income - tax_rate*pre_tax_income
    ri = round_to(net_income - r*equity, 4)
  end function computing_abs_ri

  pure function computing_ri(bgn_bvps, eps, r) result(ri)
    real(dp), intent(in) :: bgn_bvps(:), eps(:), r
    real(dp) :: ri(size(bgn_bvps))
    if (size(eps) /= size(bgn_bvps)) then
      ri = nan_value()
    else
      ri = round_to(eps - r*bgn_bvps, 4)
    end if
  end function computing_ri

  pure real(dp) function share_value_roe(roe, bgn_bvps, r, times) result(value)
    real(dp), intent(in) :: roe(:), bgn_bvps(:), r, times(:)
    if (size(bgn_bvps) == 0 .or. size(roe) /= size(bgn_bvps)) then
      value = nan_value()
    else
      value = round_to(bgn_bvps(1) + discounted_sum((roe - r)*bgn_bvps, times, r), 2)
    end if
  end function share_value_roe

  pure real(dp) function single_stage_r(roe_rate, bgn_bvps, r, g) result(value)
    real(dp), intent(in) :: roe_rate, bgn_bvps, r, g
    value = round_to(bgn_bvps + ((roe_rate - r)/(r - g))*bgn_bvps, 2)
  end function single_stage_r

  pure real(dp) function share_value_ri_multi_stage_eps(bgn_bvps, eps, r, times, premium, n) result(value)
    real(dp), intent(in) :: bgn_bvps(:), eps(:), r, times(:), premium, n
    if (size(bgn_bvps) == 0 .or. size(bgn_bvps) /= size(eps)) then
      value = nan_value()
    else
      value = round_to(bgn_bvps(1) + discounted_sum(eps - r*bgn_bvps, times, r) + premium/(1.0_dp + r)**n, 2)
    end if
  end function share_value_ri_multi_stage_eps

  pure real(dp) function share_value_ri_multi_stage_roe(roe, bgn_bv, r, times, premium, n) result(value)
    real(dp), intent(in) :: roe(:), bgn_bv(:), r, times(:), premium, n
    if (size(bgn_bv) == 0 .or. size(roe) /= size(bgn_bv)) then
      value = nan_value()
    else
      value = round_to(bgn_bv(1) + discounted_sum((roe - r)*bgn_bv, times, r) + premium/(1.0_dp + r)**n, 2)
    end if
  end function share_value_ri_multi_stage_roe

  pure real(dp) function share_value_ri_plus_pvtv(bgn_bvps, eps, r, times, persistence_factor, n) result(value)
    real(dp), intent(in) :: bgn_bvps(:), eps(:), r, times(:), persistence_factor
    integer, intent(in) :: n
    if (size(bgn_bvps) == 0 .or. size(bgn_bvps) /= size(eps) .or. n < 1 .or. n > size(eps)) then
      value = nan_value()
    else
      value = bgn_bvps(1) + discounted_sum(eps - r*bgn_bvps, times, r)
      value = value + (eps(n) - r*bgn_bvps(n))/((1.0_dp + r - persistence_factor)*(1.0_dp + r)**n)
      value = round_to(value, 3)
    end if
  end function share_value_ri_plus_pvtv

  pure real(dp) function trailing_pe_basic_eps(current_share_price, basic_eps) result(value)
    real(dp), intent(in) :: current_share_price, basic_eps
    value = round_to(current_share_price/basic_eps, 1)
  end function trailing_pe_basic_eps

  pure real(dp) function trailing_pe_diluted_eps(current_share_price, diluted_eps) result(value)
    real(dp), intent(in) :: current_share_price, diluted_eps
    value = round_to(current_share_price/diluted_eps, 1)
  end function trailing_pe_diluted_eps

  pure real(dp) function earning_yield_ep(current_share_price, ttm_diluted_eps) result(value)
    real(dp), intent(in) :: current_share_price, ttm_diluted_eps
    value = round_to(ttm_diluted_eps/current_share_price, 4)
  end function earning_yield_ep

  pure real(dp) function leading_pe_next_4qs(current_share_price, q1_eps, q2_eps, q3_eps, q4_eps) result(value)
    real(dp), intent(in) :: current_share_price, q1_eps, q2_eps, q3_eps, q4_eps
    value = round_to(current_share_price/(q1_eps + q2_eps + q3_eps + q4_eps), 1)
  end function leading_pe_next_4qs

  pure real(dp) function leading_fy1_pe(current_share_price, fy1_eps) result(value)
    real(dp), intent(in) :: current_share_price, fy1_eps
    value = round_to(current_share_price/fy1_eps, 1)
  end function leading_fy1_pe

  pure real(dp) function leading_fy2_pe(current_share_price, fy2_eps) result(value)
    real(dp), intent(in) :: current_share_price, fy2_eps
    value = round_to(current_share_price/fy2_eps, 1)
  end function leading_fy2_pe

  pure real(dp) function predicted_pe_on_csr(b0, b1, b2, b3, x1_drp, x2_beta, x3_egr) result(value)
    real(dp), intent(in) :: b0, b1, b2, b3, x1_drp, x2_beta, x3_egr
    value = round_to(b0 + b1*x1_drp + b2*x2_beta + b3*x3_egr, 1)
  end function predicted_pe_on_csr

  pure real(dp) function forward_peg(leading_pe, percent_eps_growth) result(value)
    real(dp), intent(in) :: leading_pe, percent_eps_growth
    value = round_to(leading_pe/percent_eps_growth, 2)
  end function forward_peg

  pure real(dp) function predicted_pe_by_fed_model(ten_year_bond_yield) result(value)
    real(dp), intent(in) :: ten_year_bond_yield
    value = round_to(1.0_dp/ten_year_bond_yield, 1)
  end function predicted_pe_by_fed_model

  pure real(dp) function implied_pe_by_yardeni_model(cby, b, lteg, residual_value) result(value)
    real(dp), intent(in) :: cby, b, lteg, residual_value
    value = round_to(1.0_dp/(cby - b*lteg - residual_value), 1)
  end function implied_pe_by_yardeni_model

  pure real(dp) function share_price_using_past_pe(average_mode, historical_pes, recent_eps) result(value)
    character(len=*), intent(in) :: average_mode
    real(dp), intent(in) :: historical_pes(:), recent_eps
    real(dp) :: benchmark_pe
    if (size(historical_pes) == 0) then
      value = nan_value()
      return
    end if
    if (trim(lowercase(average_mode)) == 'median') then
      benchmark_pe = median_value(historical_pes)
    else
      benchmark_pe = sum(historical_pes)/real(size(historical_pes), dp)
    end if
    value = round_to(benchmark_pe*recent_eps, 0)
  end function share_price_using_past_pe

  pure real(dp) function pe_for_pass_through_inflation(real_ror, inflation, pass_through_rate) result(value)
    real(dp), intent(in) :: real_ror, inflation, pass_through_rate
    value = round_to(1.0_dp/(real_ror + inflation*(1.0_dp - pass_through_rate)), 1)
  end function pe_for_pass_through_inflation

  pure real(dp) function terminal_value_using_pe(mode, benchmark_pe, earnings_n, payout, g, r) result(value)
    character(len=*), intent(in) :: mode
    real(dp), intent(in) :: benchmark_pe, earnings_n, payout, g, r
    if (trim(lowercase(mode)) == 'ggm') then
      value = round_to(earnings_n*payout*(1.0_dp + g)/(r - g), 2)
    else
      value = round_to(benchmark_pe*earnings_n, 2)
    end if
  end function terminal_value_using_pe

  pure real(dp) function computing_bv_per_share(total_equity, preferred_stock_mv, outstanding_common_shares) result(value)
    real(dp), intent(in) :: total_equity, preferred_stock_mv, outstanding_common_shares
    value = round_to((total_equity - preferred_stock_mv)/outstanding_common_shares, 2)
  end function computing_bv_per_share

  pure real(dp) function computing_pb(mode, bv0, current_share_price, roe, g, r) result(value)
    character(len=*), intent(in) :: mode
    real(dp), intent(in) :: bv0, current_share_price, roe, g, r
    if (trim(lowercase(mode)) == 'ggm') then
      value = round_to((roe - g)/(r - g), 2)
    else
      value = round_to(current_share_price/bv0, 2)
    end if
  end function computing_pb

  pure real(dp) function computing_ps(mode, current_share_price, payout, eps0, sales0, g, r) result(value)
    character(len=*), intent(in) :: mode
    real(dp), intent(in) :: current_share_price, payout, eps0, sales0, g, r
    if (trim(lowercase(mode)) == 'ggm') then
      value = round_to((eps0/sales0)*payout*(1.0_dp + g)/(r - g), 1)
    else
      value = round_to(current_share_price/sales0, 2)
    end if
  end function computing_ps

  pure real(dp) function computing_ev_dollar_val(common_equity_mv, preferred_stock_mv, debt_mv, cash_and_equivalents) result(value)
    real(dp), intent(in) :: common_equity_mv, preferred_stock_mv, debt_mv, cash_and_equivalents
    value = round_to(common_equity_mv + preferred_stock_mv + debt_mv - cash_and_equivalents, 2)
  end function computing_ev_dollar_val

  pure real(dp) function computing_ev_multiple(basis, ev, ebitda, sales) result(value)
    character(len=*), intent(in) :: basis
    real(dp), intent(in) :: ev, ebitda, sales
    if (trim(lowercase(basis)) == 'ebitda') then
      value = round_to(ev/ebitda, 1)
    else
      value = round_to(ev/sales, 2)
    end if
  end function computing_ev_multiple

  pure real(dp) function computing_sustainable_g(retention_rate, roe) result(value)
    real(dp), intent(in) :: retention_rate, roe
    value = round_to(retention_rate*roe, 4)
  end function computing_sustainable_g

  pure real(dp) function computing_r_with_capm(rfr, market_beta, erp) result(value)
    real(dp), intent(in) :: rfr, market_beta, erp
    value = round_to(rfr + market_beta*erp, 4)
  end function computing_r_with_capm

  pure real(dp) function computing_wacc(dollar_value_debt, dollar_value_common_equity, &
      r_debt, r_common_equity, tax_rate) result(value)
    real(dp), intent(in) :: dollar_value_debt, dollar_value_common_equity
    real(dp), intent(in) :: r_debt, r_common_equity, tax_rate
    real(dp) :: total_capital
    total_capital = dollar_value_debt + dollar_value_common_equity
    value = round_to((dollar_value_debt/total_capital)*r_debt*(1.0_dp - tax_rate) + &
      (dollar_value_common_equity/total_capital)*r_common_equity, 5)
  end function computing_wacc

  pure real(dp) function computing_r_with_ffm(rfr, market_beta, size_beta, value_beta, rmrf, smb, hml) result(value)
    real(dp), intent(in) :: rfr, market_beta, size_beta, value_beta, rmrf, smb, hml
    value = round_to(rfr + market_beta*rmrf + size_beta*smb + value_beta*hml, 3)
  end function computing_r_with_ffm

  pure real(dp) function computing_r_with_hmodel(div_not, sp_not, n, h, g_s, g_l) result(value)
    real(dp), intent(in) :: div_not, sp_not, n, h, g_s, g_l
    if (.not. ieee_is_finite(n)) then
      value = nan_value()
    else
      value = round_to((div_not/sp_not)*((1.0_dp + g_l) + h*(g_s - g_l)) + g_l, 3)
    end if
  end function computing_r_with_hmodel

end module stock_analyst
