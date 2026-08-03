! SPDX-License-Identifier: GPL-3.0-only AND LicenseRef-ISDA-CDS-Standard-Model
module creditr_cds
  use creditr_kinds, only: dp, creditr_ok, creditr_invalid_input, creditr_no_bracket, creditr_max_iter
  use creditr_dates, only: date_t, cds_dates_t, add_months, add_days, following, year_fraction, date_from_serial
  use creditr_curve, only: zero_curve_t
  implicit none
  private

  type, public :: cds_leg_values_t
    real(kind=dp) :: protection = 0.0_dp
    real(kind=dp) :: premium_dirty = 0.0_dp
    real(kind=dp) :: premium_clean = 0.0_dp
    real(kind=dp) :: accrued = 0.0_dp
    real(kind=dp) :: hazard_rate = 0.0_dp
  end type cds_leg_values_t

  public :: bootstrap_hazard, price_cds_points, spread_to_upfront_curve
  public :: upfront_to_spread_curve, cds_leg_values

contains

  subroutine cds_leg_values(curve, dates, coupon_rate, recovery, hazard, values, status)
    type(zero_curve_t), intent(in) :: curve
    type(cds_dates_t), intent(in) :: dates
    real(kind=dp), intent(in) :: coupon_rate, recovery, hazard
    type(cds_leg_values_t), intent(out) :: values
    integer, intent(out), optional :: status
    real(kind=dp) :: unit_dirty, unit_clean, accrued

    if (coupon_rate < 0.0_dp .or. recovery < 0.0_dp .or. recovery >= 1.0_dp .or. hazard < 0.0_dp) then
      if (present(status)) status = creditr_invalid_input
      return
    end if
    call fee_leg_unit(curve, dates, hazard, unit_dirty, unit_clean, accrued)
    values%protection = protection_leg(curve, dates, hazard, recovery)
    values%premium_dirty = coupon_rate * unit_dirty
    values%premium_clean = coupon_rate * unit_clean
    values%accrued = coupon_rate * accrued
    values%hazard_rate = hazard
    if (present(status)) status = creditr_ok
  end subroutine cds_leg_values

  subroutine bootstrap_hazard(curve, dates, par_spread, recovery, hazard, status)
    type(zero_curve_t), intent(in) :: curve
    type(cds_dates_t), intent(in) :: dates
    real(kind=dp), intent(in) :: par_spread, recovery
    real(kind=dp), intent(out) :: hazard
    integer, intent(out), optional :: status
    real(kind=dp) :: lo, hi, flo, fhi, mid, fmid
    integer :: iter

    if (par_spread < 0.0_dp .or. recovery < 0.0_dp .or. recovery >= 1.0_dp) then
      if (present(status)) status = creditr_invalid_input
      hazard = 0.0_dp
      return
    end if
    if (abs(par_spread) <= tiny(1.0_dp)) then
      hazard = 0.0_dp
      if (present(status)) status = creditr_ok
      return
    end if
    lo = 0.0_dp
    hi = max(1.0_dp, 4.0_dp * par_spread / max(1.0e-8_dp, 1.0_dp - recovery))
    flo = par_equation(lo, curve, dates, par_spread, recovery)
    fhi = par_equation(hi, curve, dates, par_spread, recovery)
    do while (flo * fhi > 0.0_dp .and. hi < 1000.0_dp)
      hi = 2.0_dp * hi
      fhi = par_equation(hi, curve, dates, par_spread, recovery)
    end do
    if (flo * fhi > 0.0_dp) then
      hazard = 0.0_dp
      if (present(status)) status = creditr_no_bracket
      return
    end if
    do iter = 1, 250
      mid = 0.5_dp * (lo + hi)
      fmid = par_equation(mid, curve, dates, par_spread, recovery)
      if (abs(fmid) < 1.0e-14_dp .or. abs(hi - lo) < 1.0e-13_dp) exit
      if (flo * fmid <= 0.0_dp) then
        hi = mid
        fhi = fmid
      else
        lo = mid
        flo = fmid
      end if
    end do
    hazard = 0.5_dp * (lo + hi)
    if (present(status)) status = merge(creditr_ok, creditr_max_iter, iter < 250)
  end subroutine bootstrap_hazard

  real(kind=dp) function par_equation(hazard, curve, dates, spread, recovery) result(value)
    real(kind=dp), intent(in) :: hazard, spread, recovery
    type(zero_curve_t), intent(in) :: curve
    type(cds_dates_t), intent(in) :: dates
    real(kind=dp) :: dirty, clean, accrued
    call fee_leg_unit(curve, dates, hazard, dirty, clean, accrued)
    value = protection_leg(curve, dates, hazard, recovery) - spread * clean
  end function par_equation

  subroutine price_cds_points(curve, dates, par_spread_bps, coupon_bps, recovery, clean, points, hazard, status)
    type(zero_curve_t), intent(in) :: curve
    type(cds_dates_t), intent(in) :: dates
    real(kind=dp), intent(in) :: par_spread_bps, coupon_bps, recovery
    logical, intent(in) :: clean
    real(kind=dp), intent(out) :: points
    real(kind=dp), intent(out), optional :: hazard
    integer, intent(out), optional :: status
    real(kind=dp) :: h, dirty, clean_fee, accrued, fee, protection
    integer :: local_status

    call bootstrap_hazard(curve, dates, par_spread_bps / 1.0e4_dp, recovery, h, local_status)
    if (local_status /= creditr_ok) then
      points = 0.0_dp
      if (present(hazard)) hazard = h
      if (present(status)) status = local_status
      return
    end if
    call fee_leg_unit(curve, dates, h, dirty, clean_fee, accrued)
    protection = protection_leg(curve, dates, h, recovery)
    if (clean) then
      fee = clean_fee
    else
      fee = dirty
    end if
    points = protection - coupon_bps / 1.0e4_dp * fee
    if (present(hazard)) hazard = h
    if (present(status)) status = creditr_ok
  end subroutine price_cds_points

  real(kind=dp) function spread_to_upfront_curve(curve, dates, spread_bps, coupon_bps, recovery, notional, &
      clean, status) result(upfront)
    type(zero_curve_t), intent(in) :: curve
    type(cds_dates_t), intent(in) :: dates
    real(kind=dp), intent(in) :: spread_bps, coupon_bps, recovery, notional
    logical, intent(in), optional :: clean
    integer, intent(out), optional :: status
    logical :: is_clean
    real(kind=dp) :: points
    integer :: local_status
    is_clean = .false.
    if (present(clean)) is_clean = clean
    call price_cds_points(curve, dates, spread_bps, coupon_bps, recovery, is_clean, points, status=local_status)
    upfront = notional * points
    if (present(status)) status = local_status
  end function spread_to_upfront_curve

  real(kind=dp) function upfront_to_spread_curve(curve, dates, upfront_points, coupon_bps, recovery, clean, &
      status) result(spread_bps)
    type(zero_curve_t), intent(in) :: curve
    type(cds_dates_t), intent(in) :: dates
    real(kind=dp), intent(in) :: upfront_points, coupon_bps, recovery
    logical, intent(in), optional :: clean
    integer, intent(out), optional :: status
    real(kind=dp) :: lo, hi, flo, fhi, mid, fmid, value
    logical :: is_clean
    integer :: iter, local_status
    is_clean = .false.
    if (present(clean)) is_clean = clean
    lo = 0.0_dp
    hi = 1000.0_dp
    call price_cds_points(curve, dates, lo, coupon_bps, recovery, is_clean, value, status=local_status)
    flo = value - upfront_points
    call price_cds_points(curve, dates, hi, coupon_bps, recovery, is_clean, value, status=local_status)
    fhi = value - upfront_points
    do while (flo * fhi > 0.0_dp .and. hi < 1.0e7_dp)
      hi = 2.0_dp * hi
      call price_cds_points(curve, dates, hi, coupon_bps, recovery, is_clean, value, status=local_status)
      fhi = value - upfront_points
    end do
    if (flo * fhi > 0.0_dp) then
      spread_bps = 0.0_dp
      if (present(status)) status = creditr_no_bracket
      return
    end if
    do iter = 1, 160
      mid = 0.5_dp * (lo + hi)
      call price_cds_points(curve, dates, mid, coupon_bps, recovery, is_clean, value, status=local_status)
      fmid = value - upfront_points
      if (abs(fmid) < 1.0e-13_dp .or. abs(hi - lo) < 1.0e-9_dp) exit
      if (flo * fmid <= 0.0_dp) then
        hi = mid
        fhi = fmid
      else
        lo = mid
        flo = fmid
      end if
    end do
    spread_bps = 0.5_dp * (lo + hi)
    if (present(status)) status = merge(creditr_ok, creditr_max_iter, iter < 160)
  end function upfront_to_spread_curve

  subroutine fee_leg_unit(curve, dates, hazard, dirty_pv, clean_pv, accrued_fraction)
    type(zero_curve_t), intent(in) :: curve
    type(cds_dates_t), intent(in) :: dates
    real(kind=dp), intent(in) :: hazard
    real(kind=dp), intent(out) :: dirty_pv, clean_pv, accrued_fraction
    integer :: nper, i
    type(date_t) :: unadjusted_start, unadjusted_end, acc_start, acc_end, pay_date
    real(kind=dp) :: regular, aod, value_forward, accrual

    nper = max(1, (12 * (dates%end_date%year - dates%start_date%year) + &
      dates%end_date%month - dates%start_date%month) / 3)
    regular = 0.0_dp
    aod = 0.0_dp
    accrued_fraction = 0.0_dp
    unadjusted_start = dates%start_date
    acc_start = dates%start_date

    do i = 1, nper
      unadjusted_end = add_months(unadjusted_start, 3)
      if (i == nper) unadjusted_end = dates%end_date
      pay_date = following(unadjusted_end)
      if (i == nper) then
        acc_end = add_days(unadjusted_end, 1)
      else
        acc_end = pay_date
      end if
      accrual = year_fraction(acc_start, acc_end, 'ACT/360')
      regular = regular + accrual * survival_probability(dates%trade_date, add_days(acc_end, -1), hazard) * &
        curve%discount(pay_date)
      aod = aod + accrual_on_default_period(curve, dates%trade_date, dates%stepin_date, &
        add_days(acc_start, -1), add_days(acc_end, -1), accrual, hazard)
      if (dates%stepin_date%serial() > acc_start%serial() .and. &
          dates%stepin_date%serial() <= acc_end%serial()) then
        accrued_fraction = year_fraction(acc_start, dates%stepin_date, 'ACT/360')
      end if
      unadjusted_start = unadjusted_end
      acc_start = pay_date
    end do

    value_forward = curve%forward_discount(dates%trade_date, dates%value_date)
    dirty_pv = (regular + aod) / value_forward
    clean_pv = dirty_pv - accrued_fraction
  end subroutine fee_leg_unit

  real(kind=dp) function protection_leg(curve, dates, hazard, recovery) result(pv)
    type(zero_curve_t), intent(in) :: curve
    type(cds_dates_t), intent(in) :: dates
    real(kind=dp), intent(in) :: hazard, recovery
    integer, allocatable :: timeline(:)
    integer :: n, i, s0, s1
    real(kind=dp) :: surv0, surv1, df0, df1, lambda, fwd, x, p0, value_forward

    call make_timeline(curve, dates%trade_date%serial(), dates%end_date%serial(), timeline, n)
    pv = 0.0_dp
    s0 = dates%trade_date%serial()
    do i = 1, n
      s1 = timeline(i)
      surv0 = survival_serial(dates%trade_date%serial(), s0, hazard)
      surv1 = survival_serial(dates%trade_date%serial(), s1, hazard)
      df0 = curve%discount_serial(s0)
      df1 = curve%discount_serial(s1)
      lambda = log(max(surv0, tiny(1.0_dp))) - log(max(surv1, tiny(1.0_dp)))
      fwd = log(max(df0, tiny(1.0_dp))) - log(max(df1, tiny(1.0_dp)))
      x = lambda + fwd + 1.0e-50_dp
      if (abs(x) > 1.0e-4_dp) then
        pv = pv + (1.0_dp - recovery) * lambda / x * (1.0_dp - exp(-x)) * surv0 * df0
      else
        p0 = (1.0_dp - recovery) * lambda * surv0 * df0
        pv = pv + p0 * (1.0_dp - x / 2.0_dp + x * x / 6.0_dp - x**3 / 24.0_dp + x**4 / 120.0_dp)
      end if
      s0 = s1
    end do
    value_forward = curve%forward_discount(dates%trade_date, dates%value_date)
    pv = pv / value_forward
  end function protection_leg

  real(kind=dp) function accrual_on_default_period(curve, today, stepin, start_minus_one, end_minus_one, &
      total_accrual, hazard) result(pv)
    type(zero_curve_t), intent(in) :: curve
    type(date_t), intent(in) :: today, stepin, start_minus_one, end_minus_one
    real(kind=dp), intent(in) :: total_accrual, hazard
    integer, allocatable :: timeline(:)
    integer :: n, i, sub_start, s0, s1, start_serial
    real(kind=dp) :: total_t, acc_rate, t0, t1, dt, surv0, surv1, df0, df1
    real(kind=dp) :: lambda, fwd, x, ratio, la, laf, laf2, laf3, laf4

    start_serial = start_minus_one%serial()
    sub_start = max(stepin%serial() - 1, start_serial)
    if (sub_start >= end_minus_one%serial()) then
      pv = 0.0_dp
      return
    end if
    total_t = real(end_minus_one%serial() - start_serial, dp) / 365.0_dp
    if (total_t <= 0.0_dp) then
      pv = 0.0_dp
      return
    end if
    acc_rate = total_accrual / total_t
    call make_timeline(curve, sub_start, end_minus_one%serial(), timeline, n)
    pv = 0.0_dp
    s0 = sub_start
    do i = 1, n
      s1 = timeline(i)
      t0 = real(s0 + 1 - start_serial, dp) / 365.0_dp - 0.5_dp / 365.0_dp
      t1 = real(s1 + 1 - start_serial, dp) / 365.0_dp - 0.5_dp / 365.0_dp
      dt = t1 - t0
      surv0 = survival_serial(today%serial(), s0, hazard)
      surv1 = survival_serial(today%serial(), s1, hazard)
      df0 = curve%discount_serial(s0)
      df1 = curve%discount_serial(s1)
      lambda = log(max(surv0, tiny(1.0_dp))) - log(max(surv1, tiny(1.0_dp)))
      fwd = log(max(df0, tiny(1.0_dp))) - log(max(df1, tiny(1.0_dp)))
      x = lambda + fwd + 1.0e-50_dp
      if (abs(x) > 1.0e-4_dp) then
        ratio = (surv1 / surv0) * (df1 / df0)
        pv = pv + lambda * acc_rate * surv0 * df0 * &
          ((t0 + dt / x) / x - (t1 + dt / x) / x * ratio)
      else
        la = lambda * surv0 * df0 * acc_rate * 0.5_dp
        pv = pv + la * (t0 + t1)
        laf = la * x / 3.0_dp
        pv = pv - laf * (t0 + 2.0_dp * t1)
        laf2 = laf * x * 0.25_dp
        pv = pv + laf2 * (t0 + 3.0_dp * t1)
        laf3 = laf2 * x * 0.2_dp
        pv = pv - laf3 * (t0 + 4.0_dp * t1)
        laf4 = laf3 * x / 6.0_dp
        pv = pv + laf4 * (t0 + 5.0_dp * t1)
      end if
      s0 = s1
    end do
  end function accrual_on_default_period

  subroutine make_timeline(curve, start_serial, end_serial, timeline, n)
    type(zero_curve_t), intent(in) :: curve
    integer, intent(in) :: start_serial, end_serial
    integer, allocatable, intent(out) :: timeline(:)
    integer, intent(out) :: n
    integer, allocatable :: temp(:)
    integer :: i
    allocate(temp(size(curve%node_serial) + 1))
    n = 0
    do i = 1, size(curve%node_serial)
      if (curve%node_serial(i) > start_serial .and. curve%node_serial(i) < end_serial) then
        n = n + 1
        temp(n) = curve%node_serial(i)
      end if
    end do
    n = n + 1
    temp(n) = end_serial
    allocate(timeline(n))
    timeline = temp(1:n)
  end subroutine make_timeline

  pure real(kind=dp) function survival_probability(today, dt, hazard) result(s)
    type(date_t), intent(in) :: today, dt
    real(kind=dp), intent(in) :: hazard
    s = survival_serial(today%serial(), dt%serial(), hazard)
  end function survival_probability

  pure real(kind=dp) function survival_serial(today_serial, serial, hazard) result(s)
    integer, intent(in) :: today_serial, serial
    real(kind=dp), intent(in) :: hazard
    s = exp(-hazard * real(max(0, serial - today_serial), dp) / 365.0_dp)
  end function survival_serial

end module creditr_cds
