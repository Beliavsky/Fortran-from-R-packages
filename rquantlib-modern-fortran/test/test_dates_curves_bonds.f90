! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
program test_dates_curves_bonds
  use rq_kinds, only: dp
  use rq_dates
  use rq_curves
  use rq_bonds
  use test_support
  implicit none
  type(date_t)::d,d2,weekend
  type(calendar_t)::cal
  type(schedule_t)::schedule
  type(discount_curve_t)::flat,boot,forward,zero_curve,future_curve
  type(bond_result)::bond,float_bond
  type(curve_fit_result)::ns,sv
  real(dp)::price,yield,times(10),rates(10),zeros(10),fitted(10),pv
  integer::status,i

  d=make_date(2024,2,29)
  d2=serial_to_date(date_to_serial(d))
  call assert_true(d2%year==2024.and.d2%month==2.and.d2%day==29, &
    'Date serial round trip')
  call assert_true(weekday(make_date(1970,1,1))==4,'Weekday reference')
  weekend=make_date(2026,7,25)
  call assert_true(is_weekend(weekend),'Weekend detection')
  allocate(cal%holidays(1)); cal%holidays(1)=make_date(2026,7,27)
  call assert_true(.not.is_business_day(cal%holidays(1),cal),'Holiday exclusion')
  d2=adjust_date(weekend,following,cal)
  call assert_true(d2%year==2026.and.d2%month==7.and.d2%day==28, &
    'Following adjustment across holiday')
  call assert_close(year_fraction(make_date(2024,1,1),make_date(2025,1,1), &
    'ActualActual'),1.0_dp,1.0e-14_dp,'Actual/Actual leap year')
  call assert_true(day_count(make_date(2024,1,31),make_date(2024,2,29), &
    '30/360')==29,'Thirty/360 day count')
  call make_schedule(make_date(2024,1,31),make_date(2025,1,31),3, &
    convention=unadjusted,end_of_month=.true.,schedule=schedule)
  call assert_true(size(schedule%dates)==5,'Quarterly schedule size')
  call assert_true(schedule%dates(2)%day==30,'End-of-month schedule')

  d2=advance_date(make_date(2024,1,31),1,'months',unadjusted, &
    end_of_month=.true.)
  call assert_true(d2%month==2.and.d2%day==29,'End-of-month advance')
  call assert_true(business_days_between(make_date(2024,1,1), &
    make_date(2024,1,8),include_first=.true.,include_last=.false.)==5, &
    'Business days between')

  call make_flat_curve(0.05_dp,30.0_dp,0.25_dp,flat)
  call assert_close(flat%discount(5.0_dp),exp(-0.25_dp),1.0e-13_dp, &
    'Flat curve discount')
  call assert_close(flat%zero_rate(7.0_dp),0.05_dp,1.0e-12_dp, &
    'Flat curve zero rate')

  call make_zero_curve([1.0_dp,2.0_dp,5.0_dp],[0.03_dp,0.04_dp,0.05_dp],zero_curve)
  call assert_close(zero_curve%zero_rate(2.0_dp),0.04_dp,1.0e-12_dp, &
    'Explicit zero curve')
  call bootstrap_curve([0.25_dp],[0.04_dp],future_starts=[0.25_dp,0.5_dp], &
    future_ends=[0.5_dp,0.75_dp],future_prices=[95.5_dp,95.0_dp], &
    curve=future_curve,status=status)
  call assert_true(status==0.and.future_curve%discount(0.75_dp)>0.0_dp, &
    'Futures bootstrap path')
  call bootstrap_curve([0.25_dp,0.5_dp],[0.04_dp,0.042_dp], &
    swap_maturities=[1.0_dp,2.0_dp,3.0_dp], &
    swap_rates=[0.043_dp,0.045_dp,0.047_dp],swap_frequency=2, &
    curve=boot,status=status)
  call assert_true(status==0,'Bootstrap status')
  call assert_true(all(boot%discounts>0.0_dp),'Positive bootstrapped discounts')
  call assert_true(all(boot%discounts(2:)<=boot%discounts(:size(boot%discounts)-1)), &
    'Monotone example discount curve')

  price=zero_price_by_yield(0.1478_dp,100.0_dp, &
    real(date_to_serial(make_date(1993,11,1))- &
    date_to_serial(make_date(1993,6,24)),dp)/365.0_dp)
  call zero_yield_by_price(price,100.0_dp, &
    real(date_to_serial(make_date(1993,11,1))- &
    date_to_serial(make_date(1993,6,24)),dp)/365.0_dp,yield,status)
  call assert_true(status==0,'Zero yield status')
  call assert_close(yield,0.1478_dp,1.0e-12_dp,'Zero price/yield inversion')

  price=zero_price_by_yield(0.06_dp,100.0_dp,2.0_dp,compounding='continuous')
  call zero_yield_by_price(price,100.0_dp,2.0_dp,yield,status, &
    compounding='continuous')
  call assert_close(yield,0.06_dp,1.0e-12_dp,'Continuous zero compounding')
  price=zero_price_by_yield(0.06_dp,100.0_dp,2.0_dp,compounding='simple')
  call zero_yield_by_price(price,100.0_dp,2.0_dp,yield,status, &
    compounding='simple')
  call assert_close(yield,0.06_dp,1.0e-12_dp,'Simple zero compounding')

  call fixed_rate_bond_from_yield(100.0_dp,0.04_dp,5.0_dp,2,0.03_dp,bond)
  call fixed_rate_bond_yield(100.0_dp,0.04_dp,5.0_dp,2,bond%npv, &
    yield,status)
  call assert_true(status==0,'Fixed bond yield status')
  call assert_close(yield,0.03_dp,1.0e-10_dp,'Fixed bond price/yield inversion')
  call assert_true(bond%duration>0.0_dp.and.bond%convexity>0.0_dp, &
    'Bond duration and convexity')

  pv=cashflow_present_value([1.0_dp,2.0_dp],[5.0_dp,105.0_dp],flat)
  call assert_true(pv>0.0_dp,'Cashflow present value')

  call fixed_rate_bond_from_curve(100.0_dp,0.05_dp,4.0_dp,2,flat,bond)
  call assert_close(bond%npv,100.0_dp,0.3_dp,'Par bond on flat coupon curve')

  call make_flat_curve(0.03_dp,10.0_dp,0.25_dp,forward)
  call floating_rate_bond_from_curves(100.0_dp,4.0_dp,2,forward,flat, &
    float_bond,spread=0.002_dp)
  call assert_true(float_bond%npv>0.0_dp.and.size(float_bond%cashflows)==8, &
    'Floating-rate bond')

  call floating_rate_bond_from_curves(100.0_dp,2.0_dp,2,forward,flat, &
    float_bond,gearing=1.5_dp,spread=0.01_dp,cap=0.04_dp,floor=0.02_dp)
  call assert_true(all(float_bond%cashflows(:size(float_bond%cashflows)-1)>=1.0_dp) .and. &
    all(float_bond%cashflows(:size(float_bond%cashflows)-1)<=2.0_dp), &
    'Floating coupon cap and floor')

  do i=1,10
    times(i)=0.5_dp*real(i,dp)
    rates(i)=nelson_siegel_rate(times(i),0.04_dp,-0.02_dp,0.015_dp,1.8_dp)
  end do
  call fit_nelson_siegel(times,rates,ns)
  call assert_true(ns%status==0,'Nelson-Siegel fit status')
  call assert_close(ns%sse,0.0_dp,1.0e-10_dp,'Nelson-Siegel synthetic fit')

  do i=1,10
    zeros(i)=svensson_rate(times(i),0.04_dp,-0.02_dp,0.015_dp, &
      -0.01_dp,1.2_dp,4.0_dp)
  end do
  call fit_svensson(times,zeros,sv)
  fitted=sv%fitted
  call assert_true(sv%status==0,'Svensson fit status')
  call assert_true(sqrt(sum((fitted-zeros)**2)/10.0_dp)<2.0e-4_dp, &
    'Svensson synthetic fit')
  write(*,'(a)') 'Date, curve, and bond tests passed.'
end program test_dates_curves_bonds
