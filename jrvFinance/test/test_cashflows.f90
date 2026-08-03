! SPDX-License-Identifier: GPL-2.0-or-later
program test_cashflows
  use jrvfinance, only: dp, JRV_OK, equiv_rate, npv, irr, duration, &
    annuity_pv, annuity_fv, annuity_instalment, annuity_periods, &
    annuity_rate, annuity_instalment_breakup, annuity_breakup_result
  implicit none
  real(dp) :: x, y, rate, payment, periods
  real(dp), parameter :: cf(3) = [100.0_dp, 250.0_dp, 300.0_dp]
  real(dp), parameter :: irr_cf(3) = [-600.0_dp, 300.0_dp, 400.0_dp]
  type(annuity_breakup_result) :: breakup
  integer :: status

  x = equiv_rate(0.10_dp, 12.0_dp, 2.0_dp)
  y = equiv_rate(x, 2.0_dp, 12.0_dp)
  call check(abs(y-0.10_dp) < 1.0e-12_dp, 'equivalent rate round trip')

  x = npv(cf, 0.05_dp)
  y = sum(cf/[1.05_dp, 1.05_dp**2, 1.05_dp**3])
  call check(abs(x-y) < 1.0e-12_dp, 'npv')

  rate = irr(irr_cf, status=status)
  call check(status == JRV_OK, 'irr status')
  call check(abs(npv(irr_cf, rate, cf_t=[0.0_dp,1.0_dp,2.0_dp])) < 1.0e-5_dp, 'irr root')

  x = duration(cf, 0.05_dp)
  call check(x > 2.0_dp .and. x < 3.0_dp, 'duration range')
  call check(duration(cf,0.05_dp,modified=.true.) < x, 'modified duration')

  x = annuity_pv(0.09_dp, 8.0_dp)
  payment = annuity_instalment(0.09_dp, 8.0_dp, pv=x)
  call check(abs(payment-1.0_dp) < 1.0e-11_dp, 'annuity payment round trip')

  y = annuity_fv(0.09_dp, 8.0_dp, payment)
  call check(y > x, 'annuity future value')
  periods = annuity_periods(0.09_dp, payment, x)
  call check(abs(periods-8.0_dp) < 1.0e-10_dp, 'annuity periods')

  rate = annuity_rate(360.0_dp, 450.0_dp, 50000.0_dp, &
    cf_freq=12.0_dp, comp_freq=2.0_dp, status=status)
  call check(status == JRV_OK .and. rate > 0.0_dp, 'annuity rate')
  call check(abs(annuity_pv(rate,360.0_dp,450.0_dp,cf_freq=12.0_dp, &
    comp_freq=2.0_dp)-50000.0_dp) < 0.02_dp, 'annuity rate round trip')

  breakup = annuity_instalment_breakup(0.09_dp, 8.0_dp, 10000.0_dp, period_no=5)
  call check(breakup%status == JRV_OK, 'breakup status')
  call check(abs(breakup%principal_part + breakup%interest_part - &
    annuity_instalment(0.09_dp,8.0_dp,10000.0_dp)) < 1.0e-10_dp, 'breakup identity')

  print '(a)', 'test_cashflows: PASS'
contains
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//trim(label)
      error stop 1
    end if
  end subroutine check
end program test_cashflows
