! SPDX-License-Identifier: GPL-2.0-or-later
program test_options
  use jrvfinance, only: dp, black_scholes_result, gen_bs, gen_bs_implied, JRV_OK
  implicit none
  type(black_scholes_result) :: bs
  real(dp) :: sigma
  integer :: status

  bs = gen_bs(100.0_dp,100.0_dp,0.10_dp,0.20_dp,1.0_dp,0.0_dp)
  call check(bs%status == JRV_OK, 'Black-Scholes status')
  call check(abs(bs%call-13.2696765846609_dp) < 1.0e-10_dp, 'call value')
  call check(abs(bs%put-3.75341838825683_dp) < 1.0e-10_dp, 'put value')
  call check(abs(bs%call-bs%put-(100.0_dp-100.0_dp*exp(-0.10_dp))) < &
    1.0e-11_dp, 'put-call parity')
  call check(bs%gamma > 0.0_dp .and. bs%vega > 0.0_dp, 'Greeks')

  sigma = gen_bs_implied(100.0_dp,100.0_dp,0.10_dp,bs%call,1.0_dp, &
    0.0_dp,.false.,status=status)
  call check(status == JRV_OK .and. abs(sigma-0.20_dp) < 1.0e-7_dp, &
    'call implied volatility')
  sigma = gen_bs_implied(100.0_dp,100.0_dp,0.10_dp,bs%put,1.0_dp, &
    0.0_dp,.true.,status=status)
  call check(status == JRV_OK .and. abs(sigma-0.20_dp) < 1.0e-7_dp, &
    'put implied volatility')

  bs = gen_bs(100.0_dp,110.0_dp,0.05_dp,0.0_dp,1.0_dp,0.0_dp)
  call check(bs%call >= 0.0_dp .and. bs%put >= 0.0_dp, 'zero volatility')
  bs = gen_bs(100.0_dp,90.0_dp,0.05_dp,0.2_dp,0.0_dp,0.0_dp)
  call check(abs(bs%call-10.0_dp) < 1.0e-12_dp, 'expiry payoff')

  print '(a)', 'test_options: PASS'
contains
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//trim(label)
      error stop 1
    end if
  end subroutine check
end program test_options
