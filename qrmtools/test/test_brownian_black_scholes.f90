! SPDX-License-Identifier: GPL-3.0-or-later
program test_brownian_black_scholes
  use qrmtools, only : dp, normal_cdf, r_brownian, de_browning, &
    black_scholes, black_scholes_greeks, brownian_result, greeks_result, &
    compute_returns, invert_returns, return_logarithmic, return_simple, &
    return_difference
  implicit none

  type(brownian_result) :: simulated
  type(greeks_result) :: greeks
  real(dp) :: times(3)
  real(dp) :: uniforms(2,1)
  real(dp), allocatable :: increments(:,:,:)
  real(dp), allocatable :: returns(:,:)
  real(dp), allocatable :: recovered(:,:)
  real(dp) :: prices(3,2)

  call assert_close(black_scholes(0.0_dp,100.0_dp,0.05_dp,0.2_dp, &
    100.0_dp,1.0_dp), 10.450583572185565_dp, 3.0e-11_dp)
  call assert_close(black_scholes(0.0_dp,100.0_dp,0.05_dp,0.2_dp, &
    100.0_dp,1.0_dp,.true.), 5.573526022256971_dp, 3.0e-11_dp)

  greeks = black_scholes_greeks(0.0_dp,100.0_dp,0.05_dp,0.2_dp, &
    100.0_dp,1.0_dp)
  call assert_true(greeks%ok)
  call assert_close(greeks%delta, 0.6368306511756191_dp, 3.0e-11_dp)
  call assert_close(greeks%gamma, 0.018762017345846895_dp, 3.0e-11_dp)
  call assert_close(greeks%vega, 37.52403469169379_dp, 3.0e-10_dp)

  times = [0.0_dp,0.5_dp,1.0_dp]
  uniforms(:,1) = [normal_cdf(1.0_dp), normal_cdf(-0.5_dp)]
  simulated = r_brownian(1,times,d=1,drift=[0.1_dp], &
    volatility=[0.2_dp],process_type='BM',u=uniforms)
  call assert_true(simulated%ok)
  call assert_close(simulated%paths(1,1,1),0.0_dp,1.0e-14_dp)
  call assert_close(simulated%paths(1,2,1),0.19142135623730952_dp,2.0e-10_dp)
  call assert_close(simulated%paths(1,3,1),0.17071067811865477_dp,2.0e-10_dp)

  increments = de_browning(simulated%paths,times,drift=[0.1_dp], &
    volatility=[0.2_dp],process_type='BM')
  call assert_close(increments(1,1,1),1.0_dp,2.0e-9_dp)
  call assert_close(increments(1,2,1),-0.5_dp,2.0e-9_dp)

  simulated = r_brownian(1,times,d=1,drift=[0.1_dp], &
    volatility=[0.2_dp],initial=[100.0_dp],process_type='GBM',u=uniforms)
  call assert_true(simulated%ok)
  call assert_close(simulated%paths(1,3,1),100.0_dp*exp(0.17071067811865477_dp), &
    3.0e-8_dp)

  prices = reshape([100.0_dp,105.0_dp,110.0_dp,50.0_dp,48.0_dp,52.0_dp],[3,2])
  returns = compute_returns(prices,return_logarithmic)
  recovered = invert_returns(returns,prices(1,:),return_logarithmic)
  call assert_array_close(recovered,prices,2.0e-12_dp)
  returns = compute_returns(prices,return_simple)
  recovered = invert_returns(returns,prices(1,:),return_simple)
  call assert_array_close(recovered,prices,2.0e-12_dp)
  returns = compute_returns(prices,return_difference)
  recovered = invert_returns(returns,prices(1,:),return_difference)
  call assert_array_close(recovered,prices,2.0e-12_dp)

  print '(a)', 'test_brownian_black_scholes: PASS'

contains

  subroutine assert_close(actual,expected,tolerance)
    real(dp), intent(in) :: actual,expected,tolerance
    if(abs(actual-expected)>tolerance*max(1.0_dp,abs(expected))) then
      print *, 'mismatch:',actual,expected,abs(actual-expected)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_array_close(actual,expected,tolerance)
    real(dp), intent(in) :: actual(:,:),expected(:,:),tolerance
    if(maxval(abs(actual-expected))>tolerance*max(1.0_dp,maxval(abs(expected)))) then
      error stop 1
    end if
  end subroutine assert_array_close

  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if(.not.condition) error stop 1
  end subroutine assert_true

end program test_brownian_black_scholes
