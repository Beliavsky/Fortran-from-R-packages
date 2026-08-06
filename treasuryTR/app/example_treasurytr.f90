! SPDX-License-Identifier: MIT
! Copyright (c) 2021 treasuryTR authors

program example_treasurytr
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use treasurytr, only : convexity, dp, mod_duration, prepare_yields, total_return
  implicit none

  real(dp) :: raw_percent(8), yields(8), returns(8)
  real(dp) :: wealth
  integer :: i

  raw_percent = [4.00_dp, 4.05_dp, 4.03_dp, 4.10_dp, 4.08_dp, 4.02_dp, 3.98_dp, 4.01_dp]
  yields = prepare_yields(raw_percent)
  returns = total_return(yields, maturity=10.0_dp, scale=12.0_dp, &
    source_compatible=.false.)

  print '(a)', '10-year constant-maturity Treasury return approximation'
  print '(a,f10.6)', 'initial modified duration: ', mod_duration(yields(1), 10.0_dp, .false.)
  print '(a,f10.6)', 'initial convexity:         ', convexity(yields(1), 10.0_dp, .false.)
  print '(a)', ' period      yield       return'
  do i = 1, size(yields)
    if (ieee_is_nan(returns(i))) then
      print '(i7,2x,f10.6,2x,a)', i, yields(i), 'NA'
    else
      print '(i7,2x,f10.6,2x,f10.6)', i, yields(i), returns(i)
    end if
  end do

  wealth = 1.0_dp
  do i = 2, size(returns)
    wealth = wealth * (1.0_dp + returns(i))
  end do
  print '(a,f10.6)', 'cumulative growth of $1: ', wealth
end program example_treasurytr
