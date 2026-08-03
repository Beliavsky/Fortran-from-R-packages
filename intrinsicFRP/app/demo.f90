! SPDX-License-Identifier: GPL-3.0-or-later
program demo
  use intrinsicfrp, only: dp, vector_result, hj_result, tfrp
  use intrinsicfrp, only: hj_misspecification_distance
  implicit none
  integer, parameter :: n = 180, n_assets = 6, n_factors = 2
  real(dp) :: factors(n, n_factors), returns(n, n_assets)
  real(dp) :: beta(n_assets, n_factors)
  type(vector_result) :: fit
  type(hj_result) :: hj
  integer :: i, j

  do i = 1, n
    factors(i, 1) = 0.008_dp + 0.030_dp * sin(0.07_dp * real(i, dp))
    factors(i, 2) = 0.004_dp + 0.020_dp * cos(0.11_dp * real(i, dp))
  end do
  do j = 1, n_assets
    beta(j, 1) = 0.30_dp + 0.08_dp * real(j, dp)
    beta(j, 2) = -0.18_dp + 0.05_dp * real(j, dp)
  end do
  do i = 1, n
    returns(i, :) = matmul(beta, factors(i, :)) + &
      0.01_dp * [(cos(0.13_dp * real(i + j, dp)), j = 1, n_assets)]
  end do

  call tfrp(returns, factors, fit, include_standard_errors=.true.)
  call hj_misspecification_distance(returns, factors, hj)
  write(*, '(a,*(1x,f10.6))') 'tradable factor risk premia:', fit%estimate
  write(*, '(a,*(1x,f10.6))') 'HAC standard errors:', fit%standard_errors
  write(*, '(a,f12.7)') 'squared HJ distance:', hj%squared_distance
  write(*, '(a,2(1x,f12.7))') '95% interval:', hj%lower_bound, hj%upper_bound
end program demo
