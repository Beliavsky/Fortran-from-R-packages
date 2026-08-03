! SPDX-License-Identifier: GPL-3.0-or-later
program identification_tests
  use intrinsicfrp, only: dp, i8, rank_test_result, pca_result
  use intrinsicfrp, only: iterative_kleibergen_paap_2006_beta_rank_test
  use intrinsicfrp, only: chen_fang_2019_beta_rank_test
  use intrinsicfrp, only: giglio_xiu_2021_risk_premia
  implicit none
  integer, parameter :: n = 200, n_assets = 9, n_factors = 3
  real(dp) :: factors(n, n_factors), returns(n, n_assets)
  real(dp) :: loadings(n_assets, n_factors)
  type(rank_test_result) :: kp, cf
  type(pca_result) :: gx
  integer :: i, j

  do i = 1, n
    factors(i, 1) = sin(0.051_dp * real(i, dp))
    factors(i, 2) = cos(0.079_dp * real(i, dp))
    factors(i, 3) = sin(0.137_dp * real(i, dp))
  end do
  do j = 1, n_assets
    loadings(j, 1) = 0.45_dp + 0.03_dp * real(j, dp)
    loadings(j, 2) = -0.30_dp + 0.04_dp * real(j, dp)
    loadings(j, 3) = 0.16_dp * (-1.0_dp) ** j
  end do
  do i = 1, n
    returns(i, :) = matmul(loadings, factors(i, :)) + &
      0.012_dp * [(sin(0.23_dp * real(i + j, dp)), j = 1, n_assets)]
  end do

  call iterative_kleibergen_paap_2006_beta_rank_test(returns, factors, kp)
  write(*, '(a,i0)') 'estimated factor-loading rank: ', kp%rank
  write(*, '(a,*(1x,f9.5))') 'rank-test p-values:', kp%p_values

  call chen_fang_2019_beta_rank_test(returns, factors, cf, &
    n_bootstrap=200, seed=12345_i8)
  write(*, '(a,f9.5)') 'Chen-Fang bootstrap p-value: ', cf%p_value

  call giglio_xiu_2021_risk_premia(returns, factors, gx, which_n_pca=3)
  write(*, '(a,*(1x,f10.6))') 'Giglio-Xiu risk premia:', gx%risk_premia
end program identification_tests
