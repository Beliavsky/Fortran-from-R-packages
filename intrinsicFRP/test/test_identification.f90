! SPDX-License-Identifier: GPL-3.0-or-later
program test_identification
  use intrinsicfrp, only: dp, i8, rank_test_result, pca_result
  use intrinsicfrp, only: iterative_kleibergen_paap_2006_beta_rank_test
  use intrinsicfrp, only: chen_fang_2019_beta_rank_test
  use intrinsicfrp, only: giglio_xiu_2021_risk_premia, scaled_factor_loadings
  use intrinsicfrp, only: npca_giglio_xiu_2021, npca_ahn_horenstein_2013
  implicit none
  integer, parameter :: n = 180, p = 8, k = 3
  real(dp) :: factors(n, k), returns(n, p), loading(p, k)
  real(dp), allocatable :: theta(:, :)
  real(dp) :: evals(6)
  type(rank_test_result) :: kp, cf
  type(pca_result) :: gx
  integer :: i, j, st, er, gr, npc

  do i = 1, n
    factors(i, 1) = sin(0.05_dp * real(i, dp))
    factors(i, 2) = cos(0.08_dp * real(i, dp))
    factors(i, 3) = sin(0.13_dp * real(i, dp))
  end do
  do j = 1, p
    loading(j, 1) = 0.5_dp + 0.03_dp * real(j, dp)
    loading(j, 2) = -0.4_dp + 0.04_dp * real(j, dp)
    loading(j, 3) = 0.2_dp * (-1.0_dp) ** j
  end do
  do i = 1, n
    returns(i, :) = matmul(loading, factors(i, :)) + &
      0.01_dp * [(cos(0.23_dp * real(i + j, dp)), j = 1, p)]
  end do

  call scaled_factor_loadings(returns, factors, theta, st)
  call check(st == 0 .and. size(theta, 1) == k .and. size(theta, 2) == p, &
    'scaled loadings')

  call iterative_kleibergen_paap_2006_beta_rank_test(returns, factors, kp, 0.05_dp)
  call check(kp%status == 0, 'kp status')
  call check(kp%rank >= 2 .and. kp%rank <= k, 'kp rank')
  call check(size(kp%statistics) == k .and. all(kp%p_values >= 0.0_dp), 'kp arrays')

  call chen_fang_2019_beta_rank_test(returns, factors, cf, 80, 0.05_dp, 12345_i8)
  call check(cf%status == 0, 'chen fang status')
  call check(cf%p_value >= 0.0_dp .and. cf%p_value <= 1.0_dp, 'chen fang pvalue')

  call giglio_xiu_2021_risk_premia(returns, factors, gx, which_n_pca=3)
  call check(gx%status == 0 .and. gx%n_pca == 3, 'gx fixed pca')
  call check(size(gx%risk_premia) == k, 'gx risk premia length')

  evals = [4.0_dp, 2.0_dp, 0.7_dp, 0.2_dp, 0.08_dp, 0.03_dp]
  npc = npca_giglio_xiu_2021(evals, 20, 100, 5)
  call npca_ahn_horenstein_2013(evals, 20, 100, 5, er, gr)
  call check(npc >= 1 .and. npc <= 5, 'gx selector range')
  call check(er >= 1 .and. er <= 4 .and. gr >= 1 .and. gr <= 4, 'ah selector range')

  print '(a)', 'test_identification: PASS'

contains
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a,1x,a)') 'FAIL:', trim(label)
      error stop 1
    end if
  end subroutine check
end program test_identification
