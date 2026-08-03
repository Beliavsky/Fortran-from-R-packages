! SPDX-License-Identifier: GPL-3.0-or-later
program oracle_and_screening
  use intrinsicfrp, only: dp, oracle_control, oracle_result, screening_result
  use intrinsicfrp, only: oracle_tfrp, gkr_factor_screening
  implicit none
  integer, parameter :: n = 220, n_assets = 10, n_factors = 4
  real(dp) :: factors(n, n_factors), returns(n, n_assets)
  real(dp) :: loadings(n_assets, n_factors)
  real(dp) :: penalties(7)
  type(oracle_control) :: control
  type(oracle_result) :: oracle
  type(screening_result) :: screen
  integer :: i, j

  call make_sample(factors, returns, loadings)
  penalties = [0.0_dp, 0.0005_dp, 0.001_dp, 0.002_dp, &
    0.005_dp, 0.010_dp, 0.025_dp]

  control = oracle_control()
  control%tuning_type = 'g'
  control%weighting_type = 'c'
  control%include_standard_errors = .true.
  call oracle_tfrp(returns, factors, penalties, oracle, control)
  write(*, '(a,f10.6)') 'selected penalty: ', oracle%penalty_parameter
  write(*, '(a,*(1x,f10.6))') 'Oracle TFRP:', oracle%risk_premia

  call gkr_factor_screening(returns, factors, screen, target_level=0.05_dp)
  write(*, '(a,*(1x,i0))') 'GKR selected factor indices:', &
    screen%selected_indices

contains

  subroutine make_sample(f, r, b)
    real(dp), intent(out) :: f(:, :), r(:, :), b(:, :)
    do i = 1, size(f, 1)
      f(i, 1) = 0.030_dp * sin(0.041_dp * real(i, dp))
      f(i, 2) = 0.025_dp * cos(0.063_dp * real(i, dp))
      f(i, 3) = 0.004_dp * sin(0.113_dp * real(i, dp))
      f(i, 4) = 0.003_dp * cos(0.157_dp * real(i, dp))
    end do
    do j = 1, size(r, 2)
      b(j, 1) = 0.40_dp + 0.03_dp * real(j, dp)
      b(j, 2) = -0.25_dp + 0.025_dp * real(j, dp)
      b(j, 3) = merge(0.03_dp, 0.0_dp, mod(j, 2) == 0)
      b(j, 4) = merge(-0.02_dp, 0.0_dp, mod(j, 3) == 0)
    end do
    do i = 1, size(r, 1)
      r(i, :) = 0.002_dp + matmul(b, f(i, :)) + &
        0.015_dp * [(cos(0.19_dp * real(i + j, dp)), j = 1, size(r, 2))]
    end do
  end subroutine make_sample
end program oracle_and_screening
