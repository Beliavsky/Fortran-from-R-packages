! SPDX-License-Identifier: GPL-3.0-or-later
program basic_risk_premia
  use intrinsicfrp, only: dp, vector_result, tfrp, frp, sdf_coefficients
  implicit none
  integer, parameter :: n = 240, n_assets = 8, n_factors = 3
  real(dp) :: factors(n, n_factors), returns(n, n_assets)
  real(dp) :: loadings(n_assets, n_factors)
  type(vector_result) :: fit
  integer :: i, j

  call make_sample(factors, returns, loadings)

  call tfrp(returns, factors, fit, include_standard_errors=.true.)
  write(*, '(a,*(1x,f10.6))') 'TFRP:', fit%estimate
  write(*, '(a,*(1x,f10.6))') 'SE:  ', fit%standard_errors

  call frp(returns, factors, fit, misspecification_robust=.true., &
    include_standard_errors=.true.)
  write(*, '(a,*(1x,f10.6))') 'KRS FRP:', fit%estimate

  call sdf_coefficients(returns, factors, fit, &
    misspecification_robust=.true., include_standard_errors=.true.)
  write(*, '(a,*(1x,f10.6))') 'GKR SDF:', fit%estimate

contains

  subroutine make_sample(f, r, b)
    real(dp), intent(out) :: f(:, :), r(:, :), b(:, :)
    do i = 1, size(f, 1)
      f(i, 1) = 0.010_dp + 0.035_dp * sin(0.061_dp * real(i, dp))
      f(i, 2) = 0.006_dp + 0.028_dp * cos(0.043_dp * real(i, dp))
      f(i, 3) = 0.003_dp + 0.020_dp * sin(0.097_dp * real(i, dp))
    end do
    do j = 1, size(r, 2)
      b(j, 1) = 0.35_dp + 0.04_dp * real(j, dp)
      b(j, 2) = -0.20_dp + 0.03_dp * real(j, dp)
      b(j, 3) = 0.12_dp * real(mod(j, 3) - 1, dp)
    end do
    do i = 1, size(r, 1)
      r(i, :) = 0.001_dp + matmul(b, f(i, :)) + &
        0.012_dp * [(sin(0.17_dp * real(i + 2 * j, dp)), j = 1, size(r, 2))]
    end do
  end subroutine make_sample
end program basic_risk_premia
