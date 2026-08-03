program test_bootstrap_bias
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use multiatsm_kinds, only : dp
  use multiatsm_linalg, only : spectral_radius
  use multiatsm_types, only : bootstrap_result
  use multiatsm_bootstrap, only : BOOTSTRAP_IID, BOOTSTRAP_BLOCK, resample_residuals, &
    bootstrap_var, percentile_bounds
  use multiatsm_bias, only : bias_correct_var, shrink_transition
  implicit none
  real(dp) :: data(2, 120), factor_residuals(2, 20), yield_residuals(1, 20)
  real(dp) :: draws(1, 1, 5), unstable(2, 2), target(2, 2), radius
  real(dp), allocatable :: f1(:, :), y1(:, :), f2(:, :), y2(:, :)
  real(dp), allocatable :: lower(:, :), median(:, :), upper(:, :)
  real(dp), allocatable :: corrected(:, :), corrected_sigma(:, :), shrunk(:, :)
  type(bootstrap_result) :: result
  integer :: t, info

  do t = 1, 20
    factor_residuals(:, t) = [real(t, dp), -real(t, dp)]
    yield_residuals(1, t) = 0.5_dp * real(t, dp)
  end do
  call resample_residuals(factor_residuals, yield_residuals, BOOTSTRAP_IID, f1, y1, info, seed=777)
  call check(info == 0, 'IID resampling status')
  call resample_residuals(factor_residuals, yield_residuals, BOOTSTRAP_IID, f2, y2, info, seed=777)
  call check(maxval(abs(f1 - f2)) < 1.0e-12_dp .and. maxval(abs(y1 - y2)) < 1.0e-12_dp, &
    'bootstrap reproducibility')
  call resample_residuals(factor_residuals, yield_residuals, BOOTSTRAP_BLOCK, f2, y2, info, &
    block_length=4, seed=123)
  call check(info == 0, 'block resampling status')

  draws(1, 1, :) = [5.0_dp, 1.0_dp, 4.0_dp, 2.0_dp, 3.0_dp]
  call percentile_bounds(draws, 0.0_dp, 1.0_dp, lower, median, upper, info)
  call check(info == 0, 'percentile status')
  call check(abs(lower(1, 1) - 1.0_dp) < 1.0e-12_dp .and. abs(median(1, 1) - 3.0_dp) < 1.0e-12_dp .and. &
    abs(upper(1, 1) - 5.0_dp) < 1.0e-12_dp, &
    'percentile values')

  data(:, 1) = [0.1_dp, -0.05_dp]
  do t = 2, size(data, 2)
    data(1, t) = 0.02_dp + 0.75_dp * data(1, t - 1) + 0.08_dp * data(2, t - 1) + &
      0.015_dp * sin(0.31_dp * real(t, dp))
    data(2, t) = -0.01_dp + 0.05_dp * data(1, t - 1) + 0.55_dp * data(2, t - 1) + &
      0.012_dp * cos(0.27_dp * real(t, dp))
  end do
  call bootstrap_var(data, 12, BOOTSTRAP_IID, result, info, seed=2026)
  call check(info == 0, 'VAR bootstrap status')
  call check(size(result%phi_draws, 3) == 12 .and. all(ieee_is_finite(result%phi_draws)), &
    'VAR bootstrap draws')
  call bias_correct_var(data, 4, 1, 2, 0.2_dp, corrected, corrected_sigma, info, seed=31337)
  call check(info == 0, 'bias correction status')
  call check(all(ieee_is_finite(corrected)) .and. all(ieee_is_finite(corrected_sigma)), &
    'finite bias correction')

  unstable = 0.0_dp
  unstable(1, 1) = 1.2_dp
  unstable(2, 2) = 0.7_dp
  target = 0.0_dp
  target(1, 1) = 0.8_dp
  target(2, 2) = 0.6_dp
  call shrink_transition(unstable, target, 0.95_dp, shrunk, info)
  call check(info == 0, 'shrink status')
  radius = spectral_radius(shrunk, info)
  call check(info == 0 .and. radius <= 0.950001_dp, 'shrink radius')
  print '(a)', 'test_bootstrap_bias: PASS'
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // message
      error stop 1
    end if
  end subroutine check
end program test_bootstrap_bias
