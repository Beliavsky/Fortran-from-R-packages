! SPDX-License-Identifier: GPL-3.0-or-later
program test_oracle_fgx
  use intrinsicfrp, only: dp, oracle_control, oracle_result, fgx_result, vector_result
  use intrinsicfrp, only: oracle_tfrp, oracle_soft_threshold, fgx_factors_test, tfrp
  implicit none
  integer, parameter :: n = 210, p = 10, k = 4, knew = 2
  real(dp) :: factors(n, k), returns(n, p), gross(n, p), newf(n, knew)
  real(dp) :: loading(p, k), penalties(5), zeroed(k)
  type(oracle_control) :: ctl
  type(oracle_result) :: out
  type(fgx_result) :: fgx
  type(vector_result) :: base
  integer :: i, j

  do i = 1, n
    factors(i, 1) = sin(0.031_dp * real(i, dp))
    factors(i, 2) = cos(0.057_dp * real(i, dp))
    factors(i, 3) = sin(0.097_dp * real(i, dp))
    factors(i, 4) = cos(0.149_dp * real(i, dp))
  end do
  do j = 1, p
    loading(j, 1) = 0.2_dp + 0.04_dp * real(j, dp)
    loading(j, 2) = -0.15_dp + 0.02_dp * real(j, dp)
    loading(j, 3) = merge(0.30_dp, 0.0_dp, mod(j, 2) == 0)
    loading(j, 4) = merge(-0.20_dp, 0.0_dp, mod(j, 3) == 0)
  end do
  do i = 1, n
    returns(i, :) = 0.004_dp + matmul(loading, factors(i, :)) + &
      0.02_dp * [(sin(0.21_dp * real(i + j, dp)), j = 1, p)]
  end do
  gross = 1.0_dp + returns
  newf(:, 1) = 0.7_dp * factors(:, 3) + 0.1_dp * sin([(0.33_dp * real(i, dp), i = 1, n)])
  newf(:, 2) = 0.5_dp * factors(:, 4) + 0.1_dp * cos([(0.27_dp * real(i, dp), i = 1, n)])
  penalties = [0.0_dp, 0.001_dp, 0.005_dp, 0.02_dp, 0.08_dp]

  call tfrp(returns, factors, base)
  zeroed = oracle_soft_threshold(base%estimate, spread(1.0_dp, 1, k), 0.0_dp)
  call check(maxval(abs(zeroed - base%estimate)) < 1.0e-12_dp, 'zero penalty')

  ctl = oracle_control()
  ctl%tuning_type = 'g'
  ctl%include_standard_errors = .true.
  call oracle_tfrp(returns, factors, penalties, out, ctl)
  call check(out%status == 0 .and. size(out%risk_premia) == k, 'oracle gcv')
  call check(size(out%model_score) == size(penalties), 'oracle score')

  ctl%tuning_type = 'c'
  ctl%n_folds = 5
  call oracle_tfrp(returns, factors, penalties, out, ctl)
  call check(out%status == 0, 'oracle cv')

  ctl%tuning_type = 'r'
  ctl%n_train_observations = 120
  ctl%n_test_observations = 24
  ctl%roll_shift = 12
  call oracle_tfrp(returns, factors, penalties, out, ctl)
  call check(out%status == 0, 'oracle rolling')

  call fgx_factors_test(gross, factors(:, 1:2), newf, fgx, n_folds=5)
  call check(fgx%status == 0, 'fgx status')
  call check(size(fgx%sdf_coefficients) == knew, 'fgx coefficients')
  call check(size(fgx%standard_errors) == knew, 'fgx standard errors')
  call check(all(fgx%standard_errors >= 0.0_dp), 'fgx nonnegative se')

  print '(a)', 'test_oracle_fgx: PASS'

contains
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a,1x,a)') 'FAIL:', trim(label)
      error stop 1
    end if
  end subroutine check
end program test_oracle_fgx
