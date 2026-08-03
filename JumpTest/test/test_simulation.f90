! SPDX-License-Identifier: MIT
program test_simulation
  use jumptest, only : dp, i8, JT_SUCCESS, simulation_result, lp_path, pvc_path, &
    pv2_path, sv, svj, sv1f, sv1fj, sv2f
  implicit none

  real(dp), allocatable :: price(:), variance1(:), variance2(:,:)
  real(dp) :: shocks2(3, 2), shocks3(3, 3)
  integer :: status
  type(simulation_result) :: first, second

  call lp_path(3.0_dp, 0.05_dp, 0.2_dp, 0.01_dp, 0.015_dp, 0.0001_dp, &
    [0.1_dp, -0.2_dp, 0.3_dp], [0.2_dp, -0.1_dp, 0.4_dp], &
    [0.0_dp, 0.01_dp, 0.0_dp], [3.8_dp, 4.0_dp, 2.0_dp], &
    price, variance1, status)
  call check(status == JT_SUCCESS, 'lp status')
  call check_close(price(4), 3.03083359_dp, 5.0e-9_dp, 'lp final price')
  call check_close(variance1(4), 0.20268583_dp, 5.0e-9_dp, 'lp final variance')

  shocks2 = reshape([0.2_dp, -0.1_dp, 0.5_dp, -0.3_dp, 0.4_dp, 0.2_dp], [3, 2])
  call pvc_path(3.0_dp, 0.001_dp, 0.0_dp, 0.125_dp, 0.5_dp, 0.1_dp, 0.99_dp, &
    shocks2, [0.0_dp, 0.01_dp, 0.0_dp], price, variance1, status)
  call check(status == JT_SUCCESS, 'pvc status')
  call check_close(price(4), 3.07691847_dp, 5.0e-9_dp, 'pvc final price')

  shocks3 = reshape([0.2_dp, -0.1_dp, 0.5_dp, -0.3_dp, 0.4_dp, 0.2_dp, &
    0.1_dp, -0.2_dp, 0.3_dp], [3, 3])
  call pv2_path(3.0_dp, 0.001_dp, -1.2_dp, 0.04_dp, 1.5_dp, 0.5_dp, 0.5_dp, &
    0.1_dp, 0.99_dp, 0.98_dp, 0.25_dp, shocks3, price, variance2, status)
  call check(status == JT_SUCCESS, 'pv2 status')
  call check_close(price(4), 2.96546301_dp, 5.0e-9_dp, 'pv2 final price')

  call sv1f(20, 3, first, seed=12345_i8)
  call sv1f(20, 3, second, seed=12345_i8)
  call check(first%status == JT_SUCCESS, 'SV1F status')
  call check(size(first%price) == 60, 'SV1F size')
  call check(maxval(abs(first%price - second%price)) <= tiny(1.0_dp), &
    'SV1F reproducibility')

  call svj(20, 2, first, lambda=0.0_dp, seed=2_i8)
  call check(first%status == JT_SUCCESS, 'SVJ status')
  call check(all(first%jump_count == 0), 'SVJ no jumps')
  call check(size(first%price) == 40, 'SVJ size')

  call sv(20, 2, first, seed=3_i8)
  call check(first%status == JT_SUCCESS, 'SV status')
  call check(size(first%price) == 41, 'SV preserves upstream initial price')
  call check(all(first%variance >= 0.0_dp), 'SV nonnegative variance')

  call sv1fj(20, 2, first, lambda=0.0_dp, seed=4_i8)
  call check(all(first%jump_count == 0), 'SV1FJ no jumps')

  call sv2f(20, 2, first, seed=5_i8)
  call check(first%status == JT_SUCCESS, 'SV2F status')
  call check(all(shape(first%variance) == [40, 2]), 'SV2F variance shape')

  print '(a)', 'test_simulation: PASS'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      print '(a,1x,a)', 'FAIL:', message
      error stop 1
    end if
  end subroutine check

  subroutine check_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    call check(abs(actual - expected) <= tolerance, message)
  end subroutine check_close

end program test_simulation
