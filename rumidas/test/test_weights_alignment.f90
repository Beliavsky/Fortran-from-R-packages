program test_weights_alignment
  use rumidas
  implicit none
  real(dp), allocatable :: w(:), wa(:), mat(:, :), component(:)
  real(dp) :: mv(6)
  integer :: periods(4), status

  w = beta_weights(4, 1.0_dp, 2.0_dp, status)
  call check(status == RUMIDAS_SUCCESS, 'beta status')
  call check_close(sum(w), 1.0_dp, 1.0e-13_dp, 'beta sum')
  call check_close(w(1), 0.5_dp, 1.0e-13_dp, 'beta first')
  call check_close(w(4), 0.0_dp, 1.0e-13_dp, 'beta endpoint')

  wa = exponential_almon_weights(3, 0.1_dp, -0.2_dp, status)
  call check(status == RUMIDAS_SUCCESS, 'almon status')
  call check_close(sum(wa), 1.0_dp, 1.0e-13_dp, 'almon sum')
  call check(all(wa > 0.0_dp), 'almon positivity')

  mv = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
  periods = [3, 4, 5, 6]
  call lag_matrix_from_period_index(mv, periods, 2, mat, status)
  call check(status == RUMIDAS_SUCCESS, 'lag matrix status')
  call check(all(abs(mat(:, 1) - [1.0_dp, 2.0_dp, 3.0_dp]) < 1.0e-14_dp), 'lag first column')
  call check(all(abs(mat(:, 4) - [4.0_dp, 5.0_dp, 6.0_dp]) < 1.0e-14_dp), 'lag last column')

  allocate(component(4))
  call midas_weighted_component(mat, 2, 2.0_dp, RUMIDAS_BETA_LAG, component, status)
  call check(status == RUMIDAS_SUCCESS, 'component status')
  call check(all(component > 0.0_dp), 'component positive')

  print '(a)', 'test_weights_alignment: PASS'
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) error stop message
  end subroutine check
  subroutine check_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    call check(abs(actual - expected) <= tolerance, message)
  end subroutine check_close
end program test_weights_alignment
