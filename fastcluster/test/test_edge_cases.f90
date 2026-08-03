! SPDX-License-Identifier: BSD-2-Clause
program test_edge_cases
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_positive_inf, ieee_quiet_nan
  use fastcluster
  implicit none

  real(dp) :: bad(3, 3), d(3, 3), members(3), x(3, 1)
  real(dp), allocatable :: condensed(:)
  type(hclust_result) :: result

  d = 0.0_dp
  d(2, 1) = 1.0_dp
  d(1, 2) = 1.0_dp
  d(3, 1) = 4.0_dp
  d(1, 3) = 4.0_dp
  d(3, 2) = 5.0_dp
  d(2, 3) = 5.0_dp
  members = [2.0_dp, 1.0_dp, 1.0_dp]
  call hclust_matrix(d, 'average', result, members)
  call check(result%ok(), 'weighted members status')
  call check_close(result%height(2), 13.0_dp / 3.0_dp, 1.0e-13_dp, &
    'weighted members update')

  call matrix_to_condensed(d, condensed)
  call hclust(condensed, 3, 'complete', result)
  call check(result%ok(), 'condensed clustering')
  call check_close(result%height(2), 5.0_dp, 1.0e-13_dp, 'condensed height')

  call hclust_matrix(d, 'not-a-method', result)
  call check(result%status == fc_invalid_argument, 'invalid method status')

  bad = d
  bad(1, 2) = 2.0_dp
  call hclust_matrix(bad, 'single', result)
  call check(result%status == fc_invalid_argument, 'asymmetric distance status')

  bad = d
  bad(2, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
  bad(1, 2) = bad(2, 1)
  call hclust_matrix(bad, 'single', result)
  call check(result%status == fc_nan_distance, 'NaN distance status')

  x(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp]
  call hclust_vector(x, 'ward', result, metric='manhattan')
  call check(result%status == fc_invalid_argument, 'invalid vector metric status')

  call hclust_matrix(d, 'average', result, [1.0_dp, 0.0_dp, 1.0_dp])
  call check(result%status == fc_invalid_argument, 'invalid members status')

  bad = ieee_value(0.0_dp, ieee_positive_inf)
  bad(1, 1) = 0.0_dp
  bad(2, 2) = 0.0_dp
  bad(3, 3) = 0.0_dp
  call hclust_matrix(bad, 'single', result)
  call check(result%ok(), 'infinite single-linkage status')
  call check(all(.not. ieee_is_finite(result%height)), 'infinite single-linkage heights')

  call hclust_matrix(bad, 'ward.D2', result)
  call check(result%ok(), 'infinite Ward status')
  call check(all(.not. ieee_is_finite(result%height)), 'infinite Ward heights')

  print '(a)', 'test_edge_cases: PASS'

contains

  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write (*, '(a)') 'FAIL: '//label
      error stop 1
    end if
  end subroutine check

  subroutine check_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label

    call check(abs(actual - expected) <= tolerance * max(1.0_dp, abs(expected)), label)
  end subroutine check_close

end program test_edge_cases
