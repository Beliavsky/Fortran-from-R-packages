! SPDX-License-Identifier: BSD-2-Clause
program test_distances
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use fastcluster
  implicit none

  real(dp), allocatable :: d(:, :), condensed(:), rebuilt(:, :)
  real(dp) :: x(2, 3), xb(2, 4), xm(2, 2)
  integer :: status
  character(len=:), allocatable :: message

  x(1, :) = [1.0_dp, -2.0_dp, 0.0_dp]
  x(2, :) = [4.0_dp, 2.0_dp, 0.0_dp]

  call pairwise_distances(x, 'euclidean', d, status=status, message=message)
  call check(status == fc_success, 'Euclidean status')
  call check_close(d(1, 2), 5.0_dp, 1.0e-13_dp, 'Euclidean distance')

  call pairwise_distances(x, 'maximum', d, status=status)
  call check_close(d(1, 2), 4.0_dp, 1.0e-13_dp, 'maximum distance')

  call pairwise_distances(x, 'manhattan', d, status=status)
  call check_close(d(1, 2), 7.0_dp, 1.0e-13_dp, 'Manhattan distance')

  call pairwise_distances(x, 'canberra', d, status=status)
  call check_close(d(1, 2), 2.4_dp, 1.0e-13_dp, 'Canberra distance')

  call pairwise_distances(x, 'minkowski', d, p=3.0_dp, status=status)
  call check_close(d(1, 2), 91.0_dp ** (1.0_dp / 3.0_dp), 1.0e-13_dp, &
    'Minkowski distance')

  xb(1, :) = [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
  xb(2, :) = [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
  call pairwise_distances(xb, 'binary', d, status=status)
  call check_close(d(1, 2), 2.0_dp / 3.0_dp, 1.0e-13_dp, 'binary distance')

  xm(1, :) = [0.0_dp, ieee_value(0.0_dp, ieee_quiet_nan)]
  xm(2, :) = [3.0_dp, 4.0_dp]
  call pairwise_distances(xm, 'euclidean', d, status=status)
  call check_close(d(1, 2), sqrt(18.0_dp), 1.0e-13_dp, 'missing-value scaling')

  call matrix_to_condensed(d, condensed, status=status)
  call check(status == fc_success, 'matrix_to_condensed status')
  call condensed_to_matrix(condensed, 2, rebuilt, status=status)
  call check(status == fc_success, 'condensed_to_matrix status')
  call check_close(rebuilt(1, 2), d(1, 2), 1.0e-13_dp, 'condensed round trip')

  print '(a)', 'test_distances: PASS'

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

end program test_distances
