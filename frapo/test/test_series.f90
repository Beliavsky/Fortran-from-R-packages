! SPDX-License-Identifier: GPL-3.0-or-later
program test_series
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use frapo
  implicit none

  real(dp) :: prices(3), y(4), matrix(3, 2), a(2, 2)
  real(dp), allocatable :: values(:), converted(:), trend(:), root(:, :), matrix_trend(:, :)
  integer :: status

  prices = [100.0_dp, 110.0_dp, 99.0_dp]
  values = returnseries(prices, method=returns_discrete, percentage=.false.)
  call assert_true(ieee_is_nan(values(1)), 'untrimmed return starts with NaN')
  call assert_close(values(2), 0.1_dp, 1.0e-13_dp, 'first discrete return')
  call assert_close(values(3), -0.1_dp, 1.0e-13_dp, 'second discrete return')

  values = returnseries(prices, method=returns_discrete, percentage=.false., compound=.true.)
  call assert_vector(values, [0.0_dp, 0.1_dp, -0.01_dp], 1.0e-13_dp, 'compound discrete returns')

  converted = returnconvert([1.0_dp, -2.0_dp], direction=convert_cont_to_disc, percentage=.true.)
  values = returnconvert(converted, direction=convert_disc_to_cont, percentage=.true.)
  call assert_vector(values, [1.0_dp, -2.0_dp], 1.0e-12_dp, 'return conversion inverse')

  values = capser([-2.0_dp, 0.5_dp, 4.0_dp], -1.0_dp, 2.0_dp)
  call assert_vector(values, [-1.0_dp, 0.5_dp, 2.0_dp], 0.0_dp, 'cap series')

  y = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
  trend = trdsma(y, 2, trim=.true., status=status)
  call assert_equal_int(status, frapo_ok, 'SMA status')
  call assert_vector(trend, [1.5_dp, 2.5_dp, 3.5_dp], 1.0e-14_dp, 'SMA')

  trend = trdwma(y, [0.75_dp, 0.25_dp], trim=.true.)
  call assert_vector(trend, [1.75_dp, 2.75_dp, 3.75_dp], 1.0e-14_dp, 'WMA orientation')

  trend = trdes(y, 0.5_dp)
  call assert_vector(trend, [0.5_dp, 1.25_dp, 2.125_dp, 3.0625_dp], 1.0e-14_dp, 'exponential trend')

  trend = trdhp(y, 1600.0_dp, status)
  call assert_equal_int(status, frapo_ok, 'HP status')
  call assert_vector(trend, y, 1.0e-10_dp, 'HP preserves linear series')

  trend = trdbinary([-1.0_dp, 0.0_dp, 1.0_dp])
  call assert_close(trend(1), -1.0_dp, 1.0e-14_dp, 'binary negative')
  call assert_close(trend(2), 0.0_dp, 0.0_dp, 'binary zero')
  call assert_close(trend(3), 1.0_dp, 1.0e-14_dp, 'binary positive')

  matrix = reshape([1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp, 4.0_dp, 6.0_dp], [3, 2])
  matrix_trend = trdsma(matrix, 2, trim=.true.)
  call assert_true(all(shape(matrix_trend) == [2, 2]), 'matrix trend shape')
  call assert_vector(matrix_trend(:, 2), [3.0_dp, 5.0_dp], 1.0e-14_dp, 'matrix SMA')

  a = 0.0_dp
  a(1, 1) = 4.0_dp
  a(2, 2) = 9.0_dp
  root = sqrm(a, status)
  call assert_equal_int(status, frapo_ok, 'matrix square root status')
  call assert_close(root(1, 1), 2.0_dp, 1.0e-13_dp, 'matrix square root first')
  call assert_close(root(2, 2), 3.0_dp, 1.0e-13_dp, 'matrix square root second')

  print '(a)', 'test_series: PASS'

contains

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tolerance) then
      write(*, '(a,2es24.15)') trim(label)//': ', actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_vector(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual(:), expected(:), tolerance
    character(len=*), intent(in) :: label
    if (size(actual) /= size(expected)) then
      write(*, '(a)') trim(label)//': size mismatch'
      error stop 1
    end if
    if (maxval(abs(actual - expected)) > tolerance) then
      write(*, '(a)') trim(label)//': values differ'
      write(*, '(*(es20.11,1x))') actual
      write(*, '(*(es20.11,1x))') expected
      error stop 1
    end if
  end subroutine assert_vector

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a)') trim(label)//': failed'
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_equal_int(actual, expected, label)
    integer, intent(in) :: actual, expected
    character(len=*), intent(in) :: label
    if (actual /= expected) then
      write(*, '(a,2i8)') trim(label)//': ', actual, expected
      error stop 1
    end if
  end subroutine assert_equal_int
end program test_series
