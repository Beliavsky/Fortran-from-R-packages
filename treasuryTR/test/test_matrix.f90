! SPDX-License-Identifier: MIT
! Copyright (c) 2021 treasuryTR authors

program test_matrix
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use treasurytr, only : dp, total_return, tt_success
  implicit none

  real(dp), parameter :: tol = 1.0e-14_dp
  real(dp) :: yields(5, 2), matrix_returns(5, 2), vector_returns(5)
  integer :: j, status
  character(len=160) :: message

  yields(:, 1) = [0.030_dp, 0.031_dp, 0.0305_dp, 0.032_dp, 0.0315_dp]
  yields(:, 2) = [0.045_dp, 0.044_dp, 0.0445_dp, 0.043_dp, 0.0425_dp]

  matrix_returns = total_return(yields, 7.0_dp, scale=52.0_dp, &
    source_compatible=.false., status=status, message=message)
  call assert_true(status == tt_success, 'matrix status')
  call assert_true(all(ieee_is_nan(matrix_returns(1, :))), 'matrix first row')

  do j = 1, 2
    vector_returns = total_return(yields(:, j), 7.0_dp, scale=52.0_dp, &
      source_compatible=.false.)
    call assert_close(maxval(abs(matrix_returns(2:, j) - vector_returns(2:))), &
      0.0_dp, tol, 'matrix column equivalence')
  end do

  print '(a)', 'test_matrix: PASS'

contains

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label

    if (abs(actual - expected) > tolerance) then
      print '(a,2(1x,es24.16))', trim(label), actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      print '(a)', trim(label)
      error stop 1
    end if
  end subroutine assert_true

end program test_matrix
