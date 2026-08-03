! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
program test_portfolio
  use r4good_personal_finances
  implicit none
  type(portfolio_result) :: result
  real(dp) :: mu3(3), sd3(3), corr3(3,3)
  real(dp) :: mu9(9), sd9(9), corr9(9,9)
  integer :: i

  mu3 = [0.0472_dp, 0.0504_dp, 0.0275_dp]
  sd3 = [0.1588_dp, 0.1718_dp, 0.0562_dp]
  corr3 = reshape([1.00_dp,0.87_dp,0.21_dp, 0.87_dp,1.00_dp,0.37_dp, 0.21_dp,0.37_dp,1.00_dp],[3,3])
  call optimize_portfolio(0.35_dp, mu3, sd3, corr3, result)
  call assert_true(result%status == r4gpf_success, "3-asset optimization status")
  call assert_vector_close(result%total, [0.2729395_dp, 0.0887746_dp, 0.6382859_dp], 2.0e-5_dp, "3-asset allocation")

  mu9 = [0.0468_dp,0.0501_dp,0.0505_dp,0.0540_dp,0.0269_dp,0.0288_dp,0.0190_dp,0.0329_dp,0.0250_dp]
  sd9 = [0.1542_dp,0.1795_dp,0.1671_dp,0.2142_dp,0.0379_dp,0.0581_dp,0.03138274_dp,0.0833_dp,0.0055_dp]
  corr9 = 0.0_dp
  do i = 1, 9
    corr9(i,i) = 1.0_dp
  end do
  call optimize_portfolio(0.35_dp, mu9, sd9, corr9, result)
  call assert_true(result%status == r4gpf_success, "9-asset optimization status")
  call assert_vector_close(result%total, &
       [0.2465167_dp, 0.2207499_dp, 0.2601587_dp, 0.1872452_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 0.0853295_dp, 0.0_dp], &
    3.0e-5_dp, "9-asset allocation")
  call assert_close(sum(result%total), 1.0_dp, 1.0e-12_dp, "allocation sum")
  print *, "test_portfolio: PASS"
contains
  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    if (abs(actual - expected) > tolerance) then
      print *, "FAIL: ", trim(message), actual, expected
      error stop 1
    end if
  end subroutine assert_close
  subroutine assert_vector_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual(:), expected(:), tolerance
    character(len=*), intent(in) :: message
    if (size(actual) /= size(expected) .or. maxval(abs(actual - expected)) > tolerance) then
      print *, "FAIL: ", trim(message)
      print *, actual
      print *, expected
      error stop 1
    end if
  end subroutine assert_vector_close
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      print *, "FAIL: ", trim(message)
      error stop 1
    end if
  end subroutine assert_true
end program test_portfolio
