! SPDX-License-Identifier: GPL-3.0-only
program test_daycount
  use bondvaluation, only: dp, date_from_ymd, accrued_interest, accrued_interest_result
  implicit none
  integer, parameter :: expected_days(16) = [182, 182, 182, 182, 179, 179, &
    180, 180, 179, 182, 182, 182, 179, 182, 182, 124]
  real(dp), parameter :: expected_fraction(16) = [ &
    0.49818848716221276_dp, 0.5_dp, 0.4986301369863014_dp, &
    0.4972677595628415_dp, 0.49722222222222223_dp, &
    0.49722222222222223_dp, 0.5_dp, 0.5_dp, &
    0.49722222222222223_dp, 0.4986301369863014_dp, &
    0.4986301369863014_dp, 0.5055555555555555_dp, &
    0.4904109589041096_dp, 0.5_dp, 0.5_dp, &
    0.49206349206349204_dp]
  real(dp), parameter :: expected_accrued(16) = [ &
    261.5489557601617_dp, 262.5_dp, 261.7808219178082_dp, &
    261.0655737704918_dp, 261.0416666666667_dp, &
    261.0416666666667_dp, 262.5_dp, 262.5_dp, &
    261.0416666666667_dp, 261.7808219178082_dp, &
    261.7808219178082_dp, 265.41666666666663_dp, &
    257.4657534246575_dp, 262.5_dp, 262.5_dp, &
    254.9768980174072_dp]
  type(accrued_interest_result) :: result
  integer :: i

  do i = 1, 16
    result = accrued_interest(date_from_ymd(2011, 8, 31), &
      date_from_ymd(2012, 2, 29), 5.25_dp, i, 10000.0_dp, 2, &
      date_from_ymd(2021, 8, 31), 2012, .true.)
    call assert_equal_int(result%status, 0, "day-count status")
    call assert_equal_int(result%days_accrued, expected_days(i), "day-count days")
    call assert_close(result%year_fraction_value, expected_fraction(i), 2.0e-14_dp, &
      "day-count fraction")
    call assert_close(result%accrued_interest, expected_accrued(i), 2.0e-10_dp, &
      "accrued interest")
  end do
  print '(a)', "test_daycount: PASS"

contains

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tolerance) then
      write(*, '(a,2(1x,es24.16))') trim(label)//" mismatch:", actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_equal_int(actual, expected, label)
    integer, intent(in) :: actual, expected
    character(len=*), intent(in) :: label
    if (actual /= expected) then
      write(*, '(a,2(1x,i0))') trim(label)//" mismatch:", actual, expected
      error stop 1
    end if
  end subroutine assert_equal_int

end program test_daycount
