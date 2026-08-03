! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
program test_household
  use r4good_personal_finances
  implicit none
  type(household) :: home
  type(household_member) :: older, younger
  type(household_timeline) :: timeline
  type(date_type) :: current_date
  integer :: st

  current_date = date_from_string("2020-07-15", st)
  call assert_true(st == r4gpf_success, "date parsing")
  older%name = "older"
  older%birth_date = date_from_string("1980-02-15")
  older%mode = 80.0_dp
  older%dispersion = 10.0_dp
  call older%add_event("retirement", 45.0_dp, status=st)
  younger%name = "younger"
  younger%birth_date = date_from_string("1990-07-15")
  younger%mode = 85.0_dp
  younger%dispersion = 9.0_dp
  call younger%add_event("kid", 35.0_dp, years=2.0_dp, status=st)
  call home%add_member(older, st)
  call home%add_member(younger, st)
  call build_household_timeline(home, current_date, timeline, st)
  call assert_true(st == r4gpf_success, "timeline status")
  call assert_true(timeline%n_members == 2, "timeline member count")
  call assert_true(timeline%n_periods == ceiling(home%lifespan(current_date)) + 1, "timeline period count")
  call assert_true(timeline%index(1) == 0 .and. timeline%years_left(timeline%n_periods) == 0, "timeline indices")
  call assert_close(timeline%joint_survival(1), 1.0_dp, 1.0e-12_dp, "joint survival start")
  call assert_true(member_event_is_on(older, "retirement", 45.0_dp), "open ended event")
  call assert_true(member_event_is_on(younger, "kid", 36.0_dp), "finite event")
  call assert_true(.not. member_event_is_on(younger, "kid", 37.0_dp), "event end")
  print *, "test_household: PASS"
contains
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      print *, "FAIL: ", trim(message)
      error stop 1
    end if
  end subroutine assert_true
  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    if (abs(actual - expected) > tolerance) then
      print *, "FAIL: ", trim(message), actual, expected
      error stop 1
    end if
  end subroutine assert_close
end program test_household
