! SPDX-License-Identifier: MIT
! SPDX-FileComment: Column parser tests for the Fortran readr translation.
program test_parsers
   !! Tests scalar types, inference, diagnostics, locales, and factors.
   use readr
   implicit none

   type(parse_result_type) :: parsed
   type(factor_parse_result_type) :: categorical

   parsed = parse_logical([character(len=5) :: "TRUE", "F", "NA", "bad"])
   call assert_true(.not. parsed%ok(), "logical problem")
   call assert_true(size(parsed%problem) == 1 .and. parsed%problem(1)%row == 4, "logical diagnostic")
   call assert_true(parsed%values%logical_values(1), "logical true")
   call assert_true(.not. parsed%values%logical_values(2), "logical false")
   call assert_true(parsed%values%missing(3) .and. parsed%values%missing(4), "logical missing")

   parsed = parse_integer([character(len=4) :: "-2", "+3", "4.5", "NA"])
   call assert_true(parsed%values%integer_values(1) == -2, "negative integer")
   call assert_true(parsed%values%integer_values(2) == 3, "positive integer")
   call assert_true(size(parsed%problem) == 1, "integer diagnostic")

   parsed = parse_double([character(len=5) :: "1,25", "-2,5"], settings=locale(decimal_mark=",", grouping_mark="."))
   call assert_true(maxval(abs(parsed%values%real_values - [1.25_dp, -2.5_dp])) < 1.0e-14_dp, "decimal comma")
   parsed = parse_number([character(len=12) :: "$1,234.50", "EUR -20.5"])
   call assert_true(maxval(abs(parsed%values%real_values - [1234.5_dp, -20.5_dp])) < 1.0e-12_dp, "flexible number")

   call assert_true(guess_parser([character(len=5) :: "TRUE", "FALSE"]) == "logical", "guess logical")
   call assert_true(guess_parser([character(len=2) :: "1", "20"]) == "integer", "guess integer")
   call assert_true(guess_parser([character(len=3) :: "1.5", "2"]) == "double", "guess double")
   call assert_true(guess_parser([character(len=3) :: "one", "2"]) == "character", "guess character")

   categorical = parse_factor([character(len=4) :: "low", "high", "bad"], &
      [character(len=4) :: "low", "high"], ordered=.true.)
   call assert_true(categorical%values%ordered, "ordered factor")
   call assert_true(categorical%values%missing(3), "unknown factor missing")
   call assert_true(size(categorical%problem) == 1, "factor diagnostic")

contains

   pure elemental subroutine assert_true(condition, message)
      !! Stops the test if a condition is false.
      logical, intent(in) :: condition        !! Condition that must hold.
      character(len=*), intent(in) :: message !! Failure description.

      if (.not. condition) error stop "FAIL: " // message
   end subroutine assert_true

end program test_parsers
