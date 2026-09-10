! SPDX-License-Identifier: MIT
! SPDX-FileComment: Rectangular I/O tests for the Fortran readr translation.
program test_read_write
   !! Tests typed CSV reading, explicit specifications, diagnostics, and formatting.
   use readr
   use tibble, only: tibble_character, tibble_integer, tibble_logical, tibble_real
   implicit none

   type(read_result_type) :: result, forced, roundtrip, variant
   type(col_spec_type) :: specification
   type(collector_type) :: score_collector(1)
   type(formatted_lines_type) :: formatted

   result = read_csv("data/sample.csv")
   call assert_true(result%ok(), "sample CSV diagnostics")
   call assert_true(result%data%nrow() == 3 .and. result%data%ncol() == 4, "sample dimensions")
   call assert_true(result%data%columns(1)%type_code == tibble_integer, "integer inference")
   call assert_true(result%data%columns(2)%type_code == tibble_logical, "logical inference")
   call assert_true(result%data%columns(3)%type_code == tibble_real, "real inference")
   call assert_true(result%data%columns(4)%type_code == tibble_character, "character inference")
   call assert_true(result%data%columns(3)%missing(2), "missing score")
   call assert_equal(result%data%columns(4)%character_values(2), "Smith, Bob", "quoted comma")
   call assert_equal(result%data%columns(4)%character_values(3), 'quote "inside"', "escaped quote")

   score_collector(1) = col_integer()
   specification = cols([character(len=5) :: "score"], score_collector)
   forced = read_csv("data/sample.csv", col_types=specification)
   call assert_true(.not. forced%ok(), "forced integer problem")
   call assert_true(size(forced%problem) == 2, "forced integer problem count")
   call assert_true(forced%problem(1)%column == 3 .and. forced%problem(1)%row == 2, "problem coordinates")

   variant = read_csv("data/ragged.csv")
   call assert_true(.not. variant%ok(), "ragged-row diagnostic")
   call assert_true(variant%problem(1)%row == 3, "ragged-row coordinate")
   variant = read_csv2("data/sample_csv2.csv")
   call assert_true(variant%ok(), "CSV2 diagnostics")
   call assert_true(maxval(abs(variant%data%columns(2)%real_values - [1.25_dp, -2.5_dp])) < &
      1.0e-14_dp, "CSV2 decimal mark")
   variant = read_table("data/sample_table.txt")
   call assert_true(variant%ok(), "whitespace table diagnostics")
   call assert_true(variant%data%columns(2)%type_code == tibble_logical, "whitespace table inference")

   formatted = format_csv(result%data)
   call assert_equal(formatted%line(1), "id,active,score,name", "formatted header")
   call assert_equal(formatted%line(3), '2,FALSE,NA,"Smith, Bob"', "formatted quoted row")
   call write_csv(result%data, "build/readr_roundtrip.csv")
   roundtrip = read_csv("build/readr_roundtrip.csv")
   call assert_true(roundtrip%ok(), "roundtrip diagnostics")
   call assert_true(all(roundtrip%data%columns(1)%integer_values == [1, 2, 3]), "roundtrip integer")
   call assert_equal(roundtrip%data%columns(4)%character_values(3), 'quote "inside"', "roundtrip quote")

contains

   pure elemental subroutine assert_true(condition, message)
      !! Stops the test if a condition is false.
      logical, intent(in) :: condition        !! Condition that must hold.
      character(len=*), intent(in) :: message !! Failure description.

      if (.not. condition) error stop "FAIL: " // message
   end subroutine assert_true

   pure elemental subroutine assert_equal(actual, expected, message)
      !! Stops the test if trimmed character values differ.
      character(len=*), intent(in) :: actual   !! Computed value.
      character(len=*), intent(in) :: expected !! Required value.
      character(len=*), intent(in) :: message  !! Failure description.

      call assert_true(trim(actual) == expected, message)
   end subroutine assert_equal

end program test_read_write
