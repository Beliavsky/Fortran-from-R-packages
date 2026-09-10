! SPDX-License-Identifier: MIT
! SPDX-FileComment: Tokenizer tests for the Fortran readr translation.
program test_tokenizer
   !! Tests quoting, escaped quotes, quoted newlines, comments, and ragged records.
   use readr
   implicit none

   type(token_table_type) :: tokens
   character(len=:), allocatable :: text

   text = "a,b,c" // achar(10) // "1,""two,too"",3" // achar(10) // &
      "2,""line" // achar(10) // "break"",""a""""b"""
   tokens = tokenize(text)
   call assert_true(tokens%valid, "balanced input")
   call assert_true(tokens%nrow() == 3 .and. tokens%ncol() == 3, "dimensions")
   call assert_equal(tokens%field(2, 2), "two,too", "quoted delimiter")
   call assert_equal(tokens%field(3, 2), "line" // achar(10) // "break", "quoted newline")
   call assert_equal(tokens%field(3, 3), 'a"b', "doubled quote")
   call assert_true(all(tokens%fields_per_row == 3), "field counts")

   tokens = tokenize("a,b" // achar(13) // achar(10) // "1,2" // achar(13) // achar(10))
   call assert_true(tokens%nrow() == 2, "CRLF records")
   tokens = tokenize("# ignored" // achar(10) // "a,b" // achar(10) // "1,2", comment="#")
   call assert_true(tokens%nrow() == 2, "comment row")
   tokens = tokenize("a,b" // achar(10) // "1")
   call assert_true(tokens%fields_per_row(2) == 1, "ragged count")
   tokens = tokenize("a,""unterminated")
   call assert_true(.not. tokens%valid, "unterminated quote")

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

end program test_tokenizer
