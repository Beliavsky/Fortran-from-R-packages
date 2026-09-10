! SPDX-License-Identifier: MIT
! SPDX-FileComment: Deterministic tests for the Fortran stringr translation.
program test_stringr
   !! Checks fixed matching, transformation, splitting, and vector operations.
   use stringr
   implicit none
   integer :: failures
   integer, allocatable :: indices(:)
   character(len=:), allocatable :: text
   character(len=8), allocatable :: fields(:)
   type(location_type) :: location
   type(split_result_type) :: split

   failures = 0
   call check(str_length('abc  ') == 3, 'length', failures)
   call check(str_detect('alphabet', 'pha'), 'detect', failures)
   call check(str_starts('alphabet', 'alpha') .and. str_ends('alphabet', 'bet'), 'bounds', failures)
   call check(str_count('banana', 'an') == 2, 'count', failures)
   location = str_locate('alphabet', 'pha')
   call check(location%start == 3 .and. location%end == 5, 'locate', failures)
   call check(trim(str_sub('abcdef', 2, 4)) == 'bcd', 'substring', failures)
   call check(trim(str_sub('abcdef', -3)) == 'def', 'negative substring', failures)
   call check(trim(str_to_lower('AbC')) == 'abc', 'lowercase', failures)
   call check(trim(str_to_upper('AbC')) == 'ABC', 'uppercase', failures)
   call check(trim(str_to_title('HELLO world')) == 'Hello World', 'title case', failures)
   call check(trim(str_squish('  a   b  ')) == 'a b', 'squish', failures)
   call check(str_dup('ab', 3) == 'ababab', 'duplicate', failures)
   call check(str_pad('a', 3, 'left', '0') == '00a', 'padding', failures)
   call check(str_trunc('abcdef', 5) == 'ab...', 'truncation', failures)
   text = str_replace('banana', 'an', 'XX')
   call check(text == 'bXXana', 'replace first', failures)
   text = str_replace_all('banana', 'an', 'X')
   call check(text == 'bXXa', 'replace all', failures)
   call check(str_remove('abc', 'a') == 'bc', 'remove at start', failures)
   call check(str_remove('abc', 'c') == 'ab', 'remove at end', failures)
   call check(str_like('alphabet', 'a%bet'), 'wildcard match', failures)
   call check(trim(str_to_snake('Hello world')) == 'hello_world', 'snake case', failures)
   call check(trim(str_to_camel('hello-world')) == 'helloWorld', 'camel case', failures)
   call check(trim(word('one two three', 2)) == 'two', 'word', failures)

   call check(str_flatten(['a', 'b', 'c'], ',') == 'a,b,c', 'flatten', failures)
   indices = str_order(['c', 'a', 'b'])
   call check(all(indices == [2, 3, 1]), 'order', failures)
   call check(all(str_sort(['c', 'a', 'b']) == ['a', 'b', 'c']), 'sort', failures)
   call check(all(str_unique(['a', 'b', 'a']) == ['a', 'b']), 'unique', failures)
   indices = str_which(['alpha', 'beta ', 'gamma'], 'a')
   call check(all(indices == [1, 2, 3]), 'which', failures)
   fields = str_split_1('a--b--c', '--')
   call check(size(fields) == 3 .and. all(fields == ['a       ', 'b       ', 'c       ']), 'split one', failures)
   split = str_split(['a,b  ', 'c,d,e'], ',')
   call check(all(split%counts == [2, 3]), 'split counts', failures)
   call check(trim(split%values(2, 3)) == 'e', 'split value', failures)

   if (failures > 0) error stop 'stringr tests failed'
   print '(a)', 'All stringr tests passed.'

contains

   subroutine check(condition, label, failure_count)
      !! Records a failed assertion.
      logical, intent(in) :: condition          !! Assertion result.
      character(len=*), intent(in) :: label     !! Assertion label.
      integer, intent(inout) :: failure_count   !! Accumulated failures.
      if (.not. condition) then
         print '(a)', 'FAIL: ' // label
         failure_count = failure_count + 1
      end if
   end subroutine check
end program test_stringr
