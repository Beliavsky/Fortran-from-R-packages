! SPDX-License-Identifier: MIT
! SPDX-FileComment: Deterministic tests for the Fortran tidyr translation.
program test_tidyr
   !! Checks pivots, missing values, character splitting, grids, and row expansion.
   use tidyr
   implicit none

   integer :: failures
   type(tibble_column) :: columns(3), replacements(1)
   type(tibble_type) :: data, long, wide, result

   failures = 0
   columns(1) = make_column('id', [1, 2])
   columns(2) = make_column('x', [10.0d0, 20.0d0])
   columns(3) = make_column('y', [11.0d0, 21.0d0])
   data = new_tibble(columns)

   long = pivot_longer(data, ['x', 'y'], 'variable', 'measurement')
   call check(long%nrow() == 4 .and. long%ncol() == 3, 'pivot_longer shape', failures)
   call check(all(long%columns(1)%integer_values == [1, 1, 2, 2]), 'pivot_longer identifiers', failures)
   call check(all(long%columns(3)%real_values == [10.0d0, 11.0d0, 20.0d0, 21.0d0]), &
              'pivot_longer values', failures)
   wide = pivot_wider(long, 'variable', 'measurement', ['id'])
   call check(wide%nrow() == 2 .and. wide%ncol() == 3, 'pivot_wider shape', failures)
   call check(all(wide%columns(2)%real_values == [10.0d0, 20.0d0]), 'pivot_wider x', failures)
   call check(all(wide%columns(3)%real_values == [11.0d0, 21.0d0]), 'pivot_wider y', failures)

   columns(1) = make_column('id', [1, 2, 3])
   columns(2) = make_column('value', [4, 0, 6], [.false., .true., .false.])
   columns(3) = make_column('label', ['a', 'b', 'c'])
   data = new_tibble(columns)
   replacements(1) = make_column('value', 5, 1)
   result = replace_na(data, replacements)
   call check(all(result%columns(2)%integer_values == [4, 5, 6]), 'replace_na', failures)
   result = fill(data, ['value'], 'down')
   call check(all(result%columns(2)%integer_values == [4, 4, 6]), 'fill down', failures)
   result = drop_na(data, ['value'])
   call check(result%nrow() == 2, 'drop_na', failures)

   columns(1) = make_column('id', [1, 2])
   columns(2) = make_column('code', ['a-b', 'c-d'])
   data = new_tibble(columns(:2))
   result = separate_wider_delim(data, 'code', ['left ', 'right'], '-')
   call check(result%ncol() == 3, 'separate_wider_delim shape', failures)
   call check(all(result%columns(2)%character_values == ['a', 'c']), 'separate left', failures)
   result = unite(result, 'code', ['left ', 'right'], sep='-')
   call check(all(result%columns(2)%character_values == ['a-b', 'c-d']), 'unite', failures)

   columns(1) = make_column('letter', ['b', 'a'])
   columns(2) = make_column('number', [2, 1])
   result = crossing(columns(:2))
   call check(result%nrow() == 4, 'crossing shape', failures)
   call check(all(full_seq([1, 3, 5], 2) == [1, 3, 5]), 'full_seq', failures)

   data = new_tibble([make_column('x', [10, 20])])
   result = uncount(data, [2, 1], id='copy')
   call check(all(result%columns(1)%integer_values == [10, 10, 20]), 'uncount values', failures)
   call check(all(result%columns(2)%integer_values == [1, 2, 1]), 'uncount identifiers', failures)

   if (failures > 0) error stop 'tidyr tests failed'
   print '(a)', 'All tidyr tests passed.'

contains

   subroutine check(condition, label, failure_count)
      !! Records a failed assertion.
      logical, intent(in) :: condition         !! Assertion result.
      character(len=*), intent(in) :: label    !! Assertion label.
      integer, intent(inout) :: failure_count  !! Accumulated failures.
      if (.not. condition) then
         print '(a)', 'FAIL: ' // label
         failure_count = failure_count + 1
      end if
   end subroutine check

end program test_tidyr
