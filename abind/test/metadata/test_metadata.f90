program test_metadata
   use abind, only : dp, array_value, int_vector, init_array, set_dimnames, set_dimname_names, &
      abind_arrays, asub_array, acorn_array, has_dimnames
   implicit none

   call test_join_names()
   call test_subset_names()
   call test_acorn_labels()
   print '(a)', 'All abind metadata tests passed.'

contains

   subroutine test_join_names()
      type(array_value) :: inputs(2)
      type(array_value) :: out
      integer :: status

      call init_array(inputs(1), real([1, 2], dp), [2], status)
      call set_dimnames(inputs(1), 1, ['a', 'b'], status)
      call init_array(inputs(2), real([3, 4], dp), [2], status)
      call set_dimnames(inputs(2), 1, ['c', 'd'], status)
      call abind_arrays(inputs, out, status, arg_names=['x', 'y'], hier_names='before')
      call assert_equal_int(status, 0, 'hierarchical abind status')
      call assert_true(has_dimnames(out, 1), 'hierarchical abind has names')
      call assert_equal_strings(out%dimnames(1)%values, ['x.a', 'x.b', 'y.c', 'y.d'], 'hierarchical join names')
   end subroutine test_join_names

   subroutine test_subset_names()
      type(array_value) :: x
      type(array_value) :: out
      type(int_vector) :: idx(1)
      integer :: status

      call init_array(x, real([1, 2, 3, 4, 5, 6], dp), [2, 3], status)
      call set_dimnames(x, 1, ['r1', 'r2'], status)
      call set_dimnames(x, 2, ['A ', 'B ', 'C '], status)
      call set_dimname_names(x, ['rows', 'cols'], status)
      idx(1)%values = [3, 1]
      call asub_array(x, idx, [2], out, status, drop=.false.)
      call assert_equal_int(status, 0, 'named asub status')
      call assert_equal_strings(out%dimnames(2)%values, ['C ', 'A '], 'named asub labels')
   end subroutine test_subset_names

   subroutine test_acorn_labels()
      type(array_value) :: x
      type(array_value) :: out
      integer :: status
      integer :: i

      call init_array(x, [(real(i, dp), i = 1, 12)], [4, 3], status)
      call acorn_array(x, [-2, 2], out, status)
      call assert_equal_int(status, 0, 'acorn label status')
      call assert_equal_strings(out%dimnames(1)%values, ['[3]', '[4]'], 'acorn row labels')
      call assert_equal_strings(out%dimnames(2)%values, ['[1]', '[2]'], 'acorn column labels')
   end subroutine test_acorn_labels

   subroutine assert_equal_int(actual, expected, label)
      integer, intent(in) :: actual !! Observed integer value.
      integer, intent(in) :: expected !! Required integer value.
      character(len=*), intent(in) :: label !! Test label printed when the assertion fails.

      if (actual /= expected) then
         print '(a,2(1x,i0))', trim(label), actual, expected
         error stop 1
      end if
   end subroutine assert_equal_int

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Boolean assertion condition.
      character(len=*), intent(in) :: label !! Test label printed when the assertion fails.

      if (.not. condition) then
         print '(a)', trim(label)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_equal_strings(actual, expected, label)
      character(len=*), intent(in) :: actual(:) !! Observed character values.
      character(len=*), intent(in) :: expected(:) !! Required character values.
      character(len=*), intent(in) :: label !! Test label printed when the assertion fails.
      integer :: i

      if (size(actual) /= size(expected)) then
         print '(a)', trim(label)
         error stop 1
      end if
      do i = 1, size(actual)
         if (trim(actual(i)) /= trim(expected(i))) then
            print '(a,1x,i0,2(1x,a))', trim(label), i, trim(actual(i)), trim(expected(i))
            error stop 1
         end if
      end do
   end subroutine assert_equal_strings

end program test_metadata
