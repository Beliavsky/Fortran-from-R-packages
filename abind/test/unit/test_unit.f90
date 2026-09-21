program test_unit
   use abind, only : dp, array_value, int_vector, init_array, abind_arrays, asub_array, &
      afill_array, adrop_array, acorn_array, set_dimnames
   implicit none

   call test_abind_vectors()
   call test_abind_inserted_dimensions()
   call test_abind_mixed_rank()
   call test_abind_zero_length()
   call test_asub()
   call test_adrop()
   call test_acorn()
   call test_afill()
   print '(a)', 'All abind unit tests passed.'

contains

   subroutine test_abind_vectors()
      type(array_value) :: inputs(2)
      type(array_value) :: out
      integer :: status

      call init_array(inputs(1), real([1, 2, 3, 4], dp), [4], status)
      call assert_equal_int(status, 0, 'init vector 1')
      call init_array(inputs(2), real([5, 6, 7, 8], dp), [4], status)
      call assert_equal_int(status, 0, 'init vector 2')

      call abind_arrays(inputs, out, status)
      call assert_equal_int(status, 0, 'abind vector default status')
      call assert_equal_int_array(out%dims, [8], 'abind vector default dims')
      call assert_close_array(out%data, real([1, 2, 3, 4, 5, 6, 7, 8], dp), 'abind vector default data')

      call abind_arrays(inputs, out, status, along=2.0_dp)
      call assert_equal_int(status, 0, 'abind vectors along 2 status')
      call assert_equal_int_array(out%dims, [4, 2], 'abind vectors along 2 dims')
      call assert_close_array(out%data, real([1, 2, 3, 4, 5, 6, 7, 8], dp), 'abind vectors along 2 data')
   end subroutine test_abind_vectors

   subroutine test_abind_inserted_dimensions()
      type(array_value) :: inputs(2)
      type(array_value) :: out
      integer :: status

      call init_array(inputs(1), [(real(status, dp), status = 1, 12)], [3, 4], status)
      call assert_equal_int(status, 0, 'init matrix 1')
      call init_array(inputs(2), [(real(status + 100, dp), status = 1, 12)], [3, 4], status)
      call assert_equal_int(status, 0, 'init matrix 2')

      call abind_arrays(inputs, out, status, along=1.5_dp)
      call assert_equal_int(status, 0, 'abind matrices along 1.5 status')
      call assert_equal_int_array(out%dims, [3, 2, 4], 'abind matrices along 1.5 dims')
      call assert_close_array(out%data(1:12), real([1, 2, 3, 101, 102, 103, 4, 5, 6, 104, 105, 106], dp), &
         'abind matrices along 1.5 first half')

      call abind_arrays(inputs, out, status, along=0.5_dp)
      call assert_equal_int(status, 0, 'abind matrices along 0.5 status')
      call assert_equal_int_array(out%dims, [2, 3, 4], 'abind matrices along 0.5 dims')
      call assert_close_array(out%data(1:8), real([1, 101, 2, 102, 3, 103, 4, 104], dp), &
         'abind matrices along 0.5 data')

      call abind_arrays(inputs, out, status, rev_along=0.0_dp)
      call assert_equal_int(status, 0, 'abind matrices reverse along zero status')
      call assert_equal_int_array(out%dims, [3, 4, 2], 'abind matrices reverse along zero dims')
   end subroutine test_abind_inserted_dimensions


   subroutine test_abind_mixed_rank()
      type(array_value) :: inputs(2)
      type(array_value) :: out
      integer :: status
      integer :: i

      call init_array(inputs(1), real([1, 2, 3, 4], dp), [4], status)
      call assert_equal_int(status, 0, 'init mixed-rank vector')
      call init_array(inputs(2), [(real(i + 4, dp), i = 1, 16)], [4, 4], status)
      call assert_equal_int(status, 0, 'init mixed-rank matrix')

      call abind_arrays(inputs, out, status, along=2.0_dp)
      call assert_equal_int(status, 0, 'mixed-rank cbind status')
      call assert_equal_int_array(out%dims, [4, 5], 'mixed-rank cbind dims')
      call assert_close_array(out%data(1:8), real([1, 2, 3, 4, 5, 6, 7, 8], dp), &
         'mixed-rank cbind data')

      call abind_arrays(inputs, out, status, along=1.0_dp)
      call assert_equal_int(status, 0, 'mixed-rank rbind status')
      call assert_equal_int_array(out%dims, [5, 4], 'mixed-rank rbind dims')
      call assert_close_array(out%data(1:10), real([1, 5, 6, 7, 8, 2, 9, 10, 11, 12], dp), &
         'mixed-rank rbind data')
   end subroutine test_abind_mixed_rank

   subroutine test_abind_zero_length()
      type(array_value) :: inputs(2)
      type(array_value) :: out
      real(dp), allocatable :: empty(:)
      integer :: status

      allocate(empty(0))
      call init_array(inputs(1), empty, [0], status)
      call assert_equal_int(status, 0, 'init first zero-length vector')
      call init_array(inputs(2), empty, [0], status)
      call assert_equal_int(status, 0, 'init second zero-length vector')
      call abind_arrays(inputs, out, status)
      call assert_equal_int(status, 0, 'zero-length abind status')
      call assert_equal_int_array(out%dims, [0], 'zero-length abind dims')
      call assert_equal_int(size(out%data), 0, 'zero-length abind data size')
   end subroutine test_abind_zero_length

   subroutine test_asub()
      type(array_value) :: x
      type(array_value) :: out
      type(int_vector) :: idx(2)
      integer :: status
      integer :: i

      call init_array(x, [(real(i, dp), i = 1, 24)], [2, 3, 4], status)
      call assert_equal_int(status, 0, 'init asub array')
      idx(1)%values = [1, 2]
      idx(2)%values = [2]
      call asub_array(x, idx, [1, 3], out, status)
      call assert_equal_int(status, 0, 'asub status')
      call assert_equal_int_array(out%dims, [2, 3], 'asub dropped dims')
      call assert_close_array(out%data, real([7, 8, 9, 10, 11, 12], dp), 'asub data')

      idx(1)%values = [2]
      idx(2)%values = [1, 2]
      call asub_array(x, idx, [3, 1], out, status, drop=.false.)
      call assert_equal_int(status, 0, 'asub reversed dimensions status')
      call assert_equal_int_array(out%dims, [2, 3, 1], 'asub reversed dimensions shape')
      call assert_close_array(out%data, real([7, 8, 9, 10, 11, 12], dp), 'asub reversed dimensions data')
   end subroutine test_asub

   subroutine test_adrop()
      type(array_value) :: x
      type(array_value) :: out
      integer :: status

      call init_array(x, real([10, 20, 30], dp), [1, 3, 1], status)
      call assert_equal_int(status, 0, 'init adrop array')
      call adrop_array(x, [1, 3], out, status)
      call assert_equal_int(status, 0, 'adrop status')
      call assert_equal_int_array(out%dims, [3], 'adrop dims')
      call assert_close_array(out%data, real([10, 20, 30], dp), 'adrop data')
   end subroutine test_adrop

   subroutine test_acorn()
      type(array_value) :: x
      type(array_value) :: out
      integer :: status
      integer :: i

      call init_array(x, [(real(i, dp), i = 1, 24)], [4, 3, 2], status)
      call assert_equal_int(status, 0, 'init acorn array')
      call acorn_array(x, [-2, 2, 1], out, status, addrownums=.false.)
      call assert_equal_int(status, 0, 'acorn status')
      call assert_equal_int_array(out%dims, [2, 2, 1], 'acorn dims')
      call assert_close_array(out%data, real([3, 4, 7, 8], dp), 'acorn data')
   end subroutine test_acorn

   subroutine test_afill()
      type(array_value) :: x
      type(array_value) :: value
      type(int_vector) :: selectors(3)
      integer :: status

      call init_array(x, [(0.0_dp, status = 1, 24)], [2, 3, 4], status)
      call assert_equal_int(status, 0, 'init afill destination')
      call set_dimnames(x, 2, ['A', 'B', 'C'], status)
      call assert_equal_int(status, 0, 'set afill dim 2 names')
      call set_dimnames(x, 3, ['w', 'x', 'y', 'z'], status)
      call assert_equal_int(status, 0, 'set afill dim 3 names')

      call init_array(value, real([1, 2, 3, 4, 5, 6], dp), [2, 3], status)
      call assert_equal_int(status, 0, 'init afill value')
      call set_dimnames(value, 1, ['A', 'B'], status)
      call assert_equal_int(status, 0, 'set afill value dim 1 names')
      call set_dimnames(value, 2, ['w', 'x', 'y'], status)
      call assert_equal_int(status, 0, 'set afill value dim 2 names')

      selectors(1)%values = [1, 2]
      allocate(selectors(2)%values(0), selectors(3)%values(0))
      call afill_array(x, selectors, value, status)
      call assert_equal_int(status, 0, 'afill status')
      call assert_close_array(x%data(1:6), real([1, 1, 2, 2, 0, 0], dp), 'afill first third-dimension block')
      call assert_close_array(x%data(7:12), real([3, 3, 4, 4, 0, 0], dp), 'afill second third-dimension block')
      call assert_close_array(x%data(13:18), real([5, 5, 6, 6, 0, 0], dp), 'afill third third-dimension block')
      call assert_close_array(x%data(19:24), real([0, 0, 0, 0, 0, 0], dp), 'afill untouched block')
   end subroutine test_afill

   subroutine assert_equal_int(actual, expected, label)
      integer, intent(in) :: actual !! Observed integer value.
      integer, intent(in) :: expected !! Required integer value.
      character(len=*), intent(in) :: label !! Test label printed when the assertion fails.

      if (actual /= expected) then
         print '(a,2(1x,i0))', trim(label), actual, expected
         error stop 1
      end if
   end subroutine assert_equal_int

   subroutine assert_equal_int_array(actual, expected, label)
      integer, intent(in) :: actual(:) !! Observed integer array.
      integer, intent(in) :: expected(:) !! Required integer array.
      character(len=*), intent(in) :: label !! Test label printed when the assertion fails.

      if (size(actual) /= size(expected) .or. any(actual /= expected)) then
         print '(a)', trim(label)
         print '(a,*(1x,i0))', 'actual:', actual
         print '(a,*(1x,i0))', 'expected:', expected
         error stop 1
      end if
   end subroutine assert_equal_int_array

   subroutine assert_close_array(actual, expected, label)
      real(dp), intent(in) :: actual(:) !! Observed floating-point values.
      real(dp), intent(in) :: expected(:) !! Required floating-point values.
      character(len=*), intent(in) :: label !! Test label printed when the assertion fails.

      if (size(actual) /= size(expected) .or. any(abs(actual - expected) > 1.0e-12_dp)) then
         print '(a)', trim(label)
         print '(a,*(1x,g0))', 'actual:', actual
         print '(a,*(1x,g0))', 'expected:', expected
         error stop 1
      end if
   end subroutine assert_close_array

end program test_unit
