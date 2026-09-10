! SPDX-License-Identifier: MIT
! SPDX-FileComment: Typed functional operations for the Fortran purrr translation.
module purrr_core
   !! Implements mapping, folding, filtering, detection, and side-effect walks.
   use purrr_callbacks
   implicit none
   private

   public :: accumulate, detect, detect_index, discard, every, head_while, keep
   public :: map, map2, map2_dbl, map2_int, map_dbl, map_int, modify, negate
   public :: none, reduce, some, tail_while, walk

   interface map
      module procedure map_dbl, map_int
   end interface map
   interface map2
      module procedure map2_dbl, map2_int
   end interface map2
   interface modify
      module procedure map_dbl, map_int
   end interface modify
   interface reduce
      module procedure reduce_real, reduce_integer
   end interface reduce
   interface accumulate
      module procedure accumulate_real, accumulate_integer
   end interface accumulate
   interface keep
      module procedure keep_real, keep_integer
   end interface keep
   interface discard
      module procedure discard_real, discard_integer
   end interface discard
   interface some
      module procedure some_real, some_integer
   end interface some
   interface every
      module procedure every_real, every_integer
   end interface every
   interface none
      module procedure none_real, none_integer
   end interface none
   interface detect
      module procedure detect_real, detect_integer
   end interface detect
   interface detect_index
      module procedure detect_index_real, detect_index_integer
   end interface detect_index
   interface head_while
      module procedure head_while_real, head_while_integer
   end interface head_while
   interface tail_while
      module procedure tail_while_real, tail_while_integer
   end interface tail_while
   interface walk
      module procedure walk_real, walk_integer
   end interface walk

contains

   pure function map_dbl(values, callback) result(output)
      !! Applies a pure real callback element by element.
      real(dp), intent(in) :: values(:) !! Input values.
      procedure(real_unary) :: callback !! Scalar transformation.
      real(dp), allocatable :: output(:)
      integer :: i
      allocate (output(size(values)))
      do i = 1, size(values)
         output(i) = callback(values(i))
      end do
   end function map_dbl

   pure function map_int(values, callback) result(output)
      !! Applies a pure integer callback element by element.
      integer, intent(in) :: values(:) !! Input values.
      procedure(integer_unary) :: callback !! Scalar transformation.
      integer, allocatable :: output(:)
      integer :: i
      allocate (output(size(values)))
      do i = 1, size(values)
         output(i) = callback(values(i))
      end do
   end function map_int

   pure function map2_dbl(left, right, callback) result(output)
      !! Maps a pure callback over two conformable real arrays.
      real(dp), intent(in) :: left(:)  !! First inputs.
      real(dp), intent(in) :: right(:) !! Second inputs.
      procedure(real_binary) :: callback !! Binary transformation.
      real(dp), allocatable :: output(:)
      integer :: i
      if (size(left) /= size(right)) error stop 'map2_dbl: array sizes differ'
      allocate (output(size(left)))
      do i = 1, size(left)
         output(i) = callback(left(i), right(i))
      end do
   end function map2_dbl

   pure function map2_int(left, right, callback) result(output)
      !! Maps a pure callback over two conformable integer arrays.
      integer, intent(in) :: left(:)  !! First inputs.
      integer, intent(in) :: right(:) !! Second inputs.
      procedure(integer_binary) :: callback !! Binary transformation.
      integer, allocatable :: output(:)
      integer :: i
      if (size(left) /= size(right)) error stop 'map2_int: array sizes differ'
      allocate (output(size(left)))
      do i = 1, size(left)
         output(i) = callback(left(i), right(i))
      end do
   end function map2_int

   pure real(dp) function reduce_real(values, callback, initial) result(output)
      !! Left-folds a real array.
      real(dp), intent(in) :: values(:) !! Input values.
      procedure(real_binary) :: callback !! Binary accumulator.
      real(dp), intent(in), optional :: initial !! Optional initial value.
      integer :: i
      if (present(initial)) then
         output = initial
         i = 1
      else
         if (size(values) == 0) error stop 'reduce: empty input requires initial'
         output = values(1)
         i = 2
      end if
      do while (i <= size(values))
         output = callback(output, values(i))
         i = i + 1
      end do
   end function reduce_real

   pure integer function reduce_integer(values, callback, initial) result(output)
      !! Left-folds an integer array.
      integer, intent(in) :: values(:) !! Input values.
      procedure(integer_binary) :: callback !! Binary accumulator.
      integer, intent(in), optional :: initial !! Optional initial value.
      integer :: i
      if (present(initial)) then
         output = initial
         i = 1
      else
         if (size(values) == 0) error stop 'reduce: empty input requires initial'
         output = values(1)
         i = 2
      end if
      do while (i <= size(values))
         output = callback(output, values(i))
         i = i + 1
      end do
   end function reduce_integer

   pure function accumulate_real(values, callback, initial) result(output)
      !! Returns every intermediate value of a real left fold.
      real(dp), intent(in) :: values(:) !! Input values.
      procedure(real_binary) :: callback !! Binary accumulator.
      real(dp), intent(in), optional :: initial !! Optional initial value.
      real(dp), allocatable :: output(:)
      integer :: i, offset
      offset = merge(1, 0, present(initial))
      allocate (output(size(values) + offset))
      if (present(initial)) then
         output(1) = initial
      else if (size(values) > 0) then
         output(1) = values(1)
      end if
      if (present(initial)) then
         do i = 1, size(values)
            output(i + 1) = callback(output(i), values(i))
         end do
      else
         do i = 2, size(values)
            output(i) = callback(output(i - 1), values(i))
         end do
      end if
   end function accumulate_real

   pure function accumulate_integer(values, callback, initial) result(output)
      !! Returns every intermediate value of an integer left fold.
      integer, intent(in) :: values(:) !! Input values.
      procedure(integer_binary) :: callback !! Binary accumulator.
      integer, intent(in), optional :: initial !! Optional initial value.
      integer, allocatable :: output(:)
      integer :: i, offset
      offset = merge(1, 0, present(initial))
      allocate (output(size(values) + offset))
      if (present(initial)) then
         output(1) = initial
      else if (size(values) > 0) then
         output(1) = values(1)
      end if
      if (present(initial)) then
         do i = 1, size(values)
            output(i + 1) = callback(output(i), values(i))
         end do
      else
         do i = 2, size(values)
            output(i) = callback(output(i - 1), values(i))
         end do
      end if
   end function accumulate_integer

   pure function keep_real(values, predicate) result(output)
      !! Keeps real values satisfying a predicate.
      real(dp), intent(in) :: values(:) !! Input values.
      procedure(real_predicate) :: predicate !! Selection predicate.
      real(dp), allocatable :: output(:)
      logical, allocatable :: selected(:)
      integer :: i
      allocate (selected(size(values)))
      do i = 1, size(values)
         selected(i) = predicate(values(i))
      end do
      output = pack(values, selected)
   end function keep_real

   pure function discard_real(values, predicate) result(output)
      !! Discards real values satisfying a predicate.
      real(dp), intent(in) :: values(:) !! Input values.
      procedure(real_predicate) :: predicate !! Rejection predicate.
      real(dp), allocatable :: output(:)
      logical, allocatable :: selected(:)
      integer :: i
      allocate (selected(size(values)))
      do i = 1, size(values)
         selected(i) = .not. predicate(values(i))
      end do
      output = pack(values, selected)
   end function discard_real

   pure function keep_integer(values, predicate) result(output)
      !! Keeps integer values satisfying a predicate.
      integer, intent(in) :: values(:) !! Input values.
      procedure(integer_predicate) :: predicate !! Selection predicate.
      integer, allocatable :: output(:)
      logical, allocatable :: selected(:)
      integer :: i
      allocate (selected(size(values)))
      do i = 1, size(values)
         selected(i) = predicate(values(i))
      end do
      output = pack(values, selected)
   end function keep_integer

   pure function discard_integer(values, predicate) result(output)
      !! Discards integer values satisfying a predicate.
      integer, intent(in) :: values(:) !! Input values.
      procedure(integer_predicate) :: predicate !! Rejection predicate.
      integer, allocatable :: output(:)
      logical, allocatable :: selected(:)
      integer :: i
      allocate (selected(size(values)))
      do i = 1, size(values)
         selected(i) = .not. predicate(values(i))
      end do
      output = pack(values, selected)
   end function discard_integer

   pure logical function some_real(values, predicate) result(answer)
      !! Tests whether at least one real value satisfies a predicate.
      real(dp), intent(in) :: values(:) !! Values to test.
      procedure(real_predicate) :: predicate !! Predicate.
      answer = detect_index_real(values, predicate) > 0
   end function some_real

   pure logical function every_real(values, predicate) result(answer)
      !! Tests whether every real value satisfies a predicate.
      real(dp), intent(in) :: values(:) !! Values to test.
      procedure(real_predicate) :: predicate !! Predicate.
      answer = .true.
      block
         integer :: i
         do i = 1, size(values)
            if (.not. predicate(values(i))) answer = .false.
         end do
      end block
   end function every_real

   pure logical function none_real(values, predicate) result(answer)
      !! Tests whether no real value satisfies a predicate.
      real(dp), intent(in) :: values(:) !! Values to test.
      procedure(real_predicate) :: predicate !! Predicate.
      answer = detect_index_real(values, predicate) == 0
   end function none_real

   pure logical function some_integer(values, predicate) result(answer)
      !! Tests whether at least one integer satisfies a predicate.
      integer, intent(in) :: values(:) !! Values to test.
      procedure(integer_predicate) :: predicate !! Predicate.
      answer = detect_index_integer(values, predicate) > 0
   end function some_integer

   pure logical function every_integer(values, predicate) result(answer)
      !! Tests whether every integer satisfies a predicate.
      integer, intent(in) :: values(:) !! Values to test.
      procedure(integer_predicate) :: predicate !! Predicate.
      answer = .true.
      block
         integer :: i
         do i = 1, size(values)
            if (.not. predicate(values(i))) answer = .false.
         end do
      end block
   end function every_integer

   pure logical function none_integer(values, predicate) result(answer)
      !! Tests whether no integer satisfies a predicate.
      integer, intent(in) :: values(:) !! Values to test.
      procedure(integer_predicate) :: predicate !! Predicate.
      answer = detect_index_integer(values, predicate) == 0
   end function none_integer

   pure real(dp) function detect_real(values, predicate) result(value)
      !! Returns the first matching real value, or zero when absent.
      real(dp), intent(in) :: values(:) !! Values to search.
      procedure(real_predicate) :: predicate !! Match predicate.
      integer :: i
      value = 0.0_dp
      do i = 1, size(values)
         if (predicate(values(i))) then
            value = values(i)
            return
         end if
      end do
   end function detect_real

   pure integer function detect_integer(values, predicate) result(value)
      !! Returns the first matching integer, or zero when absent.
      integer, intent(in) :: values(:) !! Values to search.
      procedure(integer_predicate) :: predicate !! Match predicate.
      integer :: i
      value = 0
      do i = 1, size(values)
         if (predicate(values(i))) then
            value = values(i)
            return
         end if
      end do
   end function detect_integer

   pure integer function detect_index_real(values, predicate) result(index)
      !! Returns the first matching real position, or zero when absent.
      real(dp), intent(in) :: values(:) !! Values to search.
      procedure(real_predicate) :: predicate !! Match predicate.
      integer :: i
      index = 0
      do i = 1, size(values)
         if (predicate(values(i))) then
            index = i
            return
         end if
      end do
   end function detect_index_real

   pure integer function detect_index_integer(values, predicate) result(index)
      !! Returns the first matching integer position, or zero when absent.
      integer, intent(in) :: values(:) !! Values to search.
      procedure(integer_predicate) :: predicate !! Match predicate.
      integer :: i
      index = 0
      do i = 1, size(values)
         if (predicate(values(i))) then
            index = i
            return
         end if
      end do
   end function detect_index_integer

   pure function head_while_real(values, predicate) result(output)
      !! Retains leading real values while a predicate is true.
      real(dp), intent(in) :: values(:) !! Input values.
      procedure(real_predicate) :: predicate !! Continuation predicate.
      real(dp), allocatable :: output(:)
      integer :: n
      n = 0
      do while (n < size(values))
         if (.not. predicate(values(n + 1))) exit
         n = n + 1
      end do
      output = values(:n)
   end function head_while_real

   pure function head_while_integer(values, predicate) result(output)
      !! Retains leading integers while a predicate is true.
      integer, intent(in) :: values(:) !! Input values.
      procedure(integer_predicate) :: predicate !! Continuation predicate.
      integer, allocatable :: output(:)
      integer :: n
      n = 0
      do while (n < size(values))
         if (.not. predicate(values(n + 1))) exit
         n = n + 1
      end do
      output = values(:n)
   end function head_while_integer

   pure function tail_while_real(values, predicate) result(output)
      !! Retains trailing real values while a predicate is true.
      real(dp), intent(in) :: values(:) !! Input values.
      procedure(real_predicate) :: predicate !! Continuation predicate.
      real(dp), allocatable :: output(:)
      integer :: first
      first = size(values) + 1
      do while (first > 1)
         if (.not. predicate(values(first - 1))) exit
         first = first - 1
      end do
      output = values(first:)
   end function tail_while_real

   pure function tail_while_integer(values, predicate) result(output)
      !! Retains trailing integers while a predicate is true.
      integer, intent(in) :: values(:) !! Input values.
      procedure(integer_predicate) :: predicate !! Continuation predicate.
      integer, allocatable :: output(:)
      integer :: first
      first = size(values) + 1
      do while (first > 1)
         if (.not. predicate(values(first - 1))) exit
         first = first - 1
      end do
      output = values(first:)
   end function tail_while_integer

   pure function negate(predicate, value) result(answer)
      !! Applies the logical negation of a real predicate.
      procedure(real_predicate) :: predicate !! Predicate to negate.
      real(dp), intent(in) :: value !! Value to test.
      logical :: answer
      answer = .not. predicate(value)
   end function negate

   subroutine walk_real(values, callback)
      !! Calls an impure procedure for each real value.
      real(dp), intent(in) :: values(:) !! Input values.
      procedure(real_walker) :: callback !! Side-effecting callback.
      integer :: i
      do i = 1, size(values)
         call callback(values(i))
      end do
   end subroutine walk_real

   subroutine walk_integer(values, callback)
      !! Calls an impure procedure for each integer value.
      integer, intent(in) :: values(:) !! Input values.
      procedure(integer_walker) :: callback !! Side-effecting callback.
      integer :: i
      do i = 1, size(values)
         call callback(values(i))
      end do
   end subroutine walk_integer

end module purrr_core
