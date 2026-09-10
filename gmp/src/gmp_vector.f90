module gmp_vector
   use, intrinsic :: iso_fortran_env, only : int64
   use gmp_bigz, only : bigz, bigz_from_int64, bigz_compare, bigz_equal, bigz_is_zero, bigz_abs, &
      bigz_add, bigz_mul, bigz_pow, bigz_fits_int64, bigz_to_int64
   use gmp_bigq, only : bigq, bigq_from_int64, bigq_add, bigq_mul, bigq_div
   use gmp_number_theory, only : factorial_z
   implicit none
   private

   public :: bigz_sum, bigz_prod, bigz_cumsum, bigz_min, bigz_max, bigz_gamma
   public :: bigq_sum, bigq_prod, bigq_cumsum, bigq_min, bigq_max, bigq_mean
   public :: bigz_unique, bigz_duplicated, bigz_diff, bigq_diff
   public :: bigq_unique, bigq_duplicated

contains

   pure function bigz_sum(values) result(total)
      type(bigz), intent(in) :: values(:) !! Integer vector whose exact sum is requested.
      type(bigz) :: total
      integer :: i

      total = bigz_from_int64(0_int64)
      do i = 1, size(values)
         total = bigz_add(total, values(i))
      end do
   end function bigz_sum

   pure function bigz_prod(values) result(product)
      type(bigz), intent(in) :: values(:) !! Integer vector whose exact product is requested.
      type(bigz) :: product
      integer :: i

      product = bigz_from_int64(1_int64)
      do i = 1, size(values)
         product = bigz_mul(product, values(i))
      end do
   end function bigz_prod

   pure function bigz_cumsum(values) result(result)
      type(bigz), intent(in) :: values(:) !! Integer vector for cumulative exact sums.
      type(bigz), allocatable :: result(:)
      type(bigz) :: running
      integer :: i

      allocate(result(size(values)))
      running = bigz_from_int64(0_int64)
      do i = 1, size(values)
         running = bigz_add(running, values(i))
         result(i) = running
      end do
   end function bigz_cumsum

   pure function bigz_min(values) result(value)
      type(bigz), intent(in) :: values(:) !! Nonempty integer vector whose minimum is requested.
      type(bigz) :: value
      integer :: i

      if (size(values) == 0) then
         value = bigz_from_int64(0_int64)
         return
      end if
      value = values(1)
      do i = 2, size(values)
         if (bigz_compare(values(i), value) < 0) value = values(i)
      end do
   end function bigz_min

   pure function bigz_max(values) result(value)
      type(bigz), intent(in) :: values(:) !! Nonempty integer vector whose maximum is requested.
      type(bigz) :: value
      integer :: i

      if (size(values) == 0) then
         value = bigz_from_int64(0_int64)
         return
      end if
      value = values(1)
      do i = 2, size(values)
         if (bigz_compare(values(i), value) > 0) value = values(i)
      end do
   end function bigz_max

   pure function bigz_gamma(x) result(value)
      type(bigz), intent(in) :: x !! Positive integer argument for gamma(x)=(x-1)!; unsupported huge indices return zero.
      type(bigz) :: value
      integer(int64) :: n

      if (.not. bigz_fits_int64(x)) then
         value = bigz_from_int64(0_int64)
         return
      end if
      n = bigz_to_int64(x)
      if (n <= 0_int64) then
         value = bigz_from_int64(0_int64)
      else
         value = factorial_z(n - 1_int64)
      end if
   end function bigz_gamma

   pure function bigq_sum(values) result(total)
      type(bigq), intent(in) :: values(:) !! Rational vector whose exact sum is requested.
      type(bigq) :: total
      integer :: i

      total = bigq_from_int64(0_int64)
      do i = 1, size(values)
         total = bigq_add(total, values(i))
      end do
   end function bigq_sum

   pure function bigq_prod(values) result(product)
      type(bigq), intent(in) :: values(:) !! Rational vector whose exact product is requested.
      type(bigq) :: product
      integer :: i

      product = bigq_from_int64(1_int64)
      do i = 1, size(values)
         product = bigq_mul(product, values(i))
      end do
   end function bigq_prod

   pure function bigq_cumsum(values) result(result)
      type(bigq), intent(in) :: values(:) !! Rational vector for cumulative exact sums.
      type(bigq), allocatable :: result(:)
      type(bigq) :: running
      integer :: i

      allocate(result(size(values)))
      running = bigq_from_int64(0_int64)
      do i = 1, size(values)
         running = bigq_add(running, values(i))
         result(i) = running
      end do
   end function bigq_cumsum

   pure function bigq_min(values) result(value)
      use gmp_bigq, only : bigq_compare
      type(bigq), intent(in) :: values(:) !! Nonempty rational vector whose minimum is requested.
      type(bigq) :: value
      integer :: i

      if (size(values) == 0) then
         value = bigq_from_int64(0_int64)
         return
      end if
      value = values(1)
      do i = 2, size(values)
         if (bigq_compare(values(i), value) < 0) value = values(i)
      end do
   end function bigq_min

   pure function bigq_max(values) result(value)
      use gmp_bigq, only : bigq_compare
      type(bigq), intent(in) :: values(:) !! Nonempty rational vector whose maximum is requested.
      type(bigq) :: value
      integer :: i

      if (size(values) == 0) then
         value = bigq_from_int64(0_int64)
         return
      end if
      value = values(1)
      do i = 2, size(values)
         if (bigq_compare(values(i), value) > 0) value = values(i)
      end do
   end function bigq_max

   pure function bigq_mean(values) result(value)
      type(bigq), intent(in) :: values(:) !! Nonempty rational vector whose exact arithmetic mean is requested.
      type(bigq) :: value

      if (size(values) == 0) then
         value = bigq_from_int64(0_int64)
      else
         value = bigq_div(bigq_sum(values), bigq_from_int64(int(size(values), int64)))
      end if
   end function bigq_mean

   pure function bigz_duplicated(values) result(is_duplicate)
      type(bigz), intent(in) :: values(:) !! Integer vector tested for earlier equal values.
      logical, allocatable :: is_duplicate(:)
      integer :: i, j

      allocate(is_duplicate(size(values)))
      is_duplicate = .false.
      do i = 2, size(values)
         do j = 1, i - 1
            if (bigz_equal(values(i), values(j))) then
               is_duplicate(i) = .true.
               exit
            end if
         end do
      end do
   end function bigz_duplicated

   pure function bigz_unique(values) result(unique_values)
      type(bigz), intent(in) :: values(:) !! Integer vector reduced to first occurrences, preserving order.
      type(bigz), allocatable :: unique_values(:)
      logical, allocatable :: duplicated(:)
      integer :: i, n

      duplicated = bigz_duplicated(values)
      n = count(.not. duplicated)
      allocate(unique_values(n))
      n = 0
      do i = 1, size(values)
         if (.not. duplicated(i)) then
            n = n + 1
            unique_values(n) = values(i)
         end if
      end do
   end function bigz_unique

   pure function bigz_diff(values, lag, differences) result(result)
      type(bigz), intent(in) :: values(:) !! Integer vector to difference repeatedly.
      integer, intent(in) :: lag !! Positive lag between subtracted elements.
      integer, intent(in) :: differences !! Positive number of differencing passes.
      type(bigz), allocatable :: result(:)
      type(bigz), allocatable :: current(:), next_values(:)
      integer :: pass, i, n

      if (lag < 1 .or. differences < 1) then
         allocate(result(0))
         return
      end if
      current = values
      do pass = 1, differences
         n = size(current) - lag
         if (n <= 0) then
            allocate(result(0))
            return
         end if
         allocate(next_values(n))
         do i = 1, n
            next_values(i) = bigz_add(current(i + lag), bigz_mul(bigz_from_int64(-1_int64), current(i)))
         end do
         call move_alloc(next_values, current)
      end do
      result = current
   end function bigz_diff

   pure function bigq_diff(values, lag, differences) result(result)
      use gmp_bigq, only : bigq_sub
      type(bigq), intent(in) :: values(:) !! Rational vector to difference repeatedly.
      integer, intent(in) :: lag !! Positive lag between subtracted elements.
      integer, intent(in) :: differences !! Positive number of differencing passes.
      type(bigq), allocatable :: result(:)
      type(bigq), allocatable :: current(:), next_values(:)
      integer :: pass, i, n

      if (lag < 1 .or. differences < 1) then
         allocate(result(0))
         return
      end if
      current = values
      do pass = 1, differences
         n = size(current) - lag
         if (n <= 0) then
            allocate(result(0))
            return
         end if
         allocate(next_values(n))
         do i = 1, n
            next_values(i) = bigq_sub(current(i + lag), current(i))
         end do
         call move_alloc(next_values, current)
      end do
      result = current
   end function bigq_diff

   pure function bigq_duplicated(values) result(is_duplicate)
      use gmp_bigq, only : bigq_equal
      type(bigq), intent(in) :: values(:) !! Rational vector tested for earlier equal values.
      logical, allocatable :: is_duplicate(:)
      integer :: i, j

      allocate(is_duplicate(size(values)))
      is_duplicate = .false.
      do i = 2, size(values)
         do j = 1, i - 1
            if (bigq_equal(values(i), values(j))) then
               is_duplicate(i) = .true.
               exit
            end if
         end do
      end do
   end function bigq_duplicated

   pure function bigq_unique(values) result(unique_values)
      type(bigq), intent(in) :: values(:) !! Rational vector reduced to first occurrences, preserving order.
      type(bigq), allocatable :: unique_values(:)
      logical, allocatable :: duplicated(:)
      integer :: i, n

      duplicated = bigq_duplicated(values)
      n = count(.not. duplicated)
      allocate(unique_values(n))
      n = 0
      do i = 1, size(values)
         if (.not. duplicated(i)) then
            n = n + 1
            unique_values(n) = values(i)
         end if
      end do
   end function bigq_unique

end module gmp_vector
