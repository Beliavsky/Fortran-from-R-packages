! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module sandwich_utils
   use sandwich_kinds, only : dp
   implicit none
   private

   public :: lowercase, count_bits, group_subset, group_vector
   public :: aggregate_rows, sort_index_integer, unique_values_integer

contains

   pure function lowercase(text) result(lower)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: i, code

      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            lower(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
   end function lowercase

   pure integer function count_bits(mask, nbits) result(nset)
      integer, intent(in) :: mask, nbits
      integer :: j

      nset = 0
      do j = 0, nbits - 1
         if (btest(mask, j)) nset = nset + 1
      end do
   end function count_bits

   subroutine group_vector(values, labels, ngroups)
      integer, intent(in) :: values(:)
      integer, allocatable, intent(out) :: labels(:)
      integer, intent(out) :: ngroups
      integer, allocatable :: representatives(:)
      integer :: n, i, g
      logical :: found

      n = size(values)
      allocate(labels(n), representatives(max(n, 1)))
      ngroups = 0
      do i = 1, n
         found = .false.
         do g = 1, ngroups
            if (values(i) == representatives(g)) then
               labels(i) = g
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            ngroups = ngroups + 1
            representatives(ngroups) = values(i)
            labels(i) = ngroups
         end if
      end do
   end subroutine group_vector

   subroutine group_subset(cluster, mask, labels, ngroups)
      integer, intent(in) :: cluster(:, :)
      integer, intent(in) :: mask
      integer, allocatable, intent(out) :: labels(:)
      integer, intent(out) :: ngroups
      integer, allocatable :: representative_rows(:)
      integer :: n, p, i, g, j
      logical :: found, equal_group

      n = size(cluster, 1)
      p = size(cluster, 2)
      allocate(labels(n), representative_rows(max(n, 1)))
      ngroups = 0
      do i = 1, n
         found = .false.
         do g = 1, ngroups
            equal_group = .true.
            do j = 1, p
               if (btest(mask, j - 1)) then
                  if (cluster(i, j) /= cluster(representative_rows(g), j)) then
                     equal_group = .false.
                     exit
                  end if
               end if
            end do
            if (equal_group) then
               labels(i) = g
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            ngroups = ngroups + 1
            representative_rows(ngroups) = i
            labels(i) = ngroups
         end if
      end do
   end subroutine group_subset

   subroutine aggregate_rows(x, labels, ngroups, sums)
      real(dp), intent(in) :: x(:, :)
      integer, intent(in) :: labels(:), ngroups
      real(dp), allocatable, intent(out) :: sums(:, :)
      integer :: i

      allocate(sums(ngroups, size(x, 2)))
      sums = 0.0_dp
      do i = 1, size(x, 1)
         sums(labels(i), :) = sums(labels(i), :) + x(i, :)
      end do
   end subroutine aggregate_rows

   subroutine sort_index_integer(values, index)
      integer, intent(in) :: values(:)
      integer, allocatable, intent(out) :: index(:)
      integer :: i, j, key

      allocate(index(size(values)))
      do i = 1, size(values)
         index(i) = i
      end do
      do i = 2, size(values)
         key = index(i)
         j = i - 1
         do while (j >= 1)
            if (values(index(j)) <= values(key)) exit
            index(j + 1) = index(j)
            j = j - 1
         end do
         index(j + 1) = key
      end do
   end subroutine sort_index_integer

   subroutine unique_values_integer(values, unique_values, inverse)
      integer, intent(in) :: values(:)
      integer, allocatable, intent(out) :: unique_values(:), inverse(:)
      integer, allocatable :: work(:)
      integer :: i, j, n_unique
      logical :: found

      allocate(work(max(size(values), 1)), inverse(size(values)))
      n_unique = 0
      do i = 1, size(values)
         found = .false.
         do j = 1, n_unique
            if (values(i) == work(j)) then
               inverse(i) = j
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            n_unique = n_unique + 1
            work(n_unique) = values(i)
            inverse(i) = n_unique
         end if
      end do
      allocate(unique_values(n_unique))
      if (n_unique > 0) unique_values = work(1:n_unique)
   end subroutine unique_values_integer

end module sandwich_utils
