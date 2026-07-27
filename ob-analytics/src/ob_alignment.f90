! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of obAnalytics.
! Copyright (C) 2015,2016 Philip Stubbings.
module ob_alignment
   use ob_kinds, only : dp, i8
   implicit none
   private
   public :: similarity_matrix_equal, similarity_matrix_time
   public :: needleman_wunsch

contains

   pure function similarity_matrix_equal(a, b) result(score)
      real(dp), intent(in) :: a(:), b(:)
      real(dp), allocatable :: score(:,:)
      integer :: i, j
      allocate(score(size(a), size(b)))
      do j = 1, size(b)
         do i = 1, size(a)
            if (a(i) == b(j)) then
               score(i,j) = 1.0_dp
            else
               score(i,j) = -1.0_dp
            end if
         end do
      end do
   end function similarity_matrix_equal

   pure function similarity_matrix_time(a_ms, b_ms, cutoff_ms) result(score)
      integer(i8), intent(in) :: a_ms(:), b_ms(:)
      real(dp), intent(in) :: cutoff_ms
      real(dp), allocatable :: score(:,:)
      integer :: i, j
      real(dp) :: distance
      allocate(score(size(a_ms), size(b_ms)))
      do j = 1, size(b_ms)
         do i = 1, size(a_ms)
            distance = real(abs(a_ms(i) - b_ms(j)), dp)
            if (distance == 0.0_dp) then
               score(i,j) = cutoff_ms
            else
               score(i,j) = cutoff_ms/distance
            end if
         end do
      end do
   end function similarity_matrix_time

   function needleman_wunsch(score, gap) result(alignment)
      real(dp), intent(in) :: score(:,:)
      real(dp), intent(in), optional :: gap
      integer, allocatable :: alignment(:,:)
      real(dp), allocatable :: f(:,:)
      integer, allocatable :: reverse_pairs(:,:)
      real(dp) :: gap_value, diagonal, up, left, tolerance
      integer :: i, j, n, m, count

      gap_value = -1.0_dp
      if (present(gap)) gap_value = gap
      n = size(score,1)
      m = size(score,2)
      allocate(f(0:n,0:m))
      f = 0.0_dp
      do i = 0, n
         f(i,0) = real(i,dp)*gap_value
      end do
      do j = 0, m
         f(0,j) = real(j,dp)*gap_value
      end do
      do i = 1, n
         do j = 1, m
            f(i,j) = max(f(i-1,j-1) + score(i,j), f(i-1,j) + gap_value, &
               f(i,j-1) + gap_value)
         end do
      end do

      allocate(reverse_pairs(2, min(n,m)))
      count = 0
      i = n
      j = m
      do while (i > 0 .or. j > 0)
         tolerance = 32.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(f(i,j)))
         if (i > 0 .and. j > 0) then
            diagonal = f(i-1,j-1) + score(i,j)
         else
            diagonal = -huge(1.0_dp)
         end if
         if (i > 0) then
            up = f(i-1,j) + gap_value
         else
            up = -huge(1.0_dp)
         end if
         if (j > 0) then
            left = f(i,j-1) + gap_value
         else
            left = -huge(1.0_dp)
         end if

         if (i > 0 .and. j > 0 .and. abs(f(i,j)-diagonal) <= tolerance) then
            count = count + 1
            reverse_pairs(:,count) = [i,j]
            i = i - 1
            j = j - 1
         else if (i > 0 .and. abs(f(i,j)-up) <= tolerance) then
            i = i - 1
         else if (j > 0 .and. abs(f(i,j)-left) <= tolerance) then
            j = j - 1
         else
            error stop 'needleman_wunsch: inconsistent backtrace'
         end if
      end do

      allocate(alignment(count,2))
      do i = 1, count
         alignment(i,:) = reverse_pairs(:,count-i+1)
      end do
   end function needleman_wunsch

end module ob_alignment
