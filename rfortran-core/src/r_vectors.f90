! SPDX-License-Identifier: MIT
module r_vectors
   use r_kinds, only : dp
   use r_status, only : r_invalid_input, r_ok
   implicit none
   private

   public :: r_difference

   interface r_difference
      module procedure r_difference_vector
      module procedure r_difference_matrix
   end interface r_difference

contains

   pure subroutine r_difference_vector(x, values, lag, differences, status)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(in), optional :: lag, differences
      integer, intent(out), optional :: status
      real(dp), allocatable :: work(:), next(:)
      integer :: d, k, nlag

      nlag = 1
      if (present(lag)) nlag = lag
      d = 1
      if (present(differences)) d = differences
      if (present(status)) status = r_ok
      if (nlag < 1 .or. d < 1 .or. size(x) <= nlag*d) then
         if (present(status)) status = r_invalid_input
         return
      end if

      allocate(work(size(x)))
      work = x
      do k = 1, d
         allocate(next(size(work) - nlag))
         next = work(1+nlag:) - work(:size(work)-nlag)
         call move_alloc(next, work)
      end do
      call move_alloc(work, values)
   end subroutine r_difference_vector

   pure subroutine r_difference_matrix(x, values, lag, differences, status)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: values(:,:)
      integer, intent(in), optional :: lag, differences
      integer, intent(out), optional :: status
      real(dp), allocatable :: work(:,:), next(:,:)
      integer :: d, k, nlag

      nlag = 1
      if (present(lag)) nlag = lag
      d = 1
      if (present(differences)) d = differences
      if (present(status)) status = r_ok
      if (nlag < 1 .or. d < 1 .or. size(x, 1) <= nlag*d .or. size(x, 2) < 1) then
         if (present(status)) status = r_invalid_input
         return
      end if

      allocate(work(size(x, 1), size(x, 2)))
      work = x
      do k = 1, d
         allocate(next(size(work, 1) - nlag, size(work, 2)))
         next = work(1+nlag:, :) - work(:size(work, 1)-nlag, :)
         call move_alloc(next, work)
      end do
      call move_alloc(work, values)
   end subroutine r_difference_matrix

end module r_vectors
