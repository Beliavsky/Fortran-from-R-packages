! rmutil computational translation
! Copyright (C) 1998-2001 J.K. Lindsey
! Copyright (C) 2026 OpenAI (modern Fortran translation)
! SPDX-License-Identifier: GPL-2.0-or-later
module rmutil_utils
   use rmutil_kinds, only : dp
   implicit none
   private
   public :: contrast_mean, group_sum, group_mean
contains
   function contrast_mean(n, contrasts) result(c)
      integer, intent(in) :: n
      logical, intent(in), optional :: contrasts
      real(dp), allocatable :: c(:,:)
      logical :: use_contrasts
      integer :: i
      if (n <= 1) error stop "contrast_mean: n must exceed 1"
      use_contrasts = .true.
      if (present(contrasts)) use_contrasts = contrasts
      if (use_contrasts) then
         allocate(c(n,n-1)); c = 0.0_dp
         do i = 1, n-1
            c(i,i) = 1.0_dp
            c(n,i) = -1.0_dp
         end do
      else
         allocate(c(n,n)); c = 0.0_dp
         do i = 1, n
            c(i,i) = 1.0_dp
         end do
      end if
   end function contrast_mean

   function group_sum(x,index,ngroup) result(ans)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: index(:)
      integer, intent(in), optional :: ngroup
      real(dp), allocatable :: ans(:)
      integer :: n, i
      if (size(x) /= size(index)) error stop "group_sum: x/index size mismatch"
      n = maxval(index)
      if (present(ngroup)) n = ngroup
      allocate(ans(n)); ans = 0.0_dp
      do i = 1, size(x)
         if (index(i) >= 1 .and. index(i) <= n) ans(index(i)) = ans(index(i)) + x(i)
      end do
   end function group_sum

   function group_mean(x,index,ngroup) result(ans)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: index(:)
      integer, intent(in), optional :: ngroup
      real(dp), allocatable :: ans(:)
      integer, allocatable :: count(:)
      integer :: n, i
      if (size(x) /= size(index)) error stop "group_mean: x/index size mismatch"
      n = maxval(index)
      if (present(ngroup)) n = ngroup
      allocate(ans(n),count(n)); ans = 0.0_dp; count = 0
      do i = 1, size(x)
         if (index(i) >= 1 .and. index(i) <= n) then
            ans(index(i)) = ans(index(i)) + x(i)
            count(index(i)) = count(index(i)) + 1
         end if
      end do
      do i = 1, n
         if (count(i) > 0) ans(i) = ans(i)/real(count(i),dp)
      end do
   end function group_mean
end module rmutil_utils
