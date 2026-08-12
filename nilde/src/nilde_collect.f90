! SPDX-License-Identifier: GPL-2.0-or-later
module nilde_collect
   use nilde_kinds, only : i8
   implicit none
   private
   public :: append_i8_column, append_int_column
contains

   subroutine append_i8_column(a, ncol, v)
      integer(i8), allocatable, intent(inout) :: a(:,:)
      integer, intent(inout) :: ncol
      integer(i8), intent(in) :: v(:)
      integer(i8), allocatable :: tmp(:,:)
      integer :: cap, newcap

      if (.not. allocated(a)) then
         allocate(a(size(v), 16))
         a = 0_i8
      else if (size(a,1) /= size(v)) then
         error stop "append_i8_column: incompatible row count"
      end if
      cap = size(a,2)
      if (ncol >= cap) then
         newcap = max(2*cap, ncol+1)
         allocate(tmp(size(v), newcap))
         tmp = 0_i8
         if (ncol > 0) tmp(:,1:ncol) = a(:,1:ncol)
         call move_alloc(tmp, a)
      end if
      ncol = ncol + 1
      a(:,ncol) = v
   end subroutine append_i8_column

   subroutine append_int_column(a, ncol, v)
      integer, allocatable, intent(inout) :: a(:,:)
      integer, intent(inout) :: ncol
      integer, intent(in) :: v(:)
      integer, allocatable :: tmp(:,:)
      integer :: cap, newcap

      if (.not. allocated(a)) then
         allocate(a(size(v), 16))
         a = 0
      else if (size(a,1) /= size(v)) then
         error stop "append_int_column: incompatible row count"
      end if
      cap = size(a,2)
      if (ncol >= cap) then
         newcap = max(2*cap, ncol+1)
         allocate(tmp(size(v), newcap))
         tmp = 0
         if (ncol > 0) tmp(:,1:ncol) = a(:,1:ncol)
         call move_alloc(tmp, a)
      end if
      ncol = ncol + 1
      a(:,ncol) = v
   end subroutine append_int_column

end module nilde_collect
