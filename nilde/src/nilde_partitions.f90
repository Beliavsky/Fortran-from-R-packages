! SPDX-License-Identifier: GPL-2.0-or-later
module nilde_partitions
   use nilde_kinds, only : i8
   use nilde_types, only : integer_solutions_t
   use nilde_collect, only : append_i8_column
   implicit none
   private
   public :: get_partitions
contains

   function get_partitions(n, m, at_most) result(res)
      integer(i8), intent(in) :: n
      integer, intent(in) :: m
      logical, intent(in), optional :: at_most
      type(integer_solutions_t) :: res
      integer(i8), allocatable :: x(:), store(:,:)
      logical :: am
      integer(i8) :: first

      if (n <= 0_i8 .or. m <= 0 .or. int(m,i8) > n) error stop "get_partitions: require 1 <= M <= n"
      am = .true.; if (present(at_most)) am = at_most
      allocate(x(m)); x = 0_i8
      first = merge(0_i8, 1_i8, am)
      res%nsol = 0
      call rec(1, n, first)
      if (res%nsol > 0) then
         allocate(res%x(m,res%nsol)); res%x = store(:,1:res%nsol)
      end if

   contains
      recursive subroutine rec(pos, rem, lo)
         integer, intent(in) :: pos
         integer(i8), intent(in) :: rem, lo
         integer(i8) :: v, vmax
         integer :: left
         left = m-pos+1
         if (pos == m) then
            if (rem < lo) return
            x(pos) = rem
            call append_i8_column(store, res%nsol, x)
            return
         end if
         if (rem < int(left,i8)*lo) return
         vmax = rem/int(left,i8)
         do v = lo, vmax
            x(pos)=v
            call rec(pos+1, rem-v, v)
         end do
      end subroutine rec
   end function get_partitions

end module nilde_partitions
